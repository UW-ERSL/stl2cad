"""Straight sides + tangent fillets for TO void outlines (closed) and outline notches (open)."""
import pickle, numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import splprep, splev
from shapely.geometry import LinearRing, LineString, Point, Polygon
import shapely
from scipy.spatial import cKDTree

def dense(C, closed=True, h=0.02):
    """Polyline resampled so that no step exceeds h (straight sides included)."""
    P=np.vstack([C,C[:1]]) if closed else np.asarray(C)
    out=[]
    for a,b in zip(P[:-1],P[1:]):
        n=max(1,int(np.ceil(np.linalg.norm(b-a)/h))); out.append(a+np.linspace(0,1,n,endpoint=False)[:,None]*(b-a))
    out.append(P[-1:]); return np.vstack(out)

def dist_to(C, X, closed=True):
    return cKDTree(dense(C,closed)).query(X)[0]

def inside_depth(O_poly, X):
    """Signed: >0 inside the outline, distance to its boundary."""
    X=np.asarray(X); d=shapely.distance(O_poly.exterior, shapely.points(X))
    return np.where(shapely.contains_xy(O_poly,X[:,0],X[:,1]),d,-d)

def denoise(ring, step=0.25, rms=0.05):
    p=np.asarray(ring.coords)[:-1]
    p=p[np.r_[True,np.linalg.norm(np.diff(p,axis=0),axis=1)>1e-6]]
    tck,_=splprep(p.T,s=len(p)*rms**2,per=1)
    n=int(LinearRing(p).length/step)
    return p, np.array(splev(np.linspace(0,1,n,endpoint=False),tck)).T

def lsq_line(P):
    c=P.mean(0); d=np.linalg.svd(P-c)[2][0]
    return c, (d if d@(P[-1]-P[0])>0 else -d)

def maxdev(P):
    a,b=P[0],P[-1]; L=np.linalg.norm(b-a)
    if L<1e-9: return np.inf
    nr=np.array([-(b-a)[1],(b-a)[0]])/L; return np.abs((P-a)@nr).max()

def runs_of(Q, tol, Lmin):
    """Greedy straight runs along an open point sequence."""
    runs=[]; i=0; n=len(Q)
    while i<n-1:
        j=i+1
        while j+1<n and maxdev(Q[i:j+2])<=tol: j+=1
        if np.linalg.norm(Q[j]-Q[i])>=Lmin: runs.append([i,j]); i=j
        else: i+=1
    return runs

AMIN=8.0; TOLM=0.4
R_MIN=1.0     # smallest corner radius (no sharp re-entrant corners in the material)
MIN_FLAT=1.0  # shortest straight piece between two fillets
TIP_LOSS=0.5  # ... reduced at acute corners so the fillet shortens the tip by at most this much
def merge(Q, runs, amin=None, tolm=None):
    amin=AMIN if amin is None else amin; tolm=TOLM if tolm is None else tolm
    changed=True
    while changed and len(runs)>2:
        changed=False
        for a in range(len(runs)-1):
            (i1,j1),(i2,j2)=runs[a],runs[a+1]
            d1=lsq_line(Q[i1:j1+1])[1]; d2=lsq_line(Q[i2:j2+1])[1]
            if np.degrees(np.arccos(np.clip(d1@d2,-1,1)))<amin:
                c,d=lsq_line(Q[i1:j2+1]); nr=np.array([-d[1],d[0]])
                if np.abs((Q[i1:j2+1]-c)@nr).max()<=tolm:
                    runs[a:a+2]=[[i1,j2]]; changed=True; break
    return runs

def side_lines(Q, runs):
    out=[]
    for i,j in runs:
        k=max(1,(j-i)//8); out.append(lsq_line(Q[i+k:j+1-k]))
    return out

def isect(l1, l2):
    (c1,d1),(c2,d2)=l1,l2; A=np.c_[d1,-d2]
    if abs(np.linalg.det(A))<1e-9: return None
    return c1+np.linalg.solve(A,c2-c1)[0]*d1

def rounded(verts, closed=True, h=0.05):
    """Dense polyline of a polygon whose vertices carry fillet radii: verts = [(x,y,r)]."""
    V=np.array([v[:2] for v in verts]); R=[v[2] for v in verts]; n=len(V); out=[]
    idx=range(n) if closed else range(n)
    for k in idx:
        if (not closed and k in (0,n-1)) or R[k]<=1e-6: out.append(V[k]); continue
        a,b,c=V[k-1],V[k],V[(k+1)%n]
        u=(b-a)/np.linalg.norm(b-a); w=(c-b)/np.linalg.norm(c-b)
        th=np.arccos(np.clip(u@w,-1,1)); t=R[k]*np.tan(th/2)
        T1=b-t*u; T2=b+t*w; bis=(w-u)/np.linalg.norm(w-u); C=b+bis*R[k]/np.cos(th/2)
        a1=np.arctan2(*(T1-C)[::-1]); s=np.sign(u[0]*w[1]-u[1]*w[0])
        m=max(3,int(th*R[k]/h))
        for f in np.linspace(0,1,m): out.append(C+R[k]*np.array([np.cos(a1+s*f*th),np.sin(a1+s*f*th)]))
    return np.array(out)

def fit_radii(Q, verts, closed):
    return _merge_flats(Q, fit_radii_raw(Q, verts, closed), closed)

def fit_radii_raw(Q, verts, closed):
    """Radius per vertex: least-squares fit of the local composite to the points near the vertex."""
    n=len(verts); V=np.array([v[:2] for v in verts]); R=[0.0]*n
    ks=range(n) if closed else range(1,n-1)
    for k in ks:
        a,b,c=V[k-1],V[k],V[(k+1)%n]
        u=(b-a)/np.linalg.norm(b-a); w=(c-b)/np.linalg.norm(c-b); th=np.arccos(np.clip(u@w,-1,1))
        if th<np.radians(2): continue
        La,Lb=np.linalg.norm(b-a),np.linalg.norm(c-b)
        fa=0.5 if (closed or k-1>0) else 0.9; fb=0.5 if (closed or k+1<n-1) else 0.9
        tmax=min(fa*La,fb*Lb); rmax=tmax/np.tan(th/2)
        win=min(0.45*La,0.45*Lb)
        if not closed and k+1==n-1: win=0.45*La          # next to a notch end on the outline: the edge leg
        if not closed and k-1==0: win=0.45*Lb            # is straight, so only the TO side limits the window
        rmax=min(rmax,0.9*win/np.tan(th/2))      # tangent points stay on the local legs
        legs=dense(np.array([b-win*u,b,b+win*w]),closed=False)
        P=Q[cKDTree(legs).query(Q)[0]<win]
        if len(P)<3: continue
        def err(r):
            loc=rounded([(*(b-win*u),0),(*b,r),(*(b+win*w),0)],closed=False,h=0.02)
            return np.sqrt(np.mean(dist_to(loc,P,closed=False)**2))
        sa=np.sin((np.pi-th)/2)                        # sine of half the interior angle
        rlim=TIP_LOSS*sa/max(1e-9,1-sa)                # radius that shortens the tip by TIP_LOSS
        R[k]=minimize_scalar(err,bounds=(min(R_MIN,rlim,rmax),rmax),method='bounded',options={'xatol':1e-3}).x
    return [(v[0],v[1],r) for v,r in zip(verts,R)]

def _merge_flats(Q, out, closed):
    # a straight piece shorter than MIN_FLAT between two fillets: the two corners become one corner
    # (intersection of the outer sides) with one fillet fitted to the points
    def tang(P,k,m):
        a,b,c=P[k-1],P[k],P[(k+1)%m]; u=(b-a)/np.linalg.norm(b-a); w=(c-b)/np.linalg.norm(c-b)
        return np.tan(np.arccos(np.clip(u@w,-1,1))/2)
    for _ in range(len(out)):
        P=np.array([v[:2] for v in out],float); RR=[v[2] for v in out]; m=len(P); hit=None
        for k in (range(m) if closed else range(1,m-2)):
            j=(k+1)%m
            if m<=3 or RR[k]<=0 or RR[j]<=0 or (not closed and j==m-1): continue
            flat=np.linalg.norm(P[j]-P[k])-RR[k]*tang(P,k,m)-RR[j]*tang(P,j,m)
            if flat<MIN_FLAT: hit=(k,j); break
        if hit is None: break
        k,j=hit; a,b=P[k-1],P[k]; c,d=P[j],P[(j+1)%m]
        X=isect((a,(b-a)/np.linalg.norm(b-a)),(c,(d-c)/np.linalg.norm(d-c)))
        if X is None or np.linalg.norm(X-b)>3*np.linalg.norm(c-b)+5: break
        new=[v for i,v in enumerate(out) if i not in (k,j)]
        new.insert(k if j>k else len(new),(X[0],X[1],0.0))
        if j<k: new=new[-1:]+new[:-1]
        out=fit_radii_raw(Q,new,closed)
    return out

def corner_vertices(Q, runs, L, k_prev, k_next):
    """Vertex (or two) between side k_prev and side k_next. A near-parallel pair or an intersection far
    from the gap between the two runs becomes a short connecting side between the run ends."""
    (c1,d1),(c2,d2)=L[k_prev],L[k_next]
    e1=Q[runs[k_prev][1]%len(Q)]; s2=Q[runs[k_next][0]%len(Q)]
    p1=c1+((e1-c1)@d1)*d1; p2=c2+((s2-c2)@d2)*d2
    ang=np.degrees(np.arccos(np.clip(d1@d2,-1,1))); p=isect(L[k_prev],L[k_next])
    reach=np.linalg.norm(p2-p1)+1.0
    if p is None or ang<8 or np.linalg.norm(p-p1)>3*reach+2 or np.linalg.norm(p-p2)>3*reach+2 \
       or (p-p1)@d1<-reach or (p2-p)@d2<-reach:
        return [(p1[0],p1[1],0.0),(p2[0],p2[1],0.0)]
    return [(p[0],p[1],0.0)]

MERGE_PARALLEL=True
def fit_closed(ring, tol=0.15, Lmin=2.0):
    raw,Q=denoise(ring)
    r0=merge(Q,runs_of(Q,tol,Lmin))
    i,j=max(r0,key=lambda r:r[1]-r[0]); Q=np.roll(Q,-(i+j)//2,axis=0)   # seam inside the longest side
    Q2=np.vstack([Q,Q[:1]])
    runs=merge(Q2,runs_of(Q2,tol,Lmin))
    if MERGE_PARALLEL: runs=merge(Q2,runs,amin=20.0,tolm=2*TOLM)
    L=side_lines(Q2,runs)
    # first and last run are the two halves of the seam side when they are collinear
    if len(L)>2:
        d=np.degrees(np.arccos(np.clip(L[0][1]@L[-1][1],-1,1)))
        if d<8: runs=[[runs[-1][0],runs[0][1]+len(Q)]]+runs[1:-1]; Q3=np.vstack([Q,Q]); L=side_lines(Q3,runs)
    QQ=np.vstack([Q,Q])
    verts=[]
    for k in range(len(L)): verts+=corner_vertices(QQ,runs,L,k-1,k)
    verts=fit_radii(Q,verts,True)
    return raw, verts

OUTLINE=None   # original outline polygon, for the two-way check of notches
BLEND_ANGLE=30.0   # notch ends meeting the outline below this angle blend into it with a tangent fillet

def outline_edge(O_poly, p):
    """Line (point, unit direction) of the original outline segment nearest to p."""
    C=np.asarray(O_poly.exterior.coords); A,B=C[:-1],C[1:]; D=B-A; L2=(D*D).sum(1)
    t=np.clip(((p-A)*D).sum(1)/np.maximum(L2,1e-12),0,1); d=np.linalg.norm(A+t[:,None]*D-p,axis=1)
    k=int(np.argmin(d)); return A[k], D[k]/np.sqrt(L2[k]), (A[k],B[k])

def fit_open(Q, O_poly, tol=0.15, Lmin=2.0, out=1.5):
    runs=merge(Q,runs_of(Q,tol,Lmin))
    if MERGE_PARALLEL: runs=merge(Q,runs,amin=20.0,tolm=2*TOLM)
    L=side_lines(Q,runs)
    verts=[]
    # end vertices: extend the first and last side outward until clear of the original outline
    def extend(c,d,start):
        s=-1 if start else 1
        p0=Q[runs[0][0]] if start else Q[runs[-1][1]]
        base=c+((p0-c)@d)*d
        T=np.arange(0,30,0.05); X=base+s*T[:,None]*d
        ok=np.nonzero(inside_depth(O_poly,X)<=-out)[0]
        return X[ok[0]] if len(ok) else X[-1]
    def end(start):
        """End of the notch: a shallow end blends into the outline edge (vertex on the edge line,
        then a point further along the edge), a steep end is extended straight out."""
        c,d=L[0] if start else L[-1]
        p0=Q[runs[0][0]] if start else Q[runs[-1][1]]
        ce,de,(ea,eb)=outline_edge(O_poly,p0)
        # angle at which the TO curve itself meets the outline: its points within 3 mm of the outline
        dep=inside_depth(O_poly,Q); idx=np.arange(len(Q))
        near=idx[(dep>0.02)&(dep<3.0)]; near=near[near<len(Q)//2] if start else near[near>=len(Q)//2]
        if len(near)>=3:
            run=near[:np.argmax(np.diff(near)>1)+1] if start and np.any(np.diff(near)>1) else near
            if not start and np.any(np.diff(near)>1): run=near[np.nonzero(np.diff(near)>1)[0][-1]+1:]
            dto=lsq_line(Q[run])[1]
            ang=np.degrees(np.arccos(np.clip(abs(dto@de),0,1)))
        else:
            ang=np.degrees(np.arccos(np.clip(abs(d@de),0,1)))
        if ang<BLEND_ANGLE:
            V=isect((c,d),(ce,de))
            if V is not None:
                s=1 if start else -1          # travel direction of the notch, reversed at the start
                de=de if (de@d)*(-s)>0 else -de
                # along the edge, at most 10 mm and never past the end of that outline segment
                reach=max((ea-V)@de,(eb-V)@de); Lx=max(0.5,min(20.0,reach))
                # tilted slightly outward: the end leaves the outline instead of lying on it (a cut
                # face lying in the part's face makes OCC imprint the face instead of cutting)
                nrm=np.array([-de[1],de[0]])
                if inside_depth(O_poly,(V+0.5*nrm)[None])[0]>0: nrm=-nrm
                X=V+de*Lx+nrm*max(0.02,0.003*Lx)
                pts=[(V[0],V[1],0.0),(X[0],X[1],0.0)]
                return pts[::-1] if start else pts
        p=extend(c,d,start); return [(p[0],p[1],0.0)]
    verts=end(True)
    for k in range(1,len(L)): verts+=corner_vertices(Q,runs,L,k-1,k)
    verts+=end(False)
    return fit_radii(Q,verts,False)

def report(raw, C, closed=True):
    g=LinearRing(C) if closed else LineString(C)
    d=dist_to(C,raw,closed)
    # the other way too: a spike of the fit lies far from every TO point
    back=cKDTree(dense(raw,closed=closed,h=0.05)).query(dense(C,closed))[0]
    if not closed: back=back[inside_depth(OUTLINE,dense(C,closed=False))>0.3] if OUTLINE is not None else back
    return d.mean(), max(d.max(), back.max() if len(back) else 0), g.is_valid if closed else g.is_simple

def circle_fit(P):
    A=np.c_[2*P,np.ones(len(P))]; b=(P**2).sum(1); x=np.linalg.lstsq(A,b,rcond=None)[0]
    c=x[:2]; r=np.sqrt(x[2]+c@c); return c,r

def fit_open_arc(Q, O_poly, out=1.5):
    """Single circular arc through the notch curve, extended past the outline at both ends.
    Returns (centre, radius, start angle, signed sweep) and the max deviation of the inner points."""
    inner=Q[inside_depth(O_poly,Q)>0.3]
    c,r=circle_fit(inner); dev=np.abs(np.linalg.norm(inner-c,axis=1)-r).max()
    ang=np.unwrap(np.arctan2(Q[:,1]-c[1],Q[:,0]-c[0])); a0,a1=ang[0],ang[-1]; s=np.sign(a1-a0)
    # arc over the inner stretch only, then out to where the circle first meets the outline
    ai=np.unwrap(np.arctan2(inner[:,1]-c[1],inner[:,0]-c[0])); ai=ai+2*np.pi*np.round((ang[0]-ai[0])/(2*np.pi))
    lo,hi=(ai.min(),ai.max()) if s>0 else (ai.max(),ai.min())
    def first_cross(a, step):
        A=a+step*np.arange(0,1500); X=c+r*np.c_[np.cos(A),np.sin(A)]
        ok=np.nonzero(inside_depth(O_poly,X)<=0)[0]; return A[ok[0]] if len(ok) else None
    a0=first_cross(lo,-s*0.001); a1=first_cross(hi,s*0.001)
    if a0 is None or a1 is None: return None, np.inf
    return (c,r,a0,a1), dev
