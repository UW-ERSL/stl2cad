import sys, time; sys.path.insert(0, '.')  # CadQuery reference of the skeleton model; run from to/
import trimesh, numpy as np, swplates, hubweb, skeleton, cadquery as cq
ns=swplates._namespace('out_CrankArm/CrankArm_exact.py')
T=trimesh.load('in/CrankArm_TO_1.stl'); T.apply_scale(25.4)
hw=hubweb.detect(ns['original'], T); assert hw['axis']==2
S=skeleton.fit(hw['shadow'], hw['hubs'], T, hw['axis'])
J=np.array(S['joints']); hubs=hw['hubs']
solids=[]
for h in hubs:
    solids.append(cq.Workplane().workplane(offset=h['lo']).center(*h['c']).circle(h['ro']).extrude(h['hi']-h['lo']))
top={}
for b in S['bars']:
    p,q=J[b['i']].copy(),J[b['j']].copy(); w,d,c=b['width'],b['depth'],b['centre']
    u=(q-p)/np.linalg.norm(q-p)
    for k,e,s in ((b['i'],p,-1),(b['j'],q,1)):
        if k in S['hub_of']: e += s*u*w/2      # reach into the hub
    L=np.linalg.norm(q-p); m=(p+q)/2; ang=np.degrees(np.arctan2(u[1],u[0]))
    solids.append(cq.Workplane().box(L,w,d).rotate((0,0,0),(0,0,1),ang).translate((m[0],m[1],c)))
    for k in (b['i'],b['j']):
        if k not in S['hub_of']:
            r,lo,hi=top.get(k,(0,1e9,-1e9)); top[k]=(max(r,w/2),min(lo,c-d/2),max(hi,c+d/2))
for k,(r,lo,hi) in top.items():
    solids.append(cq.Workplane().workplane(offset=lo).center(*J[k]).circle(r).extrude(hi-lo))
u=solids[0]
for s in solids[1:]: u=u.union(s)
for h in hubs:
    u=u.cut(cq.Workplane().workplane(offset=h['lo']-1).center(*h['c']).circle(h['ri']).extrude(h['hi']-h['lo']+2))
cq.exporters.export(u,'out_CrankArm/crank_sk.stl',tolerance=0.02,angularTolerance=0.2)
cq.exporters.export(u,'out_CrankArm/crank_sk.step')
m=trimesh.load('out_CrankArm/crank_sk.stl'); to=T
print('bars %d joints %d (hub %d)'%(len(S['bars']),len(J),len(S['hub_of'])))
print('volume: TO %.0f, skeleton %.0f mm^3 (%+.1f %%)'%(to.volume,m.volume,100*(m.volume/to.volume-1)))
lo=np.minimum(to.bounds[0],m.bounds[0])-1; pitch=1.0
shape=tuple(np.ceil((np.maximum(to.bounds[1],m.bounds[1])+1-lo)/pitch).astype(int)+1)
def vox(mesh):
    g=mesh.voxelized(pitch).fill(); idx=np.argwhere(g.matrix)
    pts=trimesh.transformations.transform_points(idx,g.transform)
    G=np.zeros(shape,bool); k=np.round((pts-lo)/pitch).astype(int); G[k[:,0],k[:,1],k[:,2]]=True; return G
A,B=vox(to),vox(m); both=(A&B).sum()
print('voxel overlap: %.0f %% of TO in model, %.0f %% of model in TO'%(100*both/A.sum(),100*both/B.sum()))
