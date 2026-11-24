"""Test stand-in for pythoncom (see win32com/client.py)."""
VT_EMPTY = 0
VT_DISPATCH = 9
VT_ARRAY = 8192
VT_R8 = 5


def PumpWaitingMessages():
    return 0
