import pyablo

reader   = pyablo.XdmfReader()

files = ("blast.xmf", "3dtest.xmf")

for file in files:
    snap = reader.readSnapshot(file)
    print(file)
    
    # t = [54,12,600]
    t = [54]
    v = snap.getCellsVolume(t)
    c = snap.getCellsCenter(t)

    print(v)
    print(c)
    print()
