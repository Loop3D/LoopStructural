from FME.supports.structured_grid import StructuredGrid


grid = StructuredGrid()
print(grid.centers())


import lavavu
import numpy as np

lv = lavavu.Viewer(background='grey')

cube = np.arange(0, 9, dtype=np.float32).reshape(3,3)

vol = lv.volume("cube", dims=cube.shape, vertices=[[0,0,0],[1,1,1]])
vol.values(cube)

colourmap = vol.colourmap("red white blue")
colourmap.colours[1][3] = 0.0 # alpha channel
cbar = vol.colourbar()

vol.control.ColourMaps()
lv.control.Panel()
vol.control('opacity')
vol.control('density')
vol.control('power')
vol.control('samples')

lv.control.show()
lv.interactive()