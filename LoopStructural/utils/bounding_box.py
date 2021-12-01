# class BoundingBox:
#     def __init__(self,origin,maximum):
#         self.bb = np.array([origin,maximum])
#         self.maximum = maximum
#         self.name_map = {'xmin':(0,0),'ymin':(0,1),'zmin':(0,2),'xmax':(1,0),'ymax':(1,1),'zmax':(1,2)
#                          'lower':(0,2),'upper':(1,2),
#                          'minx':(0,0),'miny':(0,1),'minz':(0,2),'maxx':(1,0),'maxy':(1,1),'maxz':(1,2)}
#     def get_value(self,name):
#         ix,iy = self.name_map.get(name,(-1,-1))

#         return self.bb[ix,]

#     def is_inside(self,xyz):
#         inside = np.zeros(xy.shape[0],dtype=bool)
#         inside = np.logical_and(inside,xyz[:,0]>self.bb[0,0])
#         inside = np.logical_and(inside,xyz[:,0]<self.bb[1,0])
#         inside = np.logical_and(inside,xyz[:,1]>self.bb[0,1])
#         inside = np.logical_and(inside,xyz[:,1]<self.bb[1,1])
#         inside = np.logical_and(inside,xyz[:,2]>self.bb[0,2])
#         inside = np.logical_and(inside,xyz[:,2]<self.bb[1,2])
#         return inside
