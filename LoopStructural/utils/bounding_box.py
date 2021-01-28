# class BoundingBox:
#     def __init__(self,origin,maximum):
#         self.origin = origin
#         self.maximum = maximum
#         self.name_map = {'xmin':(0,0),'ymin':(0,1),'zmin':(0,2),'xmax':(1,0),'ymax':(1,1),'zmax':(1,2)
#                          'lower':(0,2),'upper':(1,2),
#                          'minx':(0,0),'miny':(0,1),'minz':(0,2),'maxx':(1,0),'maxy':(1,1),'maxz':(1,2)}
#     def get_value(self,name):
#         ix,iy = self.name_map.get(name,(-1,-1))            
        