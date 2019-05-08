import numpy as np
class StructuralFrame():
    """
    Structural frame is a curvilinear coordinate system for representing
    the fault frame or fold frame. It is made up of 2/3 scalar fields and their 
    related direction vectors.
    """
    def __init__(self,age,mesh,z,y,x=None):
        self.age = age
        self.mesh = mesh
        self.z = z
        self.y = y
        self.compute_x = True
        self.x_cache = True
        if x is None:
            self.compute_x = True
            self.x_cache = False
        else:
            self.compute_x = False
            self.x = x
            self.x_cache = False
    def get_x_value(self,pos=None,ele=None):
        if self.compute_x:
            #print("X is analytical and there is no scalar field computed")
            return False
        return self.mesh.eval_interpolant(array=pos,prop=self.x,k=50,e=ele)
    def get_y_value(self,pos=None,ele=None):
        #return self.mesh.eval_interpolant(array=pos,prop=self.x,k=50,e=ele)
        return self.mesh.get_property_value(pos=pos,element=ele,prop=self.y)        
    def get_z_value(self,pos=None,ele=None):
        return self.mesh.get_property_value(pos=pos,element=ele,prop=self.z)
    def get_x_dir(self,pos=None,ele=None):
        #if self.x_cache == True: #if x is already computed don't be silly
        #    return self.mesh.get_element_property_gradient(pos=pos,element=ele,prop=self.x)
        if self.compute_x:
            return np.cross(self.get_z_dir(pos=pos,ele=ele),self.get_y_dir(pos=pos,ele=ele))

        return(self.mesh.get_element_property_gradient(pos=pos,element=ele,prop=self.x))  
    def get_y_dir(self,pos=None,ele=None):     
        return(self.mesh.get_element_property_gradient(pos=pos,element=ele,prop=self.y))
    def get_z_dir(self,pos=None,ele=None):
        return(self.mesh.get_element_property_gradient(pos=pos,element=ele,prop=self.z))
    def __call__(self,pos=None,ele=None):
        pass
    def add_x_to_mesh(self,name='X'):
        self.x = name
        grad = np.zeros((self.mesh.elements.shape[0],3))
        for i, e in enumerate(self.mesh.elements):
            grad[i,:] = self.get_x_dir(ele=e)
        self.mesh.property_gradients[name] = grad
        self.x_cache = True
            