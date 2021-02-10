import numpy as np
class BaseStructuredSupport:
    """

    """
    def __init__(self,
                 origin=np.zeros(3),
                 nsteps=np.array([10, 10, 10]),
                 step_vector=np.ones(3),
                 ):
        """

        Parameters
        ----------
        origin - 3d list or numpy array
        nsteps - 3d list or numpy array of ints
        step_vector - 3d list or numpy array of int
        """
        # the geometry in the mesh can be calculated from the 
        # nsteps, step vector and origin
        # we use property decorators to update these when different parts of
        # the geometry need to change
        # inisialise the private attributes
        self._nsteps = nsteps
        self._step_vector = step_vector
        self._origin = np.array(origin)  
        self.supporttype='Base'

    @property
    def nsteps(self):
        return self._nsteps
    @nsteps.setter
    def nsteps(self,nsteps):
        #if nsteps changes we need to change the step vector
        change_factor = nsteps/self.nsteps
        self._step_vector/=change_factor
        self._nsteps = nsteps
        

    @property
    def step_vector(self):
        return self._step_vector
    
    @step_vector.setter
    def step_vector(self,step_vector):
        change_factor = step_vector/self.step_vector
        newsteps = self._nsteps/change_factor
        self._nsteps =np.ceil(newsteps).astype(int)
        self._step_vector = step_vector

    @property
    def origin(self):
        return self._origin
    @origin.setter
    def origin(self,origin):
        origin = np.array(origin)
        length = self.maximum-origin
        length /= self.step_vector
        self._nsteps = np.ceil(length).astype(int)
        self._origin = origin
    @property
    def maximum(self):
        return self.origin+self.nsteps*self.step_vector
    @maximum.setter
    def maximum(self,maximum):
        """
        update the number of steps to fit new boundary
        """
        maximum = np.array(maximum,dtype=float)
        length = maximum-self.origin
        length /= self.step_vector
        self._nsteps = np.ceil(length).astype(int)
    @property
    def n_nodes(self):
        return np.product(self.nsteps)
    
    @property
    def nsteps_cells(self):
        return self.nsteps -1
    @property
    def n_elements(self):
        return np.product(self.nsteps_cells)
    def __str__(self):
        return 'LoopStructural interpolation support:  {} \n'\
                'Origin: {} {} {} \n'\
                'Maximum: {} {} {} \n'\
                'Step Vector: {} {} {} \n'\
                'Number of Steps: {} {} {} \n'\
                'Degrees of freedon {}'.format(self.supporttype,self.origin[0],self.origin[1],self.origin[2],\
                                                    self.maximum[0],self.maximum[1],self.maximum[2],\
                                                      self.step_vector[0],self.step_vector[1],self.step_vector[2],\
                                                     self.nsteps[0],self.nsteps[1],self.nsteps[2],self.n_nodes)
    @property
    def nodes(self):
        max = self.origin + self.nsteps_cells * self.step_vector
        x = np.linspace(self.origin[0], max[0], self.nsteps[0])
        y = np.linspace(self.origin[1], max[1], self.nsteps[1])
        z = np.linspace(self.origin[2], max[2], self.nsteps[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        return np.array([xx.flatten(order='F'), yy.flatten(order='F'),
                               zz.flatten(order='F')]).T