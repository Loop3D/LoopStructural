# from LoopStructural.interpolators import DiscreteInterpolator

# class NonLinearInterpolator(DiscreteInterpolator):
#     def __init__(self, support):
#         """
#         """
#         super().__init__(self, support)
#         self._non_linear_constraints = {}

#     def add_non_linear_constraint(self,):
#         self._non_linear_constraints[name] = constraint

#     def _solve(self,solver,niter=10,**kwargs):
#         A_base, B_base = feature.interpolator.build_matrix()
#         for i in range(niter):
#             for constraint in self._non_linear_constraints:
#                 ATA, ATB = constraint()
#                 A = ATA+A_base
#                 B = ATB+B_base
#         def iterative_norm(feature,niter,w=0.1):
#     values = np.zeros((niter,100,100))
#     norms = np.zeros((niter,100,100))
#     pts = np.array(np.meshgrid(np.linspace(model.bounding_box[0,0],model.bounding_box[1,0],100),np.linspace(model.bounding_box[0,1],model.bounding_box[1,1],100))).T.reshape((-1,2))
#     pts = np.hstack([pts,np.zeros((pts.shape[0],1))])
#     feature.interpolator.solve_system(solver='cg')
    
#     v = feature.evaluate_value(pts)
#     values[0,:,:] = v.reshape((100,100))
    
#     norm = np.linalg.norm(feature.evaluate_gradient(pts),axis=1)
#     norms[0,:,:] = norm.reshape((100,100))
#     print('Average gradient norm: {} std: {}'.format(np.nanmean(norm),np.nanstd(norm)))
#     A_base, B_base = feature.interpolator.build_matrix()

#     for i in range(1,niter):
#         print('Adding constant norm: {} / {}'.format(i+1,niter))
#         ATA, ATB = constant_norm_constraint(feature)
        
#         ATA*=w*((i**2+1)/(niter**2))
        
        
#         A = ATA+A_base
#         B = ATB+B_base
#         from scipy.sparse.linalg import cg, splu
 
#         soln = cg(A,B,tol=feature.builder.build_arguments['tol'])
#         feature.interpolator.c = soln[0]
#         v = feature.evaluate_value(pts)

#         values[i,:,:] = v.reshape((100,100))

#         norm = np.linalg.norm(feature.evaluate_gradient(pts),axis=1)
#         norms[i,:,:] = norm.reshape((100,100))

#         print('Average gradient norm: {} std: {}'.format(np.nanmean(norm),np.nanstd(norm)))
  
#     return values, norms
    
