# class LinearTetrahedron:
#     def __init__(self):
#         pass

#     def element_gradient(self, elements):
#         ps = self.
#         m = np.array(
#             [[(ps[:, 1, 0] - ps[:, 0, 0]), (ps[:, 1, 1] - ps[:, 0, 1]),
#               (ps[:, 1, 2] - ps[:, 0, 2])],
#              [(ps[:, 2, 0] - ps[:, 0, 0]), (ps[:, 2, 1] - ps[:, 0, 1]),
#               (ps[:, 2, 2] - ps[:, 0, 2])],
#              [(ps[:, 3, 0] - ps[:, 0, 0]), (ps[:, 3, 1] - ps[:, 0, 1]),
#               (ps[:, 3, 2] - ps[:, 0, 2])]])
#         I = np.array(
#             [[-1., 1., 0., 0.],
#              [-1., 0., 1., 0.],
#              [-1., 0., 0., 1.]])
#         m = np.swapaxes(m, 0, 2)
#         element_gradients = np.linalg.inv(m)

#         element_gradients = element_gradients.swapaxes(1, 2)
#         element_gradients = element_gradients @ I

#         return element_gradients[elements,:,:]
