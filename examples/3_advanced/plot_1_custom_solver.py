"""
Custom solver
--------------

This example shows how to use a custom solver to solve the optimization problem.
"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from scipy.sparse.linalg import lsqr

data, bb = load_claudius()
model = GeologicalModel(bb[0, :], bb[1, :])
model.set_model_data(data)


def solver(A, b):
    print("Using custom solver")
    res = lsqr(A, b)
    return res[0]


model.create_and_add_foliation("strati", solver=solver)
model.update()

##################################################################################################
# Building interpolation matrix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# There are two matrices that are built by LoopStructural. The first is the interpolation matrix
# which is used to interpolate the scalar field using only linear constraints. The second is the
# inequality matrix which is used incorporate the lower and upper bound constraints. We can
# access both matrices and solve the system of equations outside of LoopStructural and then update
# the model with the new solution.
##################################################################################################

A, b = model['strati'].interpolator.build_matrix()
print(A.shape, b.shape)
Q, bounds = model['strati'].interpolator.build_inequality_matrix()

# solve the system of equations
c = solver(A, b)

model['strati'].interpolator.c = c
