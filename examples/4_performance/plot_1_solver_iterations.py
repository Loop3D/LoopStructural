"""
Solver Iteration Limits
-----------------------
LoopStructural uses iterative solvers to solve the least squares problem.
The number of iterations can be controlled by the user. 
"""

"""
Solver Iteration Limits
-----------------------
LoopStructural uses iterative solvers to solve the least squares problem.
The number of iterations can be controlled by the user. 
The number of iterations required by LoopStructural depends on the complexity of the geometry
being modelled.
In some cases a satisfactory solution can be found with a small number of iterations.
In the example below we show the time taken to solve the problem for different numbers of iterations
and show the resulting model.
"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.visualisation import Loop3DView
import time

data, bb = load_claudius()
model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
isovalue = 220
for maxiter in [10, 100, 500, 1000, 10000]:
    start = time.time()
    model.create_and_add_foliation(
        "strati", interpolatortype="FDI", solver_kwargs={"maxiter": maxiter}
    )
    model.update()
    end = time.time()
    print(f"maxiter: {maxiter} time: {end-start} seconds")
    view = Loop3DView(model)
    view.plot_surface(model["strati"], [220])
    view.show()
    view.clear()
