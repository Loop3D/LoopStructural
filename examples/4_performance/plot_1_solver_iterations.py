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
We can see that for a lower number of iterations the solution is not satisfactory but as the number of 
iterations increases the solution beocmes satisfactory.
"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.visualisation import Loop3DView
import time

data, bb = load_claudius()
model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
isovalue = 0
images = {}
view = Loop3DView(model, off_screen=True)

for maxiter in [10, 100, 500, 1000, 10000]:
    start = time.time()
    model.create_and_add_foliation(
        "strati", interpolatortype="FDI", solver_kwargs={"maxiter": maxiter}
    )
    model.update()
    end = time.time()
    print(f"maxiter: {maxiter} time: {end-start} seconds")
    view.plot_surface(model["strati"], [0, 60, 250, 330])
    view.plot_data(model["strati"])

    images[maxiter] = view.screenshot()
    view.clear()

data = data.drop(data.index[data["val"] == 250])
data = data.drop(data.index[data["val"] == 330])
images2 = {}
model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
for maxiter in [10, 100, 500, 1000, 10000]:
    start = time.time()
    model.create_and_add_foliation(
        "strati", interpolatortype="FDI", solver_kwargs={"maxiter": maxiter}
    )
    model.update()
    end = time.time()
    print(f"maxiter: {maxiter} time: {end-start} seconds")
    view.plot_surface(model["strati"], [0, 60], name='surface')
    view.plot_data(model["strati"], name='data')
    images2[maxiter] = view.screenshot()
    view.clear()
import matplotlib.pyplot as plt

fig, ax = plt.subplots(len(images), 2, figsize=(10, 20))
for i, (maxiter) in enumerate(images.keys()):
    ax[i, 0].imshow(images[maxiter])
    ax[i, 0].set_title(f"maxiter: {maxiter}")
    ax[i, 1].imshow(images2[maxiter])
    ax[i, 1].set_title(f"maxiter: {maxiter}")
