2. Modelling Folds
-------------------
Standard implicit interpolation struggles to reproduce folded surfaces
from sparse data, because it only has a smoothness/regularisation term to
fill in between observations. These examples show how LoopStructural
instead uses a **fold frame** - a curvilinear coordinate system built
around the fold axis and axial surface - together with calculated fold
rotation angles to constrain folded and refolded (multiply-deformed)
surfaces directly from structural geology.
