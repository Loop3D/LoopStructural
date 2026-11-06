# loop-link

Send geoscience data from Python to the LoopConnect viewer.

## Install

```
pip install loop-link
```

Optional extras: `pyvista`, `loopstructural`, `pandas`.

## Token

The package needs a token. It reads the `GEOVIEWER_TOKEN` environment
variable. If this variable is not set, it uses the token that the viewer
saved.

## Example

```python
import numpy as np
from loop_link import Viewer, Points

viewer = Viewer()
xyz = np.random.rand(100, 3)
viewer.connect()
viewer.send(Points(xyz), id="random-points", name="random points")
viewer.close()
```
