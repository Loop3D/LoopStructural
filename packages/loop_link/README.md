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

## Test against the real codec

Some tests give the frames of `loop_link` to the real geoviewer codec.
The geoviewer repository is private, so these tests are optional.
Set `LOOP_LINK_GEOVIEWER_REPO` to a geoviewer checkout and install Rust (`cargo`).
The tests build a small validator program the first time that they run.
Without the variable, or without `cargo`, these tests skip.
