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

## Select and fetch

```python
from loop_link import Viewer, BlockModel

grid = BlockModel((0, 0, 0), (1, 1, 1), (4, 4, 4))
viewer = Viewer()
viewer.send(grid, id="grid")
viewer.select("grid", cells=[12, 13])
sel = viewer.selection()
print(sel.object_id, sel.cells)
watcher = viewer.watch_selection(print, interval=0.5)
watcher.stop()
result = viewer.fetch("grid")
print(result.object)
```

Select, read and fetch use the tool connection. They work while another program holds the data connection.
A viewer without geoviewer protocol 1.18 gives `MissingFeatureError`.
A fetch of a drillhole keeps the first interval set. It warns about dropped properties.

## Filter and colour

```python
from loop_link import BlockModel, ColourSettings, Viewer

grid = BlockModel((0, 0, 0), (1, 1, 1), (4, 4, 4))
viewer = Viewer()
viewer.send(grid, id="grid", colour=ColourSettings(array="grade", colour_map="viridis", range=(0.0, 1.0), scale="linear"))
result = viewer.filter_range("grid", "grade", 0.5, 1.0)
print(result.passing, result.total)
viewer.filter_values("grid", "unit", [0, 2], exclude=True)
viewer.filter_remove("grid", "grade")
viewer.set_colour_map("grid", "coolwarm")
print(viewer.colour_map_names())
```

Filter and colour use the tool connection. They work while another program holds the data connection.
The colour map names come from the viewer. A viewer without geoviewer protocol 1.18 gives `MissingFeatureError` for `colour_map_names()`.
A filter call waits for the end of the build by default. `wait=False` does not wait, `cancel` is a `threading.Event`, and `progress` gets the fraction.
The counts need a viewer of protocol 1.20. An older viewer gives a `FilterResult` with no counts.
`state` is `building` when the viewer stopped to wait (the build goes on) and `unbuilt` when the viewer has no counts for the new filters.
A value filter uses the integer ids of the categories. These are the codes of `Categorical`.

## Test against the real codec

Some tests give the frames of `loop_link` to the real geoviewer codec.
The geoviewer repository is private, so these tests are optional.
Set `LOOP_LINK_GEOVIEWER_REPO` to a geoviewer checkout and install Rust (`cargo`).
The tests build a small validator program the first time that they run.
Without the variable, or without `cargo`, these tests skip.
