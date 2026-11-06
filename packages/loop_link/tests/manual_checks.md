## M1: speed with the real headless viewer (SC-001)

Purpose: check that an object is visible in the viewer within 5 s after a reconnect.

1. In PowerShell, set the token: `$env:GEOVIEWER_TOKEN = "secret"`.
2. Start the viewer: `geoviewer-headless --port 0`.
3. Read the port number that the viewer prints.
4. Save this script. Replace `0` with that port number.

```python
import time
import numpy as np
from loop_link.objects import BlockModel
from loop_link.viewer import Viewer

PORT = 0  # replace with the port that the viewer printed
grade = np.random.default_rng(0).random(1_000_000, dtype=np.float32)
block = BlockModel(origin=(0, 0, 0), cell_sizes=(1, 1, 1), cells=(100, 100, 100),
                   cell_properties={"grade": grade})

viewer = Viewer(port=PORT, token="secret", crs="EPSG:28350")
viewer.send(block, id="grid")
viewer.release_data()

start = time.perf_counter()
while True:
    other = Viewer(port=PORT, token="secret")
    other.connect()
    found = any(entry.id == "grid" for entry in other.inventory())
    other.close()
    if found:
        break
    time.sleep(0.1)
print(f"Time until grid is in the inventory: {time.perf_counter() - start:.3f} s")
```

5. Run the script: `packages/loop_link/.venv/Scripts/python.exe <script path>`.
6. Read the time that the script prints.
7. Write the time in your test record.
8. Stop the viewer. Press Ctrl+C.

Pass: the printed time is under 5 s.

## M2: seven types in the window viewer

Purpose: check that the window viewer shows seven objects, each with the right type and property.

1. Start the viewer. Open its window.
2. Turn on the link in the viewer settings.
3. Save the token `secret` in the viewer settings.
4. Save this script and run it with `packages/loop_link/.venv/Scripts/python.exe <script path>`.

```python
from loop_link.objects import BlockModel, Drillhole, Lines, Points, Surface, TetraMesh, Vectors
from loop_link.viewer import Viewer

grade3 = [1.0, 2.0, 3.0]
tri = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
objects = {
    "points": Points(tri, {"grade": grade3}),
    "vectors": Vectors.from_directions(tri, [[0, 0, 1]] * 3, kind="normal",
                                       properties={"grade": grade3}),
    "lines": Lines(tri, [[0, 1], [1, 2]], {"grade": grade3}),
    "surface": Surface(tri, [[0, 1, 2]], {"grade": grade3}),
    "block": BlockModel((0, 0, 0), (1, 1, 1), (2, 2, 2), cell_properties={"grade": [1.0] * 8}),
    "tetra": TetraMesh(tri + [[0, 0, 1]], [[0, 1, 2, 3]], {"grade": [1.0, 2.0, 3.0, 4.0]}),
    "drillhole": Drillhole([0, 0, 0], [[0, 0, 0], [0, 0, -10]], [[0, 5], [5, 10]],
                           interval_properties={"grade": [1.0, 2.0]}),
}

viewer = Viewer(token="secret", crs="EPSG:28350")
for name, obj in objects.items():
    viewer.send(obj, id=name)
viewer.flush()
viewer.release_data()
```

5. Check that seven objects show in the object list.
6. Check that the types are: points, vectors, lines, surface, block model, tetra mesh and drillhole.
7. Select each object. Check that the property name `grade` shows in the property table.
8. Write the result in your test record.

Pass: seven objects show, each with the right type, and each property name shows in the property table.
