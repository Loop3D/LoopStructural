## M1: speed with the real headless viewer (SC-001)

Purpose: check that an object is visible in the viewer within 5 s after a reconnect.

1. In PowerShell, set the token. Type `$env:GEOVIEWER_TOKEN = "secret"` and press Enter.
2. Start the viewer. Type `geoviewer-headless --port 0` and press Enter.
3. Read the port number that the viewer prints.
4. Save this script as a file outside the repository.

```python
import sys
import time
import numpy as np
from loop_link.errors import BusyError
from loop_link.objects import BlockModel
from loop_link.viewer import Viewer

PORT = int(sys.argv[1])  # the port that the viewer printed
grade = np.random.default_rng(0).random(1_000_000, dtype=np.float32)
block = BlockModel(origin=(0, 0, 0), cell_sizes=(1, 1, 1), cells=(100, 100, 100),
                   cell_properties={"grade": grade})

start = time.perf_counter()
viewer = Viewer(port=PORT, token="secret", crs="EPSG:28350")
viewer.send(block, id="grid")
viewer.release_data()

while True:
    other = Viewer(port=PORT, token="secret")
    try:
        other.connect()
    except BusyError:
        time.sleep(0.1)  # the viewer has not yet noticed the close
        continue
    found = any(entry.id == "grid" for entry in other.inventory())
    other.close()
    if found:
        break
    time.sleep(0.1)
print(f"Time until grid is in the inventory: {time.perf_counter() - start:.3f} s")
```

5. Run the script. Type `packages/loop_link/.venv/Scripts/python.exe <script path> <port>`, with the port from step 3, and press Enter.
6. Note that a `BusyError` is normal. The script waits 0.1 s and tries the connection again until the viewer accepts it.
7. Write the time that the script prints in your test record.
8. Stop the viewer. Press Ctrl+C.

Pass: the printed time is under 5 s.

## M2: seven types in the window viewer

Purpose: check that the window viewer shows seven objects, each with the right type and property.

1. Start the geoviewer window application.
2. Turn on the link in the viewer settings.
3. Save the token `secret` in the viewer settings.
4. Save this script as a file outside the repository.

```python
from loop_link.objects import BlockModel, Drillhole, Lines, Points, Surface, TetraMesh, Vectors
from loop_link.viewer import Viewer

grade3 = [1.0, 2.0, 3.0]
tri = [[0, 0, 0], [100, 0, 0], [0, 100, 0]]  # 100 units: the default glyph size is 1 percent of the object
objects = {
    "points": Points(tri, {"grade": grade3}),
    "vectors": Vectors.from_directions(tri, [[0, 0, 1]] * 3, kind="normal",
                                       properties={"grade": grade3}),
    "lines": Lines(tri, [[0, 1], [1, 2]], {"grade": grade3}),
    "surface": Surface(tri, [[0, 1, 2]], {"grade": grade3}),
    "block": BlockModel((0, 0, 0), (50, 50, 50), (2, 2, 2), cell_properties={"grade": [1.0] * 8}),
    "tetra": TetraMesh(tri + [[0, 0, 100]], [[0, 1, 2, 3]], {"grade": [1.0, 2.0, 3.0, 4.0]}),
    "drillhole": Drillhole([0, 0, 0], [[0, 0, 90], [100, 0, 90]], [[0, 50], [50, 100]],
                           interval_properties={"grade": [1.0, 2.0]}),
}

viewer = Viewer(token="secret", crs="EPSG:28350")
for name, obj in objects.items():
    viewer.send(obj, id=name)
viewer.flush()
viewer.release_data()
```

5. Run the script. Type `packages/loop_link/.venv/Scripts/python.exe <script path>` and press Enter.
6. Check that seven objects show in the object list. Check that the types are: points, vectors, lines, surface, block model, tetra mesh and drillhole.
7. Select each object, one at a time. Check that the property name `grade` shows in the property table.
   A tetra mesh has no Data button: use the isosurface or the filter controls to check `grade`.
   The vectors show as flat discs that face up (+z).
8. Write the result in your test record.

Pass: seven objects show, each with the right type, and the property name `grade` shows in the property table for each object.

## M3: select, read and fetch with the real window viewer (stories 2 to 4)

Purpose: check that the window viewer selects, reads and fetches cells, and that errors are reported.

Before you start, start the geoviewer window application. Turn on the link in the viewer settings. Save the token `secret` in the viewer settings. Save your test scripts as files outside the repository. Run them with `packages/loop_link/.venv/Scripts/python.exe <script path>`.

1. Run a script that sends a `BlockModel` of 10 x 10 x 10 cells. Then call `Viewer.select` with cells 12, 13 and 14. Check that the viewer outlines three cells. Check that `Viewer.selection()` returns these three cells.
2. Select cells in the window with the mouse. Call `Viewer.selection()`. Check that the returned cells are the cells that you selected.
3. Call `Viewer.select` with index 5000. Check that the script raises `InvalidInputError`. Check that the selection does not change.
4. Send the block model again to replace it. Then call `Viewer.select` with the old revision. Check that the script raises `StaleRevisionError`.
5. Fetch each of the seven types with the `fetch` method. Use the script from M2 to send the objects. Compare each fetched object with the original object.
6. Fetch a drillhole that the viewer imported. Check that the script or the viewer shows a warning for the dropped log sets.
7. Start `watch_selection`. Change the selection in the window. Check that the callback runs one time.
8. Start a second program. Check that its send raises `BusyError`. While the first program runs, check that it can still select and fetch.

M3 also checks the fake viewer assumptions OP12 (which types have a selection table), OP13 (one object per selection) and OP14 (parts arrive in order). Record the time of a fetch of 1,000,000 cells.

Pass: all 8 steps give the expected result.

## M4: filter and colour with the real window viewer (stories 6 and 7)

Purpose: check that the window viewer filters cells, sets colours and reports errors.

Before you start, start the geoviewer window application. Step 7 needs a viewer of protocol 1.20 (geoviewer spec 072). The other steps work with a viewer of protocol 1.18. Turn on the link in the viewer settings. Save the token `secret` in the viewer settings. Save your test scripts as files outside the repository. Run them with `packages/loop_link/.venv/Scripts/python.exe <script path>`. Start each script with `Viewer.clear()`. Place the test objects near the origin of the scene. An object far from the origin can give rendering artefacts (see M3).

1. Send a `BlockModel` with a number property `grade`. Call `filter_range("grid", "grade", 0.5, 1.0)`. Check that only the cells with a `grade` from 0.5 to 1.0 show.
2. Send an object with a category property. Call `filter_values` with two codes. Check that only the cells with these codes show. Call `filter_values` with the same two codes and `exclude=True`. Check that only the other cells show.
3. Call `filter_enable(False)`, then `filter_remove`. Check that all cells show.
4. Call `filter_range` for a text property. Then call `filter_range` for a property that does not exist. Check that the script raises `InvalidInputError` each time. Check the text of each error.
5. Call `set_colour_property`, `set_colour_map("coolwarm")`, `set_colour_range`, `set_colour_scale("log")` and `set_object_colour`. Check the picture after each call.
6. Call the colour map, colour range and colour scale methods for a category property. Check that each call raises an error. Call `colour_map_names()`. Check that it gives the names of the colour maps of the viewer.
7. Send the same id again (a new revision). Check that the filter and the colours stay. Then do these checks for the filter wait:
   - Send an object of 1,000,000 cells.
   - Call `filter_range("grid", "grade", 0.5, 1.0)` with the default `wait=True`. Check that the call returns after the build. Check that `FilterResult.passing` and `total` are numbers, and that `passing` fits the picture.
   - Call it again with `wait=False`. Check that it returns at once.
   - Set a `threading.Event` during a build. Check that the call raises `CancelledError` and that the filter is as before.
   - Record the time of the build.
8. Send an object of 1,000,000 cells with `colour=ColourSettings(...)`. Check that the settings apply. The 0.5 s of `flush` can be too short (OP19). If the settings do not apply, write the result in your test record.
9. Call `add` for a filter that already exists (OP16). Check the result. Then call `add` and `set_range` for a new filter. Check that the new filter is on.

M4 checks OP16, OP18 and OP19. Record the time of one `set_range` on the object of 1,000,000 cells.

Pass: all 9 steps give the expected result.

## Record

2026-11-07, window viewer of the geoviewer build with selection and fetch, token from the user:

- M1: pass. The grid was in the inventory 0.43 s after the send (limit 5 s). It ran against the window viewer, not the headless viewer.
- M2: pass after one fix.
  - Drillhole `grade` was not selectable. The viewer reads drillhole interval values only from arrays with association `interval` in a log set. The encoder sent `cell`. It now sends `interval` with the log set `intervals`.
  - Tetra mesh `grade` works through the isosurface. The viewer has no Data button for a tetra mesh.
  - Vectors were too small to see on a 1 unit triangle (glyph size default is 1 percent of the object). They work on a 100 unit triangle.
  - Open viewer issue: the vector glyphs and the drillhole tubes look pixelated (geoviewer issue 181).
- M3: partly passed on 2026-11-07, window viewer `ui-fixes-build`, token from the user. A person did steps 2 and 7. The scripts did the other steps.
  - Step 1 pass: select cells 12, 13 and 14 and `selection()` returns them. The first try failed. The scene had objects at map coordinates (origin 500000, 7000000) and the test block was at (0, 0, 0). The viewer showed rendering artefacts and `selection.get` returned nothing. After a fresh scene (`clear`), it worked. Put a test object near the scene origin, or use a fresh scene.
  - Step 2 pass: the user selected 11 rows with the mouse. `selection()` returned `m3-block`, revision 2, cells `[1, 2, 3, 4, 5, 6, 7, 8, 9, 45, 302]`.
  - Step 3 pass: index 5000 gives `InvalidInputError`. The selection does not change.
  - Step 4 pass: an old revision gives `StaleRevisionError`.
  - Step 5 pass: all seven types fetch with the `grade` property.
  - Step 6 partly passed: see the fix below. The warning for dropped properties works. The warning for dropped log sets was not seen, because the imported holes seem to have one log set.
  - Step 7 not decided: the callback ran 4 times, once for each different selection, with no repeats. The record does not say how many changes the user made.
  - Step 8 partly passed: the first program selected and fetched while a second program held the data connection. A real `BusyError` was not tested.
  - OP12: points and vectors take vertices only. Lines, surface, block model and tetra mesh take cells. A drillhole takes neither kind (no table tab).
  - OP13: one object has a selection at a time. A new selection on another object replaces the old one.
  - OP14: a fetch of 1,000,000 cells took 0.07 s and arrived complete. Several parts were not tested directly.
  - Fix: a viewer-imported drillhole (`import-1/0`) failed to fetch. Text columns such as `COMMENTS` come as `referenced` data with an empty `value_map`. `loop_link` now drops such a property and warns. A fetch takes about 1.8 s for a hole with about 100 properties.
  - Open viewer gap: a tool cannot select drillhole intervals or holes. A new geoviewer spec is needed.
