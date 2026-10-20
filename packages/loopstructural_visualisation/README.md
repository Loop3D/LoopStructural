# loopstructural-visualisation

A LoopStructural interface for pyvista's Plotter class.

The viewer uses the current LoopStructural geometry and stratigraphic-column
APIs. Feature evaluation takes world coordinates; the viewer does not apply an
additional model-coordinate transform. Notebook integration supports both
PyVista's original Trame modules and the newer `trame-pyvista` package.

Compatibility checked on Windows with Python 3.13 and workspace LoopStructural
1.8.4:

| VTK | PyVista | Checks |
| --- | --- | --- |
| 9.7.1 | 0.49.1 | Core plotting, off-screen rendering, clipping, notebook UI |
| 9.6.2 | 0.48.4 | Core plotting, off-screen rendering, clipping |

CI also runs the plotting and notebook tests on Windows and Linux against the
VTK 9.6 and 9.7 release lines with the latest compatible PyVista.

From the workspace root, run the tests with:

```sh
uv sync --package loopstructuralvisualisation --extra tests --extra jupyter
uv run --no-sync pytest packages/loopstructural_visualisation/tests
```


To install `pip install loopstructuralvisualisation` or for a jupyter notebook environment (including the pyvista[jupyuter] dependencies) `pip install loopstructuralvisualisation[jupyter]`

