try:
    from trame_pyvista import jupyter
except ImportError:
    from pyvista.trame import jupyter
from .trame import initialize


jupyter.initialize = initialize
