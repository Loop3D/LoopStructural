from scipy.sparse import linalg
from scipy.sparse import coo_matrix, spdiags
from .discete_interpolator import DiscreteInterpolator
class FDI(DiscreteInterpolator):
    """
    Finite Difference Interpolator
    """
    def __init__(self, **kwargs):
        #define cartesian grid
        self.params = {}
        self.params['shape'] = 'rectangular'
        self.params['nx'] = 10
        self.params['ny'] = 10
        self.params['nz'] = 10
        self.AA = []
        self.B = []
        self.constraints = []

    def _solve(self, solver='lsqr'):
        if solver == 'lsqr':

        self.c[self.region] = self.cc_[0]
