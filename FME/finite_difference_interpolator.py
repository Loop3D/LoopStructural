from scipy.sparse import linalg
from scipy.sparse import coo_matrix, spdiags

class FDI(GeologicalInterpolator):
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
    def _setup_interpolator(self, **kwargs):


    def add_regularisation(self, w=0.1):

        return

    def add_value_ctr_points(self, w=1.0):  # for now weight all value points the same
        for p in self.p_i:
            element, flag = self.mesh.get_element(p.pos)
            if flag == False:
                print('Could not find triangle for x:%f y:%f z:%f' % (p.pos[0], p.pos[1], p.pos[2]))
                continue
            if ~np.all(self.region[element]):
                print('Could not find element for x:%f y:%f z:%f inside region' % (p.pos[0], p.pos[1], p.pos[2]))

                continue


    def add_gradient_ctr_pts(self, w=1.0):  # for now weight all gradient points the same
        for p in self.p_g:

    def add_tangent_ctr_pts(self, w=1.0):
        for p in self.p_t:


    def add_elements_gradient_orthogonal_constraint(self, elements, normals, w=1.0, B=0):

    def add_elements_gradient_constraint(self, elements, normal, w):
        if self.params['shape'] == 'rectangular':
            cols = self.region_map[np.array(self.col)]

            self.AA = coo_matrix((np.array(self.A), (np.array(self.row), \
                                                     cols)), shape=(self.c_, self.nx), dtype=float)  # .tocsr()
        if self.params['shape'] == 'square':
            cols = np.array(self.col)
            rows = np.array(self.row)
            self.AA = coo_matrix((np.array(self.A), (np.array(rows).astype(np.int64), \
                                                     np.array(cols).astype(np.int64))), dtype=float)
            d = np.zeros(self.nx)
            d += np.finfo('float').eps
            self.AA += spdiags(d, 0, self.nx, self.nx)
        B = np.array(self.B)
        self.cc_ = [0, 0, 0, 0]
        self.c = np.zeros(self.mesh.n_nodes)
        self.c[:] = np.nan
        if solver == 'lsqr':
            self.cc_ = sla.lsqr(self.AA, B)

    def _solve(self, solver='lsqr'):
        if solver == 'lsqr':

        self.c[self.region] = self.cc_[0]
