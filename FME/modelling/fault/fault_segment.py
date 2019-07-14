import numpy as np

from FME.modelling.structuralframes import StructuralFrame


class FaultSegment(StructuralFrame):
    """
    Class for representing a slip event of a fault
    """

    def __init__(self, **kwargs):  # mesh,fault_event,data,name,region):
        """
        mesh:  support for interpolation
        fault_event: the parent fault that this segment belongs to
        data: the data that are associated with this segment
        name: a name for this segment
        region: region where the interpolation occurs for this segment
        """
        super().__init__(**kwargs)
        self.fault_event = None
        if 'fault_event' in kwargs:
            self.fault_event = kwargs['fault_event']
        if 'name' in kwargs:
            self.name = kwargs['name']

    def overlapFaultSegments(self, faultsegment, region):
        # here just find the
        self.overlap = {'segment': faultsegment, 'region': region}

    def buildFrame(self, **kwargs):
        self.buildFaultFrame(**kwargs)

    def buildFaultFrame(self, solver='lsqr', **kwargs):
        """
        Build the fault frame for this segment using the solver specified, default is scipy lsqr
        
        """
        # determine region
        overlap = False
        gxxgy = 1.
        gxxgz = 1.
        gyxgz = 1.
        gxcg = 0.1
        gycg = 0.1
        gzcg = 0.1
        gxcp = 1
        gycp = 1
        gzcp = 1
        gxgcp = 1
        gygcp = 1
        gzgcp = 1
        gx = True
        gy = True
        gz = True
        if 'gx' in kwargs:
            gx = kwargs['gx']
        if 'gy' in kwargs:
            gy = kwargs['gy']
        if 'gz' in kwargs:
            gz = kwargs['gz']
        if 'gxxgy' in kwargs:
            gxxgy = kwargs['gxxgy']
        if 'gxxgz' in kwargs:
            gxxgz = kwargs['gxxgz']
        if 'gyxgz' in kwargs:
            gyxgz = kwargs['gyxgz']
        if 'gxcg' in kwargs:
            gxcg = kwargs['gxcg']
        if 'gycg' in kwargs:
            gycg = kwargs['gycg']
        if 'gzcg' in kwargs:
            gzcg = kwargs['gzcg']
        if 'gxcp' in kwargs:
            gxcp = kwargs['gxcp']
        if 'gycp' in kwargs:
            gycp = kwargs['gycp']
        if 'gzcp' in kwargs:
            gzcp = kwargs['gzcp']
        if 'gxgcp' in kwargs:
            gxgcp = kwargs['gxgcp']
        if 'gygcp' in kwargs:
            gygcp = kwargs['gygcp']
        if 'gzgcp' in kwargs:
            gzgcp = kwargs['gzgcp']
        if 'segment' in self.overlap:
            overlapsegment = self.overlap['segment']
            overlapregion = self.overlap['region']
            overlap = True
        shape = 'rectangular'
        if 'shape' in kwargs:
            shape = kwargs['shape']

        for d in self.data:
            if d['type'] == 'gx':
                self.interpolators['gx'].add_data(d['data'])
            if d['type'] == 'gy':
                self.interpolators['gy'].add_data(d['data'])
            if d['type'] == 'gz':
                self.interpolators['gz'].add_data(d['data'])
        self.interpolators['gx'].setup_interpolator(cgw=gxcg, cpw=gxcp, gpw=gxgcp)
        if overlap:
            self.mesh.indices = np.array(range(0, self.mesh.n_nodes))
            sharedindices = self.mesh.indices[overlapregion]
            node_values = self.mesh.properties[overlapsegment.interpolators['gx'].propertyname] \
                [self.mesh.indices[overlapregion]]
            As = np.zeros(node_values.shape)
            As[:] = 1
            rows = np.array(range(self.interpolators['gx'].c_, self.interpolators['gx'].c_ + len(sharedindices)))
            self.interpolators['gx'].A.extend(As.tolist())
            self.interpolators['gx'].col.extend(sharedindices.tolist())
            if shape == 'rectangular':
                self.interpolators['gx'].B.extend(node_values)
                self.interpolators['gx'].row.extend(rows.tolist())
                self.interpolators['gx'].c_ += +len(sharedindices)
            if shape == 'square':
                self.interpolators['gx'].B[overlapregion] += node_values
                self.interpolators['gx'].row.extend(sharedindices.tolist())

        if gx:
            self.interpolators['gx'].solve_system(solver=solver)
            self.mesh.update_property(self.name + '_' + 'gx', self.interpolators['gx'].c)
            self.interpolators['gy'].add_elements_gradient_orthogonal_constraint(np.arange(0, \
                                                                                           self.mesh.n_elements), \
                                                                                 self.mesh.property_gradients[
                                                                                     self.name + '_' + 'gx'], w=gxxgy)

        self.interpolators['gy'].setup_interpolator(cgw=gycg, cpw=gycp, gpw=gygcp)  # cgw=0.1)#
        if overlap:
            self.mesh.indices = np.array(range(0, self.mesh.n_nodes))
            sharedindices = self.mesh.indices[overlapregion]
            node_values = self.mesh.properties[overlapsegment.interpolators['gy'].propertyname] \
                [self.mesh.indices[overlapregion]]
            As = np.zeros(node_values.shape)
            As[:] = 1
            rows = np.array(range(self.interpolators['gy'].c_, self.interpolators['gy'].c_ + len(sharedindices)))
            self.interpolators['gy'].A.extend(As.tolist())
            self.interpolators['gy'].col.extend(sharedindices.tolist())

            if shape == 'rectangular':
                self.interpolators['gy'].B.extend(node_values)
                self.interpolators['gy'].row.extend(rows.tolist())
                self.interpolators['gy'].c_ += +len(sharedindices)
            if shape == 'square':
                self.interpolators['gy'].B[overlapregion] += node_values
                self.interpolators['gy'].row.extend(sharedindices.tolist())
        if gy:
            self.interpolators['gy'].solve_system(solver=solver)
            self.mesh.update_property(self.name + '_' + 'gy', self.interpolators['gy'].c)
            self.interpolators['gz'].add_elements_gradient_orthogonal_constraint(np.arange(0, self.mesh.n_elements) \
                                                                                 , \
                                                                                 self.mesh.property_gradients[
                                                                                     self.name + '_' + 'gx']
                                                                                 , w=gxxgz)
            self.interpolators['gz'].add_elements_gradient_orthogonal_constraint(np.arange(0, self.mesh.n_elements) \
                                                                                 , self.mesh.property_gradients[
                                                                                     self.name + '_' + 'gy'], w=gyxgz)

        self.interpolators['gz'].setup_interpolator(cgw=gzcg, cpw=gzcp, gpw=gzgcp)  # cgw=0.1)
        if overlap:
            self.mesh.indices = np.array(range(0, self.mesh.n_nodes))
            sharedindices = self.mesh.indices[overlapregion]
            node_values = self.mesh.properties[overlapsegment.interpolators['gz'].propertyname] \
                [self.mesh.indices[overlapregion]]
            As = np.zeros(node_values.shape)
            As[:] = 1
            rows = np.array(range(self.interpolators['gz'].c_, self.interpolators['gz'].c_ + len(sharedindices)))
            self.interpolators['gz'].A.extend(As.tolist())
            self.interpolators['gz'].col.extend(sharedindices.tolist())

            if shape == 'rectangular':
                self.interpolators['gz'].B.extend(node_values)
                self.interpolators['gz'].row.extend(rows.tolist())
                self.interpolators['gz'].c_ += +len(sharedindices)
            if shape == 'square':
                self.interpolators['gz'].B[overlapregion] += node_values
                self.interpolators['gz'].row.extend(sharedindices.tolist())
        if gz:
            self.interpolators['gz'].solve_system(solver=solver)
            self.mesh.update_property(self.name + '_' + 'gz', self.interpolators['gz'].c)

    def buildFaultDisplacementField(self, **kwargs):
        gx = self.mesh.properties[self.interpolators['gx'].propertyname]
        if 'faultfunction' in kwargs:
            gy = self.mesh.properties[self.interpolators['gy'].propertyname]
            gz = self.mesh.properties[self.interpolators['gz'].propertyname]
        gx0 = 0.
        gy0 = 0.
        gz0 = 0.
        gxl = 1.
        gyl = 1.
        gzl = 1.
        dm = 2.
        d = np.zeros(self.mesh.n_nodes)
        if 'gx0' in kwargs:
            gx0 = kwargs['gx0']
        if 'gy0' in kwargs:
            gy0 = kwargs['gy0']
        if 'gz0' in kwargs:
            gz0 = kwargs['gz0']

        if 'gxl' in kwargs:
            gxl = kwargs['gxl']
        if 'gyl' in kwargs:
            gyl = kwargs['gyl']
        if 'gzl' in kwargs:
            gzl = kwargs['gzl']
        if 'faultfunction' in kwargs:
            fault_function = kwargs['faultfunction']
            gxn, gyn, gzn, d = fault_function(gx, gy, gz, gx0, gy0, gz0, gxl, gyl, gzl)
        if 'add_normalised' in kwargs:
            self.mesh.update_property(self.interpolators['gx'].propertyname + '_n', gxn)
            self.mesh.update_property(self.interpolators['gy'].propertyname + '_n', gyn)
            self.mesh.update_property(self.interpolators['gz'].propertyname + '_n', gzn)

        if 'dm' in kwargs:
            dm = kwargs['dm']
        if 'faultfunction' not in kwargs:
            d[np.logical_and(~np.isnan(d), gx > 0)] = 1.
            d[np.logical_and(~np.isnan(d), gx < 0)] = 0.
        if 'gy' in kwargs:
            d[gy < kwargs['gy']] = 0.
        if 'boundary' in kwargs and 'boundary_value' in kwargs:
            d[kwargs['boundary'] < kwargs['boundary_value']] = 0
        if 'mask' in kwargs:
            d[kwargs['mask']] = 0
        d *= dm

        # d*=20.
        self.mesh.update_property(self.name + '_d', d)

    def applyFault(self, newpoints, **kwargs):  # ,region,steps=10):
        steps = 10
        region = 'everywhere'
        if 'region' in kwargs:
            region = kwargs['region']
        if 'steps' in kwargs:
            steps = kwargs['steps']
        if 'segment' in self.overlap and 'overprint' in kwargs:
            print('overprint')
            newpointstemp = np.array(self.mesh.nodes[self.mesh.regions[region]], copy=True, order='C')
            self.overlap['segment'].applyFault(newpointstemp, **kwargs)
            newgx = np.zeros(self.mesh.n_nodes)
            newgy = np.zeros(self.mesh.n_nodes)
            newgz = np.zeros(self.mesh.n_nodes)
            self.mesh.update_property(self.name + '_gxold', self.mesh.properties[self.name + '_gx'])
            self.mesh.update_property(self.name + '_gyold', self.mesh.properties[self.name + '_gy'])
            self.mesh.update_property(self.name + '_gzold', self.mesh.properties[self.name + '_gz'])
            newgx[:] = np.nan
            newgy[:] = np.nan
            newgz[:] = np.nan

            newgx[self.mesh.regions[region]] = self.mesh.eval_interpolant(newpointstemp, self.name + '_gx', k=100)
            newgy[self.mesh.regions[region]] = self.mesh.eval_interpolant(newpointstemp, self.name + '_gy', k=100)
            newgz[self.mesh.regions[region]] = self.mesh.eval_interpolant(newpointstemp, self.name + '_gz', k=100)

            self.mesh.update_property(self.name + '_gx', newgx)
            self.mesh.update_property(self.name + '_gy', newgy)
            self.mesh.update_property(self.name + '_gz', newgz)
            self.buildFaultDisplacementField()
            # mesh2.update_property('gx1',mesh.eval_interpolant(mesh2.nodes,'gx1',k=100))
            boundaries = np.zeros(self.mesh.properties[self.name + '_gx'].shape)
            boundaries[:] = 0

            for i, e in enumerate(self.mesh.elements):
                if np.all(self.mesh.properties[self.name + '_gx'][e] > 0.):
                    continue
                if np.all(self.mesh.properties[self.name + '_gx'][e] < 0.):
                    continue
                if np.any(np.isnan(self.mesh.properties[self.name + '_gx'][e])):
                    continue
                #                 c = 0
                #                 self.mesh.property_gradients[self.name+'_gy'][i][:] = np.array([0.,0.,0.])
                #                 for n in self.mesh.neighbours[i]:
                #                     if np.all(self.mesh.properties[self.name+'_gx'][self.mesh.elements[n]] > 0.):
                #                         self.mesh.property_gradients[self.name+'_gy'][i] += self.mesh.property_gradients[self.name+'_gy'][n]
                #                         c+=1
                #                 #        print(c)
                #                 self.mesh.property_gradients[self.name+'_gy'][i] /=c
                for n in e:
                    # print(mesh2.properties['gx1'][n])
                    if self.mesh.properties[self.name + '_gx'][n] > 0:
                        boundaries[n] = 1
            self.mesh.properties['bound'] = boundaries
            self.mesh.save()
        for i in range(steps):
            g = self.mesh.eval_gradient(newpoints, self.interpolators['gy'].propertyname, k=500)
            g /= np.linalg.norm(g, axis=1)[:, None]
            g *= (1. / steps) * self.mesh.properties[self.name + '_d'][self.mesh.regions[region], None]
            g[np.any(np.isnan(g), axis=1)] = np.zeros(3)
            newpoints = newpoints + g
        return newpoints
