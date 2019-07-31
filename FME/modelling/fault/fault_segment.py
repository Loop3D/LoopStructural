import numpy as np

class FaultSegment:
    """
    Class for representing a slip event of a fault
    """

    def __init__(self, faultframe,**kwargs):  # mesh,fault_event,data,name,region):
        """
        constructor
        :param faultframe: StructuralFrame object with three coordinates
        :param kwargs:
        """
        self.faultframe = faultframe
        self.displacement = 1
        self.gy_max = 9999
        self.gy_min = -9999
        self.faultfunction = None
        if 'name' in kwargs:
            self.name = kwargs['name']
        if 'faultfunction' in kwargs:
            self.faultfunction = kwargs['faultfunction']
        if 'displacement' in kwargs:
            self.displacement = kwargs['displacement']
        if 'gy_min' in kwargs:
            self.gy_min = kwargs['gy_min']
        if 'gy_max' in kwargs:
            self.gy_max = kwargs['gy_max']

    def evaluate(self, locations):
        return self.faultframe.features[0].evaluate_value(locations) > 0

    def apply_fault_to_support(self, support):  # ,region,steps=10):
        steps = 10


        # gx = self.faultframe.gx()
        # displacement = 1.
        # self.d = np.zeros(support.n_nodes)
        # if 'faultfunction' in kwargs:
        #     gy = self.faultframe.gy()
        #     gz = self.faultframe.gz()
        #     fault_function = kwargs['faultfunction']
        #     gxn, gyn, gzn, self.d = fault_function(gx, gy, gz)
        # if 'displacement' in kwargs:
        #     displacement = kwargs['displacement']
        # if 'faultfunction' not in kwargs:


        region = 'everywhere'
        #get all of the points for the support that are on the hanging wall or part
        #of a support element (tetra/cube) that is on the hanging wall
        # hw_p, hw_m = support.get_connected_nodes_for_mask()
        buffer = (self.faultframe.features[0].support.max_property_value() -
                self.faultframe.features[0].support.min_property_value())*.1
        hw_m = self.faultframe.get_values(0) >= -buffer
        gy_mask = np.logical_and(self.faultframe.get_values(1)>self.gy_min,
                                 self.faultframe.get_values(1) < self.gy_max)
        hw_p = support.mesh.nodes[hw_m]

        # get all of the points for the support that are on the foot wall or part
        # of a support element (tetra/cube) that is on the foot wall
        # fw_p, fw_m = support.get_connected_nodes_for_mask(self.faultframe.get_values(0) <= 5)
        fw_m = self.faultframe.get_values(0) <= buffer
        fw_p = support.mesh.nodes[fw_m]
        self.d_fw = np.zeros(fw_m.shape)
        self.d_hw = np.zeros(hw_m.shape)

        gx = self.faultframe.get_values(0)
        gy = self.faultframe.get_values(1)
        gz = self.faultframe.get_values(2)
        if self.faultfunction is None:
            self.d_hw[hw_m] = 1.
            self.d_fw[fw_m] = 0.
        if self.faultfunction is not None:
            self.d_hw[hw_m] = self.faultfunction.hw(gx[hw_m],gy[hw_m],gz[hw_m])
            self.d_fw[fw_m] = self.faultfunction.fw(gx[fw_m],gy[fw_m],gz[fw_m])
            # self.d_hw[np.logical_and(gx>1,gx<-1)] = 0
            self.d_hw[np.logical_and(gy>1,gy<-1)] = 0
            self.d_hw[np.logical_and(gz>1,gz<-1)] = 0
            # self.d_fw[np.logical_and(gx>1,gx<-1)] = 0
            self.d_fw[np.logical_and(gy>1,gy<-1)] = 0
            self.d_fw[np.logical_and(gz>1,gz<-1)] = 0
        self.d_hw *= self.displacement
        support.mesh.update_property('fw_d',self.d_fw)
        support.mesh.update_property('hw_d',self.d_hw)

        self.d_fw *= self.displacement
        hw_n = hw_p
        fw_n = fw_p
        for i in range(steps):
            fw_g = self.faultframe.evaluate_gradient(fw_p, 1)
            fw_g_m = np.linalg.norm(fw_g, axis=1)
            fw_g[fw_g_m > 0] /= fw_g_m[fw_g_m > 0, None]
            fw_g *= (1. / steps) * self.d_fw[fw_m, None]
            fw_g[np.any(np.isnan(fw_g), axis=1)] = np.zeros(3)

            hw_g = self.faultframe.evaluate_gradient(hw_p, 1)
            hw_g_m = np.linalg.norm(hw_g, axis=1)
            hw_g[hw_g_m > 0] /= hw_g_m[hw_g_m > 0, None]
            hw_g *= (1. / steps) * self.d_hw[hw_m, None]
            hw_g[np.any(np.isnan(hw_g), axis=1)] = np.zeros(3)
            hw_n += hw_g
            fw_n += fw_g
        return hw_n, fw_n, hw_m, fw_m

    def apply_fault_to_data(self,data):
        for d in data:
            for i in range(steps):
                g = self.faultframe.gy(d.pos, True)
                g /= np.linalg.norm(g, axis=1)[:, None]
                g *= (1. / steps) * self.d[:, None]
                g[np.any(np.isnan(g), axis=1)] = np.zeros(3)
                d.pos = d.pos + g
        return
