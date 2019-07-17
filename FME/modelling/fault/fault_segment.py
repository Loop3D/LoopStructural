import numpy as np

from FME.modelling.structural_frame import StructuralFrame


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
        if 'name' in kwargs:
            self.name = kwargs['name']
        if 'faultfunction' in kwargs:
            self.faultfunction = kwargs['faultfunction']
        if 'displacement' in kwargs:
            self.displacement = kwargs['displacement']
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
        hw = support.get_connected_nodes_for_mask(self.faultframe.gx() >= 0)
        # get all of the points for the support that are on the foot wall or part
        # of a support element (tetra/cube) that is on the foot wall
        fw = support.get_connected_nodes_for_mask(self.faultframe.gx() <= 0)
        self.d_fw = np.zeros(fw.shape)
        self.d_hw = np.zeros(hw.shape)
        self.d_hw[hw] = 1.
        self.d_hw *= self.displacement

        self.d_fw[fw] = 0.
        self.d_fw *= self.displacement
        for i in range(steps):
            fw_g = self.faultframe.gy(fw,True)
            fw_g /= np.linalg.norm(fw_g, axis=1)[:, None]
            fw_g *= (1. / steps) * self.d_fw[:, None]
            fw_g[np.any(np.isnan(fw_g), axis=1)] = np.zeros(3)

            hw_g = self.faultframe.gy(hw, True)
            hw_g /= np.linalg.norm(hw_g, axis=1)[:, None]
            hw_g *= (1. / steps) * self.d_hw[:, None]
            hw_g[np.any(np.isnan(hw_g), axis=1)] = np.zeros(3)


            hw_n = hw + g
            fw_n = fw + g
        return hw_n, fw_n
    def apply_fault_to_data(self,data):
        for d in data:
            for i in range(steps):
                g = self.faultframe.gy(d.pos, True)
                g /= np.linalg.norm(g, axis=1)[:, None]
                g *= (1. / steps) * self.d[:, None]
                g[np.any(np.isnan(g), axis=1)] = np.zeros(3)
                d.pos = d.pos + g
        return
