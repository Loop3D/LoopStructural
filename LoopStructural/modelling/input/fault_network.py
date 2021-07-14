class FaultNetwork:
    def __init__(self,faults):
        """A fault network is a basic graph structure that 
        can return the faults for building a geological model

        Parameters
        ----------
        faults : list
            list of fault names
        """
        self.faults = faults
        self.fault_edge_count = np.zeros(len(faults),dtype=int)
        self.fault_edges = dict(zip(faults,np.arange(len(faults),dtype=int)))

        # connections 
        self.connections = {}

    def add_connection(self,fault1,fault2):
        """fault 1 is younger than fault2

        Parameters
        ----------
        fault1 : string
            name of younger fault
        fault2 : string
            name of older fault
        """
        self.connections[fault1] = fault2
        self.fault_edge_count[self.fault_edges[fault1]] +=1
        self.fault_edge_count[self.fault_edges[fault2]] +=1

    def get_fault_iterators(self):
        """
        Returns
        -------
        iterators : list
            list of fault iterators
        """
        fault_idxs = np.where(self.fault_edge_count < 2)
        return [fault for fault in faults[fault_idx] FaultNetworkIter(fault)]
        
class FaultNetworkIter:
    def __init__(self,faultname,fault_network):
        self.faultname = faultname
        self.fault_network = fault_network
    def __next__(self):
        if faultname in self.fault_network.connections:
            return self.fault_network.connections[faultname]
        else:
            return None