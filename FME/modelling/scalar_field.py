class ScalarField:
    """
    Generic support to represent a scalar field
    """
    def __init__(self, support)
        self.support = support

    def evaluate_value(self,evaluation_points):

        pass


    def evaluate_gradient(self,evaluation_points):


class TetrahedralMeshScalarField(ScalarField):
    """
    Support to represent a scalar field using a tetrahedral mesh
    """
    def __init__(self, mesh, **kwargs):
        super().__init__(self,mesh)
        self.support = mesh

        if 'propertyname' in kwargs:
            self.property_name = kwargs['propertyname']
        if 'nodevalues' in kwargs:
            self.node_vaues = kwargs['nodevalues']
            self.support.update_property(self.property_name, self.node_values)

    def evaluate_value(self,evaluation_points):
        self.support.evaluate_property(self.property_name,evaluation_points)

    def evaluate_gradient(self,evaluation_points):
        return self.support.evaluate_gradient(self.property_name, evaluation_points)
