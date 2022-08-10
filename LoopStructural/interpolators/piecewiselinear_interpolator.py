"""
Piecewise linear interpolator
"""
import numpy as np
from LoopStructural.utils import getLogger

logger = getLogger(__name__)
try:
    from LoopStructural.interpolators._cython.dsi_helper import cg, fold_cg
except:
    from LoopStructural.interpolators._python.dsi_helper import cg, fold_cg

    logger.warning("Using python code not compiled cython code, this will be slower")
from LoopStructural.interpolators import DiscreteInterpolator
from LoopStructural.interpolators import InterpolatorType

from LoopStructural.utils.helper import get_vectors


class PiecewiseLinearInterpolator(DiscreteInterpolator):
    """ """

    def __init__(self, support):
        """
        Piecewise Linear Interpolator
        Approximates scalar field by finding coefficients to a piecewise linear
        equation on a tetrahedral mesh. Uses constant gradient regularisation.

        Parameters
        ----------
        support - TetMesh
            interpolation support
        """

        self.shape = "rectangular"
        DiscreteInterpolator.__init__(self, support)
        # whether to assemble a rectangular matrix or a square matrix
        self.interpolator_type = InterpolatorType.PIECEWISE_LINEAR
        self.support = support

        self.interpolation_weights = {
            "cgw": 0.1,
            "cpw": 1.0,
            "npw": 1.0,
            "gpw": 1.0,
            "tpw": 1.0,
            "ipw": 1.0,
        }
        self.__str = "Piecewise Linear Interpolator with %i unknowns. \n" % self.nx
        self.type = InterpolatorType.PIECEWISE_LINEAR

    def __str__(self):
        return self.__str

    def copy(self):
        return PiecewiseLinearInterpolator(self.support)

    def _setup_interpolator(self, **kwargs):
        """
        Searches through kwargs for any interpolation weights and updates
        the dictionary.
        Then adds the constraints to the linear system using the
        interpolation weights values
        Parameters
        ----------
        kwargs -
            interpolation weights

        Returns
        -------

        """
        # can't reset here, clears fold constraints
        # self.reset()
        logger.info("Setting up PLI interpolator for %s" % self.propertyname)
        for key in kwargs:
            if "regularisation" in kwargs:
                self.interpolation_weights["cgw"] = 0.1 * kwargs["regularisation"]
            self.up_to_date = False
            self.interpolation_weights[key] = kwargs[key]
        if self.interpolation_weights["cgw"] > 0.0:
            self.up_to_date = False
            self.add_constant_gradient(
                self.interpolation_weights["cgw"],
                direction_feature=kwargs.get("direction_feature", None),
                direction_vector=kwargs.get("direction_vector", None),
            )
            logger.info(
                "Using constant gradient regularisation w = %f"
                % self.interpolation_weights["cgw"]
            )
        logger.info(
            "Added %i gradient constraints, %i normal constraints,"
            "%i tangent constraints and %i value constraints"
            "to %s" % (self.n_g, self.n_n, self.n_t, self.n_i, self.propertyname)
        )
        self.add_gradient_constraints(self.interpolation_weights["gpw"])
        self.add_norm_constraints(self.interpolation_weights["npw"])
        self.add_value_constraints(self.interpolation_weights["cpw"])
        self.add_tangent_constraints(self.interpolation_weights["tpw"])
        self.add_interface_constraints(self.interpolation_weights["ipw"])

    def add_constant_gradient(
        self, w=0.1, direction_vector=None, direction_feature=None
    ):
        """
        Add the constant gradient regularisation to the system

        Parameters
        ----------
        w (double) - weighting of the cg parameter

        Returns
        -------

        """
        if direction_feature is not None:
            direction_vector = direction_feature.evaluate_gradient(
                self.support.barycentre
            )
        if direction_vector is not None:
            if direction_vector.shape[0] == 1:
                # if using a constant direction, tile array so it works for cg calc
                direction_vector = np.tile(
                    direction_vector, (self.support.barycentre.shape[0], 1)
                )
        if direction_vector is not None:
            logger.info("Running constant gradient")
            elements_gradients = self.support.get_element_gradients(
                np.arange(self.ntetra)
            )
            if elements_gradients.shape[0] != direction_vector.shape[0]:
                logger.error(
                    "Cannot add directional CG, vector field is not the correct length"
                )
                return
            region = self.region.astype("int64")

            neighbours = self.support.get_neighbours()
            elements = self.support.get_elements()
            idc, c, ncons = fold_cg(
                elements_gradients,
                direction_vector,
                neighbours.astype("int64"),
                elements.astype("int64"),
                self.nodes,
            )

            idc = np.array(idc[:ncons, :])
            A = np.array(c[:ncons, :])
            B = np.zeros(c.shape[0])
            gi = np.zeros(self.support.n_nodes)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx)
            idc = gi[idc]
            outside = ~np.any(idc == -1, axis=1)

            # w/=A.shape[0]
            self.add_constraints_to_least_squares(
                A[outside, :],
                B[outside],
                idc[outside, :],
                w=w,
                name="direction_regularisation",
            )
        else:
            logger.info("Running constant gradient")
            elements_gradients = self.support.get_element_gradients(
                np.arange(self.support.ntetra)
            )
            region = self.region.astype("int64")

            neighbours = self.support.get_neighbours()
            elements = self.support.get_elements()
            idc, c, ncons, area = cg(
                elements_gradients,
                neighbours.astype("int64"),
                elements.astype("int64"),
                self.support.nodes,
                region.astype("int64"),
            )

            idc = np.array(idc[:ncons, :])
            A = np.array(c[:ncons, :])
            area = np.array(area[:ncons])
            B = np.zeros(A.shape[0])
            gi = np.zeros(self.support.n_nodes)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx)
            idc = gi[idc]
            outside = ~np.any(idc == -1, axis=1)
            # w/=A.shape[0]
            self.add_constraints_to_least_squares(
                A[outside, :], B[outside], idc[outside, :], w=w, name="regularisation"
            )

    def add_gradient_constraints(self, w=1.0):
        """
        Adds gradient constraints to the least squares system with a weight
        defined by w
        Parameters
        ----------
        w : either numpy array of length number of

        Returns
        -------
        Notes
        -----
        Gradient constraints add a constraint that the gradient of the
        implicit function should
        be orthogonal to the strike vector and the dip vector defined by the
        normal.
        This does not control the direction of the gradient and therefore
        requires at least two other
        value constraints OR a norm constraint for the interpolant to solve.
        """

        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            (
                vertices,
                element_gradients,
                tetras,
                inside,
            ) = self.support.get_element_gradient_for_location(points[:, :3])
            # e, inside = self.support.elements_for_array(points[:, :3])
            # nodes = self.support.nodes[self.support.elements[e]]
            vecs = vertices[:, 1:, :] - vertices[:, 0, None, :]
            vol = np.abs(np.linalg.det(vecs)) / 6
            norm = np.linalg.norm(points[:, 3:6], axis=1)
            points[:, 3:6] /= norm[:, None]
            element_gradients /= norm[:, None, None]
            # d_t *= vol[:,None,None]
            strike_vector, dip_vector = get_vectors(points[:, 3:6])
            A = np.einsum("ji,ijk->ik", strike_vector, element_gradients)

            # A *= vol[:, None]

            gi = np.zeros(self.support.n_nodes).astype(int)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx).astype(int)
            # w /= 3
            idc = gi[tetras]
            B = np.zeros(idc.shape[0])
            outside = ~np.any(idc == -1, axis=1)
            w *= points[:, 6]
            self.add_constraints_to_least_squares(
                A[outside, :], B[outside], idc[outside, :], w=w, name="gradient strike"
            )
            A = np.einsum("ji,ijk->ik", dip_vector, element_gradients)
            # A *= vol[:, None]
            self.add_constraints_to_least_squares(
                A[outside, :], B[outside], idc[outside, :], w=w, name="gradient dip"
            )

    def add_norm_constraints(self, w=1.0):
        """
        Extracts the norm vectors from the interpolators p_n list and adds
        these to the implicit
        system

        Parameters
        ----------
        w : double
            weighting of the norm constraints in a least squares system

        Returns
        -------
        Notes
        -----
        Controls the direction and magnitude of the norm of the scalar field
        gradient.
        This constraint can conflict with value constraints if the magnitude
        of the vector doesn't
        match with the value constraints added to the implicit system.
        """

        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            (
                vertices,
                element_gradients,
                tetras,
                inside,
            ) = self.support.get_element_gradient_for_location(points[:, :3])
            # e, inside = self.support.elements_for_array(points[:, :3])
            # nodes = self.support.nodes[self.support.elements[e]]
            vol = np.zeros(element_gradients.shape[0])
            vecs = vertices[:, 1:, :] - vertices[:, 0, None, :]
            vol = np.abs(np.linalg.det(vecs)) / 6
            d_t = element_gradients
            d_t[inside, :, :] *= vol[inside, None, None]
            # add in the element gradient matrix into the inte
            idc = np.tile(tetras[:, :], (3, 1, 1))
            idc = idc.swapaxes(0, 1)
            # idc = self.support.elements[e]
            gi = np.zeros(self.support.n_nodes).astype(int)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx).astype(int)
            idc = gi[idc]
            outside = ~np.any(idc == -1, axis=2)
            outside = outside[:, 0]
            w = points[:, 6] * w
            points[inside, 3:6] *= vol[inside, None]
            # w /= 3

            self.add_constraints_to_least_squares(
                d_t[outside, :, :], points[outside, 3:6], idc[outside], w=w, name="norm"
            )

    def add_value_constraints(self, w=1.0):  # for now weight all value points the same
        """
        Adds value constraints to the least squares system

        Parameters
        ----------
        w : double
            weight

        Returns
        -------

        """

        # get elements for points
        points = self.get_value_constraints()
        if points.shape[0] > 1:
            vertices, c, tetras, inside = self.support.get_element_for_location(
                points[:, :3]
            )
            # calculate volume of tetras
            vecs = vertices[inside, 1:, :] - vertices[inside, 0, None, :]
            vol = np.abs(np.linalg.det(vecs)) / 6
            A = c[inside]
            A *= vol[:, None]
            idc = tetras[inside, :]
            # now map the index from global to region create array size of mesh
            # initialise as np.nan, then map points inside region to 0->nx
            gi = np.zeros(self.support.n_nodes).astype(int)
            gi[:] = -1

            gi[self.region] = np.arange(0, self.nx)
            idc = gi[idc]
            outside = ~np.any(idc == -1, axis=1)
            w *= points[inside, 4]
            self.add_constraints_to_least_squares(
                A[outside, :],
                points[inside, :][outside, 3] * vol[:],
                idc[outside, :],
                w=w,
                name="value",
            )

    def add_interface_constraints(
        self, w=1.0
    ):  # for now weight all value points the same
        """
        Adds a constraint that defines all points with the same 'id' to be the same value
        Sets all P1-P2 = 0 for all pairs of points

        Parameters
        ----------
        w : double
            weight

        Returns
        -------

        """
        points = self.get_interface_constraints()
        if points.shape[0] > 1:
            vertices, c, tetras, inside = self.support.get_element_for_location(
                points[:, :3]
            )
            # calculate volume of tetras
            vecs = vertices[inside, 1:, :] - vertices[inside, 0, None, :]
            vol = np.abs(np.linalg.det(vecs)) / 6
            A = c[inside]
            A *= vol[:, None]
            idc = tetras[inside, :]
            for unique_id in np.unique(
                points[np.logical_and(~np.isnan(points[:, 3]), inside), 3]
            ):
                mask = points[inside, 3] == unique_id
                ij = np.array(
                    np.meshgrid(
                        np.arange(0, A[mask, :].shape[0]),
                        np.arange(0, A[mask, :].shape[0]),
                    )
                ).T.reshape(-1, 2)
                interface_A = np.hstack(
                    [A[mask, :][ij[:, 0], :], -A[mask, :][ij[:, 1], :]]
                )
                # interface_A = interface_A.reshape((interface_A.shape[0]*interface_A.shape[0],A.shape[1]))
                interface_idc = np.hstack(
                    [idc[mask, :][ij[:, 0], :], idc[mask, :][ij[:, 1], :]]
                )

                # now map the index from global to region create array size of mesh
                # initialise as np.nan, then map points inside region to 0->nx
                gi = np.zeros(self.support.n_nodes).astype(int)
                gi[:] = -1

                gi[self.region] = np.arange(0, self.nx)
                interface_idc = gi[interface_idc]
                # interface_idc = np.tile(interface_idc,(interface_idc.shape[0],1)).reshape(interface_A.shape)#flatten()
                outside = ~np.any(interface_idc == -1, axis=1)
                self.add_constraints_to_least_squares(
                    interface_A[outside, :],
                    np.zeros(interface_A[outside, :].shape[0]),
                    interface_idc[outside, :],
                    w=w,
                    name="interface_{}".format(unique_id),
                )

    def add_gradient_orthogonal_constraints(self, points, vector, w=1.0, B=0):
        """
        constraints scalar field to be orthogonal to a given vector

        Parameters
        ----------
        position
        normals
        w
        B

        Returns
        -------

        """
        if points.shape[0] > 0:
            (
                vertices,
                element_gradients,
                tetras,
                inside,
            ) = self.support.get_element_gradient_for_location(points[:, :3])
            # e, inside = self.support.elements_for_array(points[:, :3])
            # nodes = self.support.nodes[self.support.elements[e]]
            norm = np.linalg.norm(vector, axis=1)
            mask = np.isnan(norm)
            mask = ~np.logical_or(mask, norm == 0)
            vector[mask, :] /= norm[mask, None]
            vecs = vertices[:, 1:, :] - vertices[:, 0, None, :]
            # vol = np.abs(np.linalg.det(vecs)) / 6
            element_gradients[mask, :, :] /= norm[mask, None, None]

            A = np.einsum("ij,ijk->ik", vector, element_gradients)

            # A *= vol[:, None]

            gi = np.zeros(self.support.n_nodes).astype(int)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx).astype(int)
            idc = gi[tetras]
            B = np.zeros(idc.shape[0]) + B
            outside = ~np.any(idc == -1, axis=1)
            self.add_constraints_to_least_squares(
                A[outside, :],
                B[outside],
                idc[outside, :],
                w=w,
                name="gradient orthogonal",
            )
