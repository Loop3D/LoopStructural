.. _interpolation_options:

Implicit Interpolators
======================
The implicit functions have no known analytical solution which means that they need to be approximated from the observations that are provided.
The implicit function is approximated using a weighted combination of basis functions: 

.. math:: f(x,y,z) = \sum^N_{i=0} v_i \varphi_i(X) 

Geological observations including the location of contacts and orientation of stratigraphy can be converted into mathematical expressions:

* Observations constraining the value of the scalar field 

.. math:: f(x,y,z) = v

* Observations constraining an form surface or interface between stratigraphy 

.. math:: \sum^N_{i=0} \sum^N_{j=i} f(x_i,y_i,z_i) - f(x_j,y_j,z_j) = 0

* Observations constraining the magnitude of the gradient norm of the scalar field   

.. math:: \nabla f(x,y,z) = \textbf{n}

* Observations constraining a vector to be orthogonal to the scalar field

.. math:: \nabla f(x,y,z) \cdot \textbf{t} = 0

**To solve the implicit system there needs to be a minimum of two different value observations, or a value observation and a gradient normal constraint. **

LoopStructural has two types of interpolation algorithms that can be used.
1. Discrete interpolation where the implicit function is approximated by basis functions on a predefined support
2. Data supported interpolation where the implicit function is approximated using basis functions located on the data points

Discrete Interpolation
-----------------------
There are two approaches that can be used for discrete interpolation. 
The first approach, **Piecewise Linear Interpolation**, uses a linear function on a tetrahedral mesh, in LoopStructural a tetrahedral mesh is generated from a catesian grid.
This means that the tetrahdral mesh is limited to object geometries that can be stored on a regular grid.

The basis of the piece wise linear interpolation algorithm is the linear tetrahedron where the property within the tetrahdron is interpolated using a linear function:

 .. math:: \phi (x,y,z) = a + bx + cy + dz

 This can be expressed by the values at the nodes (0-3):
 
 .. math::
 
    \begin{split}
    \phi_0 = a + bx_0 + cy_0 + dz_0 \\
    \phi_1 = a + bx_1 + cy_1 + dz_1 \\
    \phi_2 = a + bx_2 + cy_2 + dz_2 \\
    \phi_3 = a + bx_3 + cy_3 + dz_3 \\
    \end{split}

The shape functions can be used to incorporate geological observations into the approximation by determining the tetrahedron that a data point is inside.
The constraint type can then be represented as a linear equation for these vertices and the shape parameter of the tetrahedron.

A regularisation constraint is used to solve the implicit system. 
The constant gradient regularisation minimises the difference in the gradient between neighbouring tetrahdron.
The regularisation constraint was first presented by Tobias Frank et al., 2007 in Computers and Geoscience. 
The same method has been extended by Guillaume Caumon et al., in 2013. 
The regularisation term is added by adding the following constraint for every pair of neighbouring tetrahedrons in the mesh.

.. math:: \nabla\phi_{T1} - \nabla\phi_{T2} = 0

Additional parameters can be specified to the interpolator including:
  .. list-table:: Piecewise Linear Interpolator (PLI)
      :widths: 25 75
      :header-rows: 1

      * - Parameter
        - Options
      * - cpw
        - Weighting of the control points default = 1.0
      * - npw
        - Weighting of gradient norm control points default = 1.0 
      * - gpw
        - Weighting of gradient control points default = 1.0
      * - ipw
        - weighting of interface constraints default = 1.0 
      * - regularisation
        - weighting of the regularisation constraint default = 1.0
      * - cgw
        - weighting of the constant gradient regularisation
        


Finite Difference Interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Instead of using the P1 finite elements we represent the property using trilinear interpolation within a cubic element.
The property is interpolated using the 8 vertices of a cube and a linear interpolation along each axis. 

The local coordinates are determined by finding the relative location of the point within a cell:

 .. math:: \xi, \eta, \zeta 


.. math:: 
    \begin{split}
    N_1 = \frac{1}{8}(1-\xi)(1-\eta)(1-\zeta) \\
    N_2 = \frac{1}{8}(1+\xi)(1-\eta)(1-\zeta) \\
    N_3 = \frac{1}{8}(1+\xi)(1+\eta)(1-\zeta) \\
    N_4 = \frac{1}{8}(1-\xi)(1+\eta)(1-\zeta) \\
    N_5 = \frac{1}{8}(1-\xi)(1-\eta)(1+\zeta) \\
    N_6 = \frac{1}{8}(1+\xi)(1-\eta)(1-\zeta) \\
    N_7 = \frac{1}{8}(1+\xi)(1+\eta)(1+\zeta) \\
    N_8 = \frac{1}{8}(1-\xi)(1+\eta)(1+\zeta) \\
    \end{split}

We use the regularisation constraints defined by Modest Ikarama which minimises the second derivative of the implicit function.

.. math::
    \frac{\partial^2}{\partial_{xx}}+\frac{\partial^2}{\partial_{yy}}+\frac{\partial^2}{\partial_{zz}}+2\frac{\partial^2}{\partial_{xz}}+2\frac{\partial^2}{\partial_{xy}}+2 \frac{\partial^2}{\partial_{zy}} = 0

Additional parameters can be specified to the interpolator including:
  .. list-table:: Finite Difference Interpolator (FDI)
      :widths: 25 75
      :header-rows: 1

      * - Parameter
        - Options
      * - cpw
        - Weighting of the control points default = 1.0
      * - npw
        - Weighting of gradient norm control points default = 1.0 
      * - gpw
        - Weighting of gradient control points default = 1.0
      * - ipw
        - weighting of interface constraints default = 1.0 
      * - regularisation
        - weighting of the regularisation constraint default = 1.0
      * - operators
        - a dictionary of numpy arrays that can be used as masks for finite difference approximation
        
Data Supported Interpolation
-----------------------------
LoopStructural provides a wrapper to the SurfE c++ library developed  Natural Resources Canada (Geological Survey of Canada) by Michael Hillier, Eric de Kemp, and Ernst Schetselaar for the purposes of 3D structural geological modelling particularly in sparse data environments.
SurfE can be used for interpolating a GeologicalFeature in a GeologicalModel by specifying the parameter
.. code-block::
    interpolatortype = 'surfe'
    
Additional parameters can be specified to the interpolator including:
  .. list-table:: Surfe parameters
      :widths: 25 75
      :header-rows: 1

      * - Parameter
        - Options
      * - method
        - single_surface, Laujaunie, horizons
      * - kernel
        - r3, 
      * - greedy
        - tuple (interface misfit, angular misfit)
      * - poly_order
        - integer
      * - radius
        - double, radius for kernel if it uses a shape parameter
      * - anisotropy
        - boolean, whether to use global anisotropy
      
:code`interpolatoron of scalar distance fields and potential fields by choosing either the provides three methods for surface estimation.
can be installed by running the written by Michael Hillier from the 
SurfE
