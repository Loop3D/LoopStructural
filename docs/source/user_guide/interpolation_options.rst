Implicit Interpolation
======================
The implicit functions have no known analytical solution which means that they need to be approximated from the observations that are provided.
The implicit function is approximated using a weighted combination of basis functions: 

.. math:: f(x,y,z) = \sum^N_{i=0} v_i \varphi_i(X) 

LoopStructural has two types of interpolation algorithms that can be used.
1. Discrete interpolation where the implicit function is approximated by basis functions on a regular grid
2. Data supported interpolation where the implicit function is approximated using basis functions located on the data points

Discrete Interpolation
-----------------------
There are two approaches that can be used for discrete interpolation. 
The first approach uses a linear function on a tetrahedral mesh. 

Piecewise Linear Interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

A regularisation constraint is used to solve the implicit system. 
The constant gradient regularisation minimises the difference in the gradient between neighbouring tetrahdron.
The regularisation constraint was first presented by Tobias Frank et al., 2007 in Computers and Geoscience. 
The same method has been extended by Guillaume Caumon et al., in 2013. 

The interpolant is assembled by combining all of the observations of the gradient, 

Additional parameters can be specified to the interpolator including:
  .. list-table:: Surfe parameters
      :widths: 25 75
      :header-rows: 1

      * - Parameter
        - Options
      * - cpw
        - 
      * - npw
        - 
      * - gpw
        - 
      * - ipw
        - integer
      * - regularisation
        - 
      * - operators
        - 
        


Finite Difference Interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Instead of using the P1 finite elements we represent the property using trilinear interpolation within a cubic element.
The property is interpolated using the 8 vertices of a cube and a linear interpolation along each axis. 

Finite differences are usedthe trilinear 

Additional parameters can be specified to the interpolator including:
  .. list-table:: Surfe parameters
      :widths: 25 75
      :header-rows: 1

      * - Parameter
        - Options
      * - cpw
        - 
      * - npw
        - 
      * - gpw
        - 
      * - ipw
        - integer
      * - regularisation
        - 
      * - operators
        - 
        
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
