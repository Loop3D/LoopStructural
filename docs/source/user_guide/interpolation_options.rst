.. _interpolation_options:

Implicit Interpolators
======================
The aim of implicit modelling is to find a function which represents the distance away from a reference geological horizon.
There is no known analytical solution for this problem which means that they need to be approximated from the observations that are provided.
We can approximate the implicit function using a weighted combination of basis functions: 

.. math:: f(x,y,z) = \sum^N_{i=0} v_i \varphi_i(X) 

Geological observations including the location of contacts and orientation of stratigraphy can be converted into mathematical expressions:

* Observations constraining the value of the scalar field 

.. math:: f(x,y,z) = v

* Observations constraining an form surface or interface between stratigraphic units 

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
        
Solving discrete system
-----------------------

The discrete interpolation problems are an over determined system of equations:

.. math::
  A \cdot x = b

Where A is a rectangular sparse matrix and a row of A represents the nodes in the discrete support. 
A is over determined.
We are looking to find the solution vector x, this can be done by using least squares:

.. math::
  A.T \cdot A \cdot x = A.T \cdot b

There are many different algorithms that can be used for solving this problem and with different advantages and use cases. 
LoopStructural can be used with many different solvers and can be used with custom solvers if required.
There are two main families of solvers that are available for solving sparse linear algeba: direct and iterative. 
Direct solvers usually try to find the inverse or pseudo inverse of the matrix. 

Direct solvers used in LoopStructural are:

* lu decomposition 'using the scipy implementation'<https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.SuperLU.html>
* cholesky decomposition using sksparse library *only on linux*




Iterative solvers:

* conjugate gradient solver - using scipy
* algorithmic multigrid solver - using pyamg
* lsqr - this solver uses the rectangular matrix directly, therefore does not require computing A.T A and A.T B. Uses scipy

.. code-block:
  model.create_and_add_foliation('my_foliation',solver='lu')
  model.create_and_add_foliation('my_foliation',solver='chol')
  model.create_and_add_foliation('my_foliation',solver='cg')
  model.create_and_add_foliation('my_foliation',solver='lsqr')
  model.create_and_add_foliation('my_foliation',solver='pyamg')

Using an external solver:

You can also pass a function that solves

.. math::
   A.T \cdot A \cdot x = A.T \cdot B

to LoopStructural if you want to use another solver by using the `external` keyword . 

.. code-block:
  def mysolver(A,B):
    from scipy.sparse.linalg import gmres
    x = gmres(A,B)
    return x[0]
  model.create_and_add_foliation('my_foliation',solver='external',external=mysolver)
 
 The solution to the least squares problem will be stored in the interpolator object and can be easily accessed:

.. code-block:
  model.create_and_add_foliation('my_foliation',solver='pyamg')
  pyamg_solution = model['my_foliation'].interpolator.c[model['my_foliation'].interpolator.region]

Note that when interpolating a subset of a mesh using a region, LoopStructural will fill the unused nodes in the interpolation support with nan.
To return only the values that are related to the solver make sure you mask using the region.  

The choice of solver is somewhat dependent on the model you are creating. 
If you have a small model and a relatively large amount of memory on your computer a direct solver may be the most appropriate.

If your model is large, or computer memory is a limitation an iterative solver is the best choice.
However, iterative solvers can suffer from poor convergence if the matrix is poorly conditioned (this is possible when modeling irregular geometries).
Solving the iterative problem can be sped up by using a preconditioner for the matrix. 

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
      

