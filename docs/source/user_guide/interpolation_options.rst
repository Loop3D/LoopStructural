Implicit Interpolation
======================
LoopStructural has two types of interpolation algorithms that can be used.
1. Discrete interpolation where the implicit function is approximated on a regular grid
2. Data supported interpolation where the implicit function is approximated using kernels located on the data points

Discrete Interpolation
-----------------------
There are two approaches that can be used for discrete interpolation. 
The first approach uses a linear function on a tetrahedral mesh. 
Piecewise Linear Interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finite Difference Interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data Supported Interpolation
-----------------------------
LoopStructural provides a wrapper to the SurfE c++ library developed  Natural Resources Canada (Geological Survey of Canada) by Michael Hillier, Eric de Kemp, and Ernst Schetselaar for the purposes of 3D structural geological modelling particularly in sparse data environments.
SurfE can be used for interpolating a GeologicalFeature in a GeologicalModel by specifying the parameter
.. code-block::
    interpolatortype = 'surfe'
    
Additional parameters can be specified to the interpolator including:
:code`interpolatoron of scalar distance fields and potential fields by choosing either the provides three methods for surface estimation.
can be installed by running the written by Michael Hillier from the 
SurfE
