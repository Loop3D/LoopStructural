Input data
==========
To create a geological model the data needs to be formatted in a table the following headings:

* X - x component of the cartesian coordinates
* Y - y component of the cartesian coordinates
* Z - z component of the cartesian coordinates
* feature_name - unique name of the geological feature being modelled
* val - value observations of the scalar field
* interface - unique identifier for an inteface containing similar scalar field values
* nx - x component of the gradient norm
* ny - y component of the gradient norm
* nz - z component of the gradient norm
* gx - x component of a gradient constraint
* gy - y component of a gradient constraint
* gz - z component of a gradient constraint
* tx - x component of a gradient tangent constraint
* ty - y component of a gradient tangent constraint
* tz - z component of a gradient tangent constraint
* coord - coordinate of the structural frame data point is used for

Value constraints
------------------
Value constraints set the value of the implicit function to equal the value. 
This should represent the distance from the reference horizon. 
It is important to consider the range in value of this field and the coordinates of the model.
If the range in value of the scalar field is very large relative to the model coordinates, it would be expected that the resulting model will be very steep.

Interface constraints
---------------------
Interface constraints are points which should have the same value of the implicit function.
This is comparable to the way data are interpreted using the Lajaunnie approach. 
A warning, when using interface constraints there is no gradient information implicitly defined by the value of the points.
This means that unless there are sufficient gradient normal constraints located between interfaces the scalar field may not appear as expected.
Another consideration is that using discrete approaches the implicit function is solved in a least squares sense, this means that all of the constraints minimised equally.
Due to the implementation of interface constraints, these tend to be weighted higher in the system of equations due to there being more pairs of points than points.
As a result, it is recommended to increase the weighting of the normal constraints by at least the average number of points in an interface.

Gradient norm constraints
-------------------------
Sets the partial derivatives of the implicit function to equal the components of the gradient norm vector.
This has a similar limitation as the value constraints, the magnitude of the vector needs to within the correct scale of the model.
This means the normal vector should be scaled in the same way as the model coordinates. 

Gradient constraints
--------------------
Gradient constraints find a pair of vectors that are orthogonal to the normal vector to the field. 
This this can be the strike vector and the dip vector.
Two constraints are added into the implicit functions constraining the scalar field to be orthogonal to these vectors.
This constraint will control the gradient of the scalar field without controlling the polarity or the length of the gradient norm.
For this reason these constraints are not as strong as the gradient norm constraint.


Tangent constraints
-------------------
Constrains a vector to be orthogonal to the scalar field by adding the following constraint:

.. math:: f(x,y,z) \cdot \textbf{t} = 0

This constraint does not control the vaue of the scalar field and is quite weak. 
It is not recommended for use, however it could be used for integrating form lines, instead of using interface constraints.
