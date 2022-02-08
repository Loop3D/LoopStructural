About
======
LoopStructural is the 3D geological modelling library for Loop (Loop3d.org). 
The development of LoopStructural is lead by Lachlan Grose as an ARC (LP170100985) post-doc at Monash University. 
Laurent Ailleres and Gautier Laurent have made significant contributions to the conceptual design and integration of geological concepts into the geological modelling workflow. 
Roy Thomson and Yohan de Rose have contributed to the implementation and integration of LoopStructural into the Loop workflow.

Background of Implicit Surface Modelling
-----------------------------------------
The interpolation algorithms behind LoopStructural are implementations of research and published work from the RING team (formerly Gocad research group), Monah University and external libraries written by Michael Hillier.
Within LoopStructural we have implemented three discrete interpolation algorithms for implicit 3D modelling:

* A Piecewise Linear Interpolation algorithm where the interpolation is performed on a tetraheral mesh minimising the second derivative of the implicit funtion between neighbouring tetrahedron.
This interpolation algorithm is an implementation of the discrete smooth interpolation algorithm included in Gocad-Skua, which is heavily based on the work of Frank et al., 2007 and Caumon et al., 2007.

* This interpolation algorithm can also be framed by minimising the second derivatives using finite different on a regular cartesian grid which has been presented by Irakarama et al.,

* Within the Piecewise Linear framework we have also implemented additional constraints for modelling the geometry of folded surfaces. To do this we build a structural frame, where the structural frame characterises the axial foliation and fold axis of the fold. We modify the regularisation constraint for the folded surfaces so that the regularisation only occurs orthogonal to the fold axis and axial surface.We also use the geometry of the fold looking down plunge to add in additional constraints on the folded surface geometry. The fold constraints were first presented by Laurent et al., 2016 and the characterisation of the geometry from datasets was introduced by Grose et al., 2017


References
----------
Caumon, G., Antoine, C. and Tertois, A.: Building 3D geological surfaces from field data using implicit surfaces Field data and the need for interpretation, , 1–6, 2007.
Frank, T., Tertois, A.-L. L. and Mallet, J.-L. L.: 3D-reconstruction of complex geological interfaces from irregularly distributed and noisy point data, Comput. Geosci., 33(7), 932–943, doi:10.1016/j.cageo.2006.11.014, 2007.
Grose, L., Laurent, G., Aillères, L., Armit, R., Jessell, M. and Caumon, G.: Structural data constraints for implicit modeling of folds, J. Struct. Geol., 104, 80–92, doi:10.1016/j.jsg.2017.09.013, 2017.
Hillier, M. J., Schetselaar, E. M., de Kemp, E. A. and Perron, G.: Three-Dimensional Modelling of Geological Surfaces Using Generalized Interpolation with Radial Basis Functions, Math. Geosci., 46(8), 931–953, doi:10.1007/s11004-014-9540-3, 2014.
Laurent, G., Ailleres, L., Grose, L., Caumon, G., Jessell, M. and Armit, R.: Implicit modeling of folds and overprinting deformation, Earth Planet. Sci. Lett., 456, 26–38, doi:10.1016/j.epsl.2016.09.040, 2016.
Mallet, J.-L.: Elements of Mathematical Sedimentary Geology: the GeoChron Model, , 1–4, doi:10.3997/9789073834811, 2014.

