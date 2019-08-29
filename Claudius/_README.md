# Implicit benchmark data set - Claudius carbonate build-ups, Carnarvon Basin, W Australia. 

*Data set interpreted and prepared in July-August 2019 by Guillaume.Caumon@univ-lorraine.fr, 
to whom questions should be addressed.*

This data set includes a Mesozoic series overlying Upper Triasic carbonate
build-ups in the Carnarvon Basin, offshore Western Australia. The data were
picked on the Claudius3D seismic survey acquired by Westerngeco, courtesy of
Geoscience Australia. The Z coordinate was obtained by multiplying the two-way time by -3. 

## Contents

This data set contains four horizons and a fault. The picked horizon points were
decimated both randomly and manually so as to generate data gaps. All horizons
are sampled by points and by one polygonal line each, so that 2D implicit 
implementations can also be tested. 
 * ClaudiusSeismic.png is a screenshot of the seismic section and of the four
   interpreted horizons, for loading QC. North correponds to increasing X, East to
   increasing Y. Z is positive upwards. 
 * [ABCD]Section.csv correspond to the interpreted horizons A, B, C and D on the seismic
   crossline 2399. 
 * Dips.csv are points with orientation data (gradient vector), picked from a
   few seismic locations. 
 * [ABCD]Points.csv are point clouds of varying density sampling the four
   horizons A, B, C and D.  
 * Fault.[ts,dxf]: the interpreted fault surface affecting Horizon D, in Gocad TSurf 
   and AutoCAD DXF formats.

## Output

Please rename the `Results_` as `Results_YOURNAME` directory. 
Please include any relevant scripts and descriptions to explain the parameters
used for the implicit method (including software name and version).
 
In case you need them, the directory already contains 
 * Box.vo: An empty Voxet (Cartesian Grid) in Gocad ASCII format
 * Solid.vo: An empty conformable tetrahedral mesh in Gocad ASCII format. 

These volumetric supports can be used to visualize (and possibly compute) 
the implicit scalar field representing the stratigraphy. 

### Scalar Field

Suggested size of the axis-aligned modeling box: 
Origin: 548800 7816600 -8400
Maximum: 552500 7822000 -11010

Suggested resolution: 100m x 100m x -90m (grid size 38 x 55 x 30)

NB: An empty Gocad Voxet with this geometry is already in the results directory. 

You can also provide the values in ASCII or binary format. 
Please just indicate the grid geometry, element size (if binary) and fast/slow axes. 


### Surfaces 

Please use DXF, Gocad and/or .obj format. 

