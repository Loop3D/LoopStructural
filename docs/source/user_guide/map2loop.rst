The class method (GeologicalModel.from_map2loop_directory(m2l_directory, **kwargs)) creates an instance of a  GeologicalModel from a root map2loop output directory. 
There are two steps to this function the first is a pre-processing step that prepares the extracted datasets for LoopStructural. This involves combining all observations (fault surface location, fault orientation, stratigraphic orientation and stratigraphic contacts) into the one pandas DataFrame. The cumulative thickness is calculated by iterating over the stratigraphic column for each group and starting with the base of the group as 0 and adding the thickness of the lower unit. This function returns a python dictionary which is structured as shown in Table 5. The processed data dictionary can then be passed to a utility function that:
1.	initialises a GeologicalModel for the specified bounding box (map area chosen for map2loop and depth specified).
2.	associates the pandas DataFrame to the model
3.	iterate through the fault-stratigraphy table
•	adding stratigraphic packages that are associated with no faults
•	add faults in reverse order, fault intersections are manually, and fault geometries are only interpolated within their object aligned bounding box to reduce computational to ensure that 
•	The base of each group is assigned as an unconformity for any older stratigraphic units or faults that occur below it
4.	associate stratigraphic column to geological model
 
A parameters dictionary can be passed to the GeologicalModel.from_map2loop_directory function which allows for the interpolator parameters (listed in Table 2) for the model to be specified. This can be passed using two python dictionaries: foliation_params contains the parameters for modelling any stratigraphic groups and fault_params contains parameters for the faults. 
