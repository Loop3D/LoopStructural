from .process_data import ProcessInputData

class ProcessM2LDirectory:
    def __init__(self,directory):
        groups = pd.read_csv(m2l_directory + '/tmp/all_sorts_clean.csv', index_col=0)
        orientations = pd.read_csv(m2l_directory + '/output/orientations_clean.csv')
        contacts = pd.read_csv(m2l_directory + '/output/contacts_clean.csv')
        fault_orientations = None
        fault_locations = None
        fault_dimensions = None
        try:
            fault_orientations = pd.read_csv(m2l_directory + '/output/fault_orientations.csv')
            fault_orientations.rename(columns={'formation':'fault_name'},inplace=True)
            fault_orientations['strike'] = fault_orientations['DipDirection'] + 90
        except:
            logger.error("Can't read fault orientation")
        try:
            fault_locations = pd.read_csv(m2l_directory + '/output/faults.csv')
            fault_locations.rename(columns={'formation':'fault_name'},inplace=True)
        except:
            logger.error("Can't read fault locations")
        try:
            fault_dimensions = pd.read_csv(m2l_directory + '/output/fault_dimensions.csv',index_col='Fault')
            fault_dimensions.rename(columns={'Fault':'fault_name','HorizontalRadius':'InfluenceDistance','VerticalRadius':'SlipDistance','InfluenceDistance':'ExtentDistance'},inplace=True)

        except:
            logger.error("Can't read fault dimensions")

        fault_graph = networkx.read_gml(m2l_directory + '/tmp/fault_network.gml')
        
        try:
            fault_displacements = pd.read_csv(m2l_directory + '/output/fault_displacements3.csv')
            fault_dimensions['displacement'] = np.nan
            fault_dimensions['downthrow'] = np.nan

            for fname in fault_dimensions.index:
                fault_dimensions.loc[fname,'displacement'] = fault_displacements.loc[fault_displacements['fname']==fname,
                                                                            'vertical_displacement'].median()
        except:
            logger.error("Can't read fault displacement setting fault displacement to 0 ")
            fault_dimensions['displacement'] = 0
        intrusions = []
        fault_stratigraphy = []
        stratigraphic_order = [list(groups['code'])]
        thicknesses = None
        try:
            formation_thickness = pd.read_csv(m2l_directory+'/output/formation_summary_thicknesses.csv')
            thicknesses = dict(zip(list(formation_thickness['formation']),list(formation_thickness['thickness median'])))
        except:
            logger.error("Can't open thickness file, not using thicknesses")

        ProcessInputData.__init__( 
                    contacts, 
                    orientations, 
                    stratigraphic_order,
                    thicknesses=thicknesses,
                    fault_orientations=fault_orientations,
                    fault_locations=fault_locations,
                    fault_dimensions=fault_dimensions,
                    fault_graph=fault_graph
                    )
    

