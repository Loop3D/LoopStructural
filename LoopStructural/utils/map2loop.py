import pandas as pd
import numpy as np
import os
import logging
import networkx
from scipy.stats import truncnorm

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

def process_map2loop(m2l_directory, flags={}):
    """
    Extracts relevant information from map2loop outputs

    Parameters
    ----------
    m2l_directory : string
        absolute path to a directory containing map2loop outputs
    Returns
    -------
    m2l_data : dict
        a dictionary containing the extracted and collated data
    """
    gradient = flags.get('gradient',False)
    vector_scale = flags.get('vector_scale',None)
    tangents = pd.read_csv(m2l_directory + '/tmp/raw_contacts.csv')
    groups = pd.read_csv(m2l_directory + '/tmp/all_sorts_clean.csv', index_col=0)
    contact_orientations = pd.read_csv(m2l_directory + '/output/orientations_clean.csv')
    # formation_thickness = pd.read_csv)
    contacts = pd.read_csv(m2l_directory + '/output/contacts_clean.csv')
    displacements = pd.read_csv(m2l_directory + '/output/fault_displacements3.csv')
    fault_orientations = pd.read_csv(m2l_directory + '/output/fault_orientations.csv')
    fault_locations = pd.read_csv(m2l_directory + '/output/faults.csv')
    fault_fault_relations = pd.read_csv(m2l_directory + '/output/fault-fault-relationships.csv')
    fault_strat_relations = pd.read_csv(m2l_directory + '/output/group-fault-relationships.csv')
    fault_dimensions = pd.read_csv(m2l_directory + '/output/fault_dimensions.csv')
    fault_graph = networkx.read_gml(m2l_directory + '/tmp/fault_network.gml')

    # whether to simulate thickness variations
    thickness_probabilities = flags.get('thickness_probabilities',False)
    orientation_probabilities = flags.get('orientation_probabilities',False)
    fault_orientation_probabilities = flags.get('fault_orientation_probabilities',False)
    fault_slip_vector = flags.get('fault_slip_vector',None)
    ## read supergroups file
    supergroups = {}
    sgi = 0
    contact_orientations.loc[contact_orientations['polarity']==0,'polarity'] = -1
    try:
        with open(m2l_directory + '/tmp/super_groups.csv') as    f:
            for l in f:
                for g in l.split(','):
                    g = g.replace('-','_').replace(' ','_')
                    if g.find('\n') > 0:
                        g = g[:g.find('\n')]
                    supergroups[g] = 'supergroup_{}'.format(sgi)
                    if g == '\n':
                        sgi += 1
                        break
    except:
        for g in groups['group'].unique():
            supergroups[g] = g
    try:
        supergroups.pop('\n')
    except KeyError:
        pass
    
    ## read bounding box
    bb = pd.read_csv(m2l_directory+'/tmp/bbox.csv')

    ## read fault intersection graph
    fault_intersection_angles = {}
    with open(m2l_directory+'/graph/fault-fault-intersection.txt','r') as f:
        for l in f:
            fault_id = int(l.split(',')[1])
            fault_intersection_angles['Fault_{}'.format(fault_id)] = []
            intersections = l[l.find('{')+1:l.find('}')]
            intersections = intersections.replace('(','')
            intersections = intersections.replace(')','')
            intersection_data = intersections.split(',')
            i = 0
            while i < len(intersection_data):
                fault_intersection_angles['Fault_{}'.format(fault_id)].append(['Fault_{}'.format(int(intersection_data[i])),intersection_data[i+1],float(intersection_data[i+2])])
                i+=3


    # process tangent data to be tx, ty, tz
    tangents['tz'] = 0
    tangents['tx'] = tangents['lsx']
    tangents['ty'] = tangents['lsy']
    tangents.drop(['angle', 'lsx', 'lsy'], inplace=True, axis=1)

    # convert azimuth and dip to gx, gy, gz
    

    # calculate scalar field values
    i = 0
    thickness = {}
    max_thickness = 0
    thickness_file = flags.get('thickness',m2l_directory + '/output/formation_summary_thicknesses.csv')
    with open(thickness_file) as file:
        for l in file:
            if i>=1:
                linesplit = l.split(',')
                if thickness_probabilities:
                    std = float(linesplit[2])
                    mean = float(linesplit[1])
                    if  np.isnan(std):
                        std = float(linesplit[1])
                    a = (0-mean) / std
                    b = 100.
                    thickness[linesplit[0]] = (truncnorm.rvs(a,b,size=1)*std)[0]+mean
                else:
                    thickness[linesplit[0]] = float(linesplit[1])
                
                # normalise the thicknesses
                if float(linesplit[1]) > max_thickness:
                    max_thickness=float(linesplit[1])
    #             print(l.split(',')[1])
            i+=1
    # for k in thickness.keys():
    #     thickness[k] /= max_thickness
    if vector_scale is None:
        vector_scale = max_thickness

    if gradient:
        from LoopStructural.utils.helper import strike_dip_vector
        contact_orientations['strike'] = contact_orientations['azimuth'] - 90
        contact_orientations['gx'] = np.nan
        contact_orientations['gy'] = np.nan
        contact_orientations['gz'] = np.nan
        contact_orientations[['gx', 'gy', 'gz']] = strike_dip_vector(contact_orientations['strike'],
                                                                    contact_orientations['dip'])*max_thickness   
    if np.sum(contact_orientations['polarity']==0) >0 and np.sum(contact_orientations['polarity']==-1)==0:
        # contact_orientations['polarity']+=1
        contact_orientations.loc[contact_orientations['polarity']==0]=-1
    if not gradient:
        from LoopStructural.utils.helper import strike_dip_vector
        contact_orientations['strike'] = contact_orientations['azimuth'] - 90
        contact_orientations['nx'] = np.nan
        contact_orientations['ny'] = np.nan
        contact_orientations['nz'] = np.nan
        contact_orientations[['nx', 'ny', 'nz']] = strike_dip_vector(contact_orientations['strike'],
                                                                    contact_orientations['dip'])*vector_scale *contact_orientations['polarity'].to_numpy()[:,None]
    contact_orientations.drop(['strike', 'dip', 'azimuth'], inplace=True, axis=1)
    # with open(m2l_directory + '/output/formation_summary_thicknesses.csv') as file:

    # thickness = {}
    # for f in formation_thickness['formation'].unique():
    #     thickness[f] = formation_thickness[formation_thickness['formation'] == f]['thickness'])

    strat_val = {}
    stratigraphic_column = {}
    unit_id = 0
    val = {}
    for i in groups['group number'].unique():
        g = supergroups[groups.loc[groups['group number'] == i, 'group'].iloc[0]]
        if g not in stratigraphic_column:
            stratigraphic_column[g] = {}
            val[g] = 0

        for c, colour in zip(groups.loc[groups['group number'] == i, 'code'],groups.loc[groups['group number'] == i, 'colour']):
            strat_val[c] = np.nan
            if c in thickness:
                stratigraphic_column[g][c] = {'max': val[g], 'min': val[g] - thickness[c], 'id': unit_id, 'colour':colour}
                unit_id += 1
                strat_val[c] = val[g] - thickness[c]
                val[g] -= thickness[c]
    group_name = None
    for g, i in stratigraphic_column.items():
        if len(i) ==0:
            for gr, sg in supergroups.items():
                if sg == g:
                    group_name = gr
                    break
            try:
                if group_name is None:
                    continue
                c=groups.loc[groups['group']==group_name,'code'].to_numpy()[0]
                strat_val[c] = 0
                stratigraphic_column[g] = {c:{'min':0,'max':9999,'id':unit_id}}
                unit_id+=1
                group_name = None
            except:
                print('Couldnt process {}'.format(g))
    #whether to use thickness or interface
    use_thickness = flags.get('use_thickness',True)
    if use_thickness:
        contacts['val'] = np.nan
        for o in strat_val:
            contacts.loc[contacts['formation'] == o, 'val'] = strat_val[o]
    if use_thickness == False:
        contacts['interface'] = np.nan
        interface_val = 0
        for u in contacts['formation'].unique():
            contacts.loc[contacts['formation'] == u,'interface'] = interface_val
            interface_val+=1
    tangents['feature_name'] = tangents['group']
    contact_orientations['feature_name'] = None
    contacts['feature_name'] = None
    for g in groups['group'].unique():
        val = 0
        for c in groups.loc[groups['group'] == g, 'code']:
            contact_orientations.loc[contact_orientations['formation'] == c, 'feature_name'] = supergroups[g]
            contacts.loc[contacts['formation'] == c, 'feature_name'] = supergroups[g]
    displacements['dip_dir'] = np.nan
    for fname in fault_orientations['formation'].unique():
        displacements.loc[displacements['fname'] == fname, 'dip_dir'] = np.mean(
            fault_orientations.loc[fault_orientations['formation'] == fname, 'DipDirection'])
    max_displacement = {}
    downthrow_dir = {}
    fault_locations['val'] = 0
    fault_locations['coord'] = 0
    fault_orientations['coord'] = 0
    fault_orientations['gx'] = np.nan
    fault_orientations['gy'] = np.nan
    fault_orientations['gz'] = np.nan

    stratigraphic_column['faults'] = {}
    for f in displacements['fname'].unique():
        fault_centers = np.zeros(6)
        normal_vector = np.zeros(3)
        strike_vector = np.zeros(3)
        slip_vector = np.zeros(3)

        fault_edges = np.zeros((2,3))
        fault_tips = np.zeros((2,3))
        fault_depth = np.zeros((2,3))
        displacements_numpy = displacements.loc[
            displacements['fname'] == f, ['vertical_displacement', 'downthrow_dir', 'dip_dir','X','Y']].to_numpy()
        # index = np.argmax(np.abs(displacements_numpy[:, 0]), )
        index = np.argsort(np.abs(displacements_numpy[:, 0]))[len(np.abs(displacements_numpy[:, 0]))//2]
        
        max_displacement[f] = displacements_numpy[
            index, 0]
        downthrow_dir[f] = displacements_numpy[index,[1,3,4]]
        if displacements_numpy[index,1] == 1.0:
            logger.info("Estimating downthrow direction using fault intersections")
            # fault_intersection_angles[f]
        if np.abs(displacements_numpy[index, 1] - displacements_numpy[index, 2]) > 90:
            # fault_orientations.loc[fault_orientations['formation'] == f, ['gx','gy','gy']]=-fault_orientations.loc[fault_orientations['formation'] == f, ['gx','gy','gy']]
            fault_orientations.loc[fault_orientations['formation'] == f, 'DipDirection'] -= 180#displacements_numpy[
                # index, 1]
        # find the middle of the fault as the mean of the line, average dip direction and the influence distance
        fault_centers[:3] = np.mean(fault_locations.loc[fault_locations['formation']==f,['X','Y','Z']],axis=0)
        fault_centers[3] = np.mean(fault_orientations.loc[fault_orientations['formation']==f,['DipDirection']])
        fault_centers[4] = fault_dimensions.loc[fault_dimensions['Fault']==f,'InfluenceDistance']
        fault_centers[5] = fault_dimensions.loc[fault_dimensions['Fault']==f,'HorizontalRadius']
        if type(fault_slip_vector) == dict:
            fault_slip = fault_slip_vector[f]
        elif type(fault_slip_vector) == np.array:
            fault_slip = fault_slip_vector
        else:
            fault_slip = np.array([0.,0.,-1.]) 
        stratigraphic_column['faults'][f] = {'FaultCenter':fault_centers[:3],
                                            'FaultDipDirection':fault_centers[3],
                                            'InfluenceDistance':fault_dimensions.loc[fault_dimensions['Fault']==f,'InfluenceDistance'].to_numpy(),
                                            'HorizontalRadius':fault_dimensions.loc[fault_dimensions['Fault']==f,'HorizontalRadius'].to_numpy(),
                                            'VerticalRadius':fault_dimensions.loc[fault_dimensions['Fault']==f,'VerticalRadius'].to_numpy(),
                                            'FaultSlip':fault_slip,
                                            'FaultNorm':normal_vector}

        if 'colour' in fault_dimensions.columns:
            stratigraphic_column['faults'][f]['colour'] = fault_dimensions.loc[fault_dimensions['Fault']==f,'colour'].to_numpy()
        normal_vector[0] = np.sin(np.deg2rad(fault_centers[3]))
        normal_vector[1] = np.cos(np.deg2rad(fault_centers[3]))
        strike_vector[0] = normal_vector[1]
        strike_vector[1] = -normal_vector[0]
        slip_vector[2]=1
        fault_edges[0,:] = fault_centers[:3]+normal_vector*fault_centers[4]
        fault_edges[1,:] = fault_centers[:3]-normal_vector*fault_centers[4]
        fault_tips[0,:] = fault_centers[:3]+strike_vector*fault_centers[5]
        fault_tips[1,:] = fault_centers[:3]-strike_vector*fault_centers[5]
        # fault_depth[0,:] = fault_centers[:3]+slip_vector*fault_centers[5]
        # fault_depth[1,:] = fault_centers[:3]-slip_vector*fault_centers[5]
        # fault_locations.loc[len(fault_locations),['X','Y','Z','formation','val','coord']] = [fault_edges[0,0],fault_edges[0,1],fault_edges[0,2],f,1,0]
        # fault_locations.loc[len(fault_locations),['X','Y','Z','formation','val','coord']] = [fault_edges[1,0],fault_edges[1,1],fault_edges[1,2],f,-1,0]
        # fault_locations.loc[len(fault_locations),['X','Y','Z','formation','val','coord']] = [fault_tips[0,0],fault_tips[0,1],fault_tips[0,2],f,1,2]
        # fault_locations.loc[len(fault_locations),['X','Y','Z','formation','val','coord']] = [fault_tips[1,0],fault_tips[1,1],fault_tips[1,2],f,-1,2]
        # # add strike vector to constraint fault extent
        # fault_orientations.loc[len(fault_orientations),['X','Y','Z','formation','DipDirection','coord']] = [fault_centers[0],fault_centers[1],fault_centers[2],f, fault_centers[3]-90,2]
        # fault_orientations.loc[len(fault_orientations),['X','Y','Z','formation','dip','coord']] = [fault_centers[0],fault_centers[1],fault_centers[2],f, 0,2]

        # print('downthro',displacements_numpy[index, 1])
        
    fault_orientations['strike'] = fault_orientations['DipDirection'] - 90
    fault_orientations[['gx', 'gy', 'gz']] = strike_dip_vector(fault_orientations['strike'], fault_orientations['dip'])
    for g in groups['group'].unique():
        groups.loc[groups['group']==g,'group'] = supergroups[g]
    # fault_orientations['strike'] = fault_orientations['DipDirection'] - 90
    # fault_orientations['gx'] = np.nan
    # fault_orientations['gy'] = np.nan
    # fault_orientations['gz'] = np.nan

    fault_orientations.drop(['strike', 'DipDirection', 'dip', 'DipPolarity'], inplace=True, axis=1)
    fault_orientations['feature_name'] = fault_orientations['formation']
    fault_locations['feature_name'] = fault_locations['formation']

    data = pd.concat([tangents, contact_orientations, contacts, fault_orientations, fault_locations])
    data.reset_index(inplace=True)
    
    return {'data': data,
            'groups': groups,
            'max_displacement': max_displacement,
            'fault_fault': fault_fault_relations,
            'stratigraphic_column': stratigraphic_column,
            'bounding_box':bb,
            'strat_va':strat_val,
            'downthrow_dir':downthrow_dir,
            'fault_graph':fault_graph,
            'fault_intersection_angles':fault_intersection_angles,
            'thickness':thickness}

def build_model(m2l_data, 
                evaluate = True, 
                skip_faults = False, 
                unconformities = False, 
                fault_params = None, 
                foliation_params = None, 
                rescale = True,
                skip_features = [],
                **kwargs):
    """[summary]

    [extended_summary]

    Parameters
    ----------
    m2l_data : dict
        [description]
    skip_faults : bool, optional
        [description], by default False
    fault_params : dict, optional
        [description], by default None
    foliation_params : dict, optional
        [description], by default None

    Returns
    -------
    [type]
        [description]
    """
    from LoopStructural import GeologicalModel


    boundary_points = np.zeros((2, 3))
    boundary_points[0, 0] = m2l_data['bounding_box']['minx']
    boundary_points[0, 1] = m2l_data['bounding_box']['miny']
    boundary_points[0, 2] = m2l_data['bounding_box']['lower']
    boundary_points[1, 0] = m2l_data['bounding_box']['maxx']
    boundary_points[1, 1] = m2l_data['bounding_box']['maxy']
    boundary_points[1, 2] = m2l_data['bounding_box']['upper']

    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :], rescale=rescale)
    # m2l_data['data']['val'] /= model.scale_factor
    model.set_model_data(m2l_data['data'])
    if not skip_faults:
        faults = []
        for f in m2l_data['max_displacement'].keys():
            if model.data[model.data['feature_name'] == f].shape[0] == 0:
                continue
            if f in skip_features:
                continue
            fault_id = f
            overprints = []
            try:
                overprint_id = m2l_data['fault_fault'][m2l_data['fault_fault'][fault_id] == 1]['fault_id'].to_numpy()
                for i in overprint_id:
                    overprints.append(i)
                logger.info('Adding fault overprints {}'.format(f))
            except:
                logger.info('No entry for %s in fault_fault_relations' % f)
        #     continue
            fault_center = m2l_data['stratigraphic_column']['faults'][f]['FaultCenter']
            fault_influence = m2l_data['stratigraphic_column']['faults'][f]['InfluenceDistance']
            fault_extent = m2l_data['stratigraphic_column']['faults'][f]['HorizontalRadius']
            fault_vertical_radius = m2l_data['stratigraphic_column']['faults'][f]['VerticalRadius']
            fault_slip_vector=m2l_data['stratigraphic_column']['faults'][f]['FaultSlip'] 
            faults.append(model.create_and_add_fault(f,
                                                    -m2l_data['max_displacement'][f],
                                                    faultfunction='BaseFault',
                                                    fault_slip_vector=fault_slip_vector,
                                                    fault_center=fault_center,
                                                    fault_extent=fault_extent,
                                                    fault_influence=fault_influence,
                                                    fault_vectical_radius=fault_vertical_radius,
                                                    # overprints=overprints,
                                                    **fault_params,
                                                    )
                        )
    # for f in m2l_data['fault_intersection_angles']:
    #     if f in m2l_data['max_displacement'].keys():
    #         f1_norm = m2l_data['stratigraphic_column']['faults'][f]['FaultNorm']
    #         for intersection in m2l_data['fault_intersection_angles'][f]:
    #             if intersection[0] in m2l_data['max_displacement'].keys():
    #                 f2_norm = m2l_data['stratigraphic_column']['faults'][intersection[0]]['FaultNorm']
    #                 if intersection[2] < 30 and np.dot(f1_norm,f2_norm)>0:
    #                     logger.info('Adding splay {} to {}'.format(intersection[0],f))
    #                     if model[f] is None:
    #                         logger.error('Fault {} does not exist, cannot be added as splay')
    #                     elif model[intersection[0]] is None:
    #                         logger.error('Fault {} does not exist')
    #                     else:
    #                         model[intersection[0]].builder.add_splay(model[f])

    #                 else:
    #                     logger.info('Adding abut {} to {}'.format(intersection[0],f))
    #                     model[intersection[0]].add_abutting_fault(model[f])
    faults = m2l_data.get('fault_graph',None)
    if faults:
        for f in faults.nodes:
            f1_norm = m2l_data['stratigraphic_column']['faults'][f]['FaultNorm']
            for e in faults.edges(f):
                data = faults.get_edge_data(*e)
                f2_norm =  m2l_data['stratigraphic_column']['faults'][e[1]]['FaultNorm']
                

                if float(data['angle']) < 30 and np.dot(f1_norm,f2_norm)>0:
                    if model[f] is None or model[e[1]] is None:
                        logger.error('Fault {} does not exist, cannot be added as splay')
                    elif model[e[1]] is None:
                        logger.error('Fault {} does not exist')
                    else:
                        region = model[e[1]].builder.add_splay(model[f])
                        model[e[1]].splay[model[f].name] = region
                else:
                    if model[f] is None or model[e[1]] is None:
                        continue

                    logger.info('Adding abut {} to {}'.format(e[1],f))
                    model[e[1]].add_abutting_fault(model[f])
    ## loop through all of the groups and add them to the model in youngest to oldest.
    group_features = []
    for i in np.sort(m2l_data['groups']['group number'].unique()):
        g = m2l_data['groups'].loc[m2l_data['groups']['group number'] == i, 'group'].unique()[0]
        group_features.append(model.create_and_add_foliation(g,
                                                            **foliation_params))
        # if the group was successfully added (not null) then lets add the base (0 to be unconformity)
        if group_features[-1] and unconformities:
            model.add_unconformity(group_features[-1], 0)
    model.set_stratigraphic_column(m2l_data['stratigraphic_column'])
    if evaluate:
        model.update(verbose=True)
    return model