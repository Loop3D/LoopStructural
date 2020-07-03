import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)

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
    tangents = pd.read_csv(m2l_directory + '/tmp/raw_contacts.csv')
    groups = pd.read_csv(m2l_directory + '/tmp/all_sorts.csv', index_col=0)
    contact_orientations = pd.read_csv(m2l_directory + '/output/orientations.csv')
    # formation_thickness = pd.read_csv)
    contacts = pd.read_csv(m2l_directory + '/output/contacts4.csv')
    displacements = pd.read_csv(m2l_directory + '/output/fault_displacements3.csv')
    fault_orientations = pd.read_csv(m2l_directory + '/output/fault_orientations.csv')
    fault_locations = pd.read_csv(m2l_directory + '/output/faults.csv')
    fault_fault_relations = pd.read_csv(m2l_directory + '/output/fault-fault-relationships.csv')
    fault_strat_relations = pd.read_csv(m2l_directory + '/output/group-fault-relationships.csv')
    supergroups = {}
    sgi = 0
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
    

    bb = pd.read_csv(m2l_directory+'/tmp/bbox.csv')

    # process tangent data to be tx, ty, tz
    tangents['tz'] = 0
    tangents['tx'] = tangents['lsx']
    tangents['ty'] = tangents['lsy']
    tangents.drop(['angle', 'lsx', 'lsy'], inplace=True, axis=1)

    # convert azimuth and dip to gx, gy, gz
    from LoopStructural.utils.helper import strike_dip_vector
    contact_orientations['strike'] = contact_orientations['azimuth'] - 90
    contact_orientations['nx'] = np.nan
    contact_orientations['ny'] = np.nan
    contact_orientations['nz'] = np.nan
    contact_orientations[['nx', 'ny', 'nz']] = strike_dip_vector(contact_orientations['strike'],
                                                                 contact_orientations['dip'])
    contact_orientations.drop(['strike', 'dip', 'azimuth'], inplace=True, axis=1)

    # calculate scalar field values
    i = 0
    thickness = {}
    with open(m2l_directory + '/output/formation_summary_thicknesses.csv') as file:
        for l in file:
            if i>1:
                linesplit = l.split(',')
                thickness[linesplit[0]] = float(linesplit[1])
    #             print(l.split(',')[1])
            i+=1
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

        for c in groups.loc[groups['group number'] == i, 'code']:
            strat_val[c] = np.nan
            if c in thickness:
                stratigraphic_column[g][c] = {'min': val[g], 'max': val[g] + thickness[c], 'id': unit_id}
                unit_id += 1
                strat_val[c] = val[g]
                val[g] += thickness[c]
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
    contacts['val'] = np.nan
    for o in strat_val:
        contacts.loc[contacts['formation'] == o, 'val'] = strat_val[o]

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
    for f in displacements['fname'].unique():
        displacements_numpy = displacements.loc[
            displacements['fname'] == f, ['vertical_displacement', 'downthrow_dir', 'dip_dir']].to_numpy()
        # index = np.argmax(np.abs(displacements_numpy[:, 0]), )
        index = np.argsort(np.abs(displacements_numpy[:, 0]))[len(np.abs(displacements_numpy[:, 0]))//2]

        max_displacement[f] = displacements_numpy[
            index, 0]
        if displacements_numpy[index, 1] - displacements_numpy[index, 2] > 90:
            fault_orientations.loc[fault_orientations['formation'] == fname, 'DipDirection'] = displacements_numpy[
                index, 1]
        # .loc[displacements['fname'] == f,'vertical_displacement'].max()
    for g in groups['group'].unique():
        groups.loc[groups['group']==g,'group'] = supergroups[g]
    fault_orientations['strike'] = fault_orientations['DipDirection'] - 90
    fault_orientations['gx'] = np.nan
    fault_orientations['gy'] = np.nan
    fault_orientations['gz'] = np.nan

    fault_orientations[['gx', 'gy', 'gz']] = strike_dip_vector(fault_orientations['strike'], fault_orientations['dip'])
    fault_orientations.drop(['strike', 'DipDirection', 'dip', 'DipPolarity'], inplace=True, axis=1)
    fault_orientations['feature_name'] = fault_orientations['formation']

    fault_locations['val'] = 0
    fault_locations['feature_name'] = fault_locations['formation']


    data = pd.concat([tangents, contact_orientations, contacts, fault_orientations, fault_locations])
    data.reset_index()

    return {'data': data,
            'groups': groups,
            'max_displacement': max_displacement,
            'fault_fault': fault_fault_relations,
            'stratigraphic_column': stratigraphic_column,
            'bounding_box':bb,
            'strat_va':strat_val}

def build_model(m2l_data, skip_faults = False, unconformities=False, fault_params = None, foliation_params=None, rescale = True):
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
    model.set_model_data(m2l_data['data'])
    if not skip_faults:
        faults = []
        for f in m2l_data['max_displacement'].keys():
            if model.data[model.data['feature_name'] == f].shape[0] == 0:
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
            faults.append(model.create_and_add_fault(f,
                                                    -m2l_data['max_displacement'][f],
                                                    faultfunction='BaseFault',
                                                    overprints=overprints,
                                                    **fault_params,
                                                    )
                        )

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
    return model