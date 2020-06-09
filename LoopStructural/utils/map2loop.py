import pandas as pd
import numpy as np

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
    formation_thickness = pd.read_csv(m2l_directory + '/output/formation_thicknesses.csv')
    contacts = pd.read_csv(m2l_directory + '/output/contacts4.csv')
    displacements = pd.read_csv(m2l_directory + '/output/fault_displacements3.csv')
    fault_orientations = pd.read_csv(m2l_directory + '/output/fault_orientations.csv')
    fault_locations = pd.read_csv(m2l_directory + '/output/faults.csv')
    fault_fault_relations = pd.read_csv(m2l_directory + '/output/fault-fault-relationships.csv')
    fault_strat_relations = pd.read_csv(m2l_directory + '/output/group-fault-relationships.csv')
    supergroups = {}
    sgi = 0
    with open(m2l_directory + '/tmp/super_groups.csv') as f:
        for l in f:
             for g in l.split(','):
                g = g.replace('-','_').replace(' ','_')
                if g.find('\n') > 0:
                    g = g[:g.find('\n')]
                supergroups[g] = 'supergroup_{}'.format(sgi)
                if g == '\n':
                    sgi += 1
                    break

    bb = pd.read_csv(m2l_directory+'/tmp/bbox.csv')

    # process tangent data to be tx, ty, tz
    tangents['tz'] = 0
    tangents['tx'] = tangents['lsx']
    tangents['ty'] = tangents['lsy']
    tangents.drop(['angle', 'lsx', 'lsy'], inplace=True, axis=1)

    # convert azimuth and dip to gx, gy, gz
    from LoopStructural.utils.helper import strike_dip_vector
    contact_orientations['strike'] = contact_orientations['azimuth'] - 90
    contact_orientations['gx'] = np.nan
    contact_orientations['gy'] = np.nan
    contact_orientations['gz'] = np.nan
    contact_orientations[['gx', 'gy', 'gz']] = strike_dip_vector(contact_orientations['strike'],
                                                                 contact_orientations['dip'])
    contact_orientations.drop(['strike', 'dip', 'azimuth'], inplace=True, axis=1)

    # calculate scalar field values
    thickness = {}
    for f in formation_thickness['formation'].unique():
        thickness[f] = np.mean(formation_thickness[formation_thickness['formation'] == f]['thickness'])

    strat_val = {}
    stratigraphic_column = {}
    unit_id = 0

    for i in groups['group number'].unique():
        g = supergroups[groups.loc[groups['group number'] == i, 'group'].iloc[0]]
        if g not in stratigraphic_column:
            stratigraphic_column[g] = {}
            val = 0
        #         print(groups.loc[groups['group number']==i,'code'])
        for c in groups.loc[groups['group number'] == i, 'code'][::-1]:
            # roups.loc[groups['group number']==i
            strat_val[c] = np.nan
            if c in thickness:
                stratigraphic_column[g][c] = {'min': val, 'max': val + thickness[c], 'id': unit_id}
                unit_id += 1
                strat_val[c] = val
                val += thickness[c]
    contacts['val'] = np.nan
    for o in strat_val:
        contacts.loc[contacts['formation'] == o, 'val'] = strat_val[o]

    tangents['type'] = tangents['group']
    contact_orientations['type'] = None
    contacts['type'] = None
    for g in groups['group'].unique():
        val = 0
        for c in groups.loc[groups['group'] == g, 'code']:
            contact_orientations.loc[contact_orientations['formation'] == c, 'type'] = supergroups[g]
            contacts.loc[contacts['formation'] == c, 'type'] = supergroups[g]
    displacements['dip_dir'] = np.nan
    for fname in fault_orientations['formation'].unique():
        displacements.loc[displacements['fname'] == fname, 'dip_dir'] = np.mean(
            fault_orientations.loc[fault_orientations['formation'] == fname, 'DipDirection'])
    max_displacement = {}
    for f in displacements['fname'].unique():
        displacements_numpy = displacements.loc[
            displacements['fname'] == f, ['vertical_displacement', 'downthrow_dir', 'dip_dir']].to_numpy()
        index = np.argmax(np.abs(displacements_numpy[:, 0]), )
        max_displacement[f] = displacements_numpy[
            index, 0]
        if displacements_numpy[index, 1] - displacements_numpy[index, 2] > 90:
            fault_orientations.loc[fault_orientations['formation'] == fname, 'DipDirection'] = displacements_numpy[
                index, 1]
        # .loc[displacements['fname'] == f,'vertical_displacement'].max()
    fault_orientations['strike'] = fault_orientations['DipDirection'] - 90
    fault_orientations['gx'] = np.nan
    fault_orientations['gy'] = np.nan
    fault_orientations['gz'] = np.nan

    fault_orientations[['gx', 'gy', 'gz']] = strike_dip_vector(fault_orientations['strike'], fault_orientations['dip'])
    fault_orientations.drop(['strike', 'DipDirection', 'dip', 'DipPolarity'], inplace=True, axis=1)
    fault_orientations['type'] = fault_orientations['formation']

    fault_locations['val'] = 0
    fault_locations['type'] = fault_locations['formation']


    data = pd.concat([tangents, contact_orientations, contacts, fault_orientations, fault_locations])
    data.reset_index()

    return {'data': data,
            'groups': groups,
            'max_displacement': max_displacement,
            'fault_fault': fault_fault_relations,
            'stratigraphic_column': stratigraphic_column,
            'bounding_box':bb}