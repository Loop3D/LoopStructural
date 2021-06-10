import numpy as np
from ._plane_fit import PlaneFitFeature
from LoopStructural.utils import getLogger
logger = getLogger(__name__)
def displacement_missfit(d,
                            fault_slip,
                            fault_center, 
                            fault_influence, 
                            fault_extent, 
                            fault_vertical_radius, 
                            fault_normal,
                            fname,
                            model,
                            fault_params,
                            strati_name,
                            view=None):
    ## create the fault in the model
    if view is not None:
        view.clear()
    model.create_and_add_fault(fname,
                                d,
                                faultfunction='BaseFault',
                                fault_slip_vector=fault_slip,
                                fault_center=fault_center,
                                fault_extent=fault_extent,
                                fault_influence=fault_influence,
                                fault_vectical_radius=fault_vertical_radius,
#                                 overprints=overprints,
                                **fault_params,
                                )
    # run interpolator
    model[fname].builder.update()
    # determine which points are inside the fault volume
    v = model[fname].faultframe.evaluate_value(model.data[['X','Y','Z']].to_numpy())
    np.all(np.logical_and(v > -1,v<1),axis=1)
    mask = model[fname].inside_volume(model.data[['X','Y','Z']].to_numpy())
    
    data2 = model.data.copy()

    # mask only data associated with the faulted strati
    mask = np.logical_and(mask,model.data['feature_name']==strati_name)
    if view is not None:
        view.add_isosurface(model[fname],value=0)
    value_data = data2.loc[mask,'val'].to_numpy()
    # gx = model[fname][0].evaluate_value(data2.loc[np.logical_and(~np.isnan(data2['val']),mask),['X','Y','Z']].to_numpy(),inplace=False))
    normals =data2.loc[np.logical_and(np.isnan(data2['val']),mask),
                            ['nx','ny','nz']].to_numpy()
    if view is not None:
        view.add_vector_field(model[fname],
                            locations=data2.loc[mask,['X','Y','Z']].to_numpy()
                            )
        view.add_value_data(model.rescale(data2.loc[mask,['X','Y','Z']].to_numpy(),inplace=False),
                            data2.loc[mask,['val']].to_numpy(),name='pts_before',
                            vmin=np.nanmin(value_data),
                            vmax=np.nanmax(value_data),cmap='rainbow')
        view.add_orientation_disks(model.rescale(data2.loc[np.logical_and(np.isnan(data2['val']),mask),
                            ['X','Y','Z']].to_numpy(),inplace=False),
                                   normals,
                                   name='normals')
    data2.loc[mask,['X','Y','Z']] = model[fname].apply_to_points(data2.loc[mask,['X','Y','Z']].to_numpy())
    tmp_pts = data2.loc[mask,['X','Y','Z','val','nx','ny','nz']]
#     tmp_pts = tmp_pts.loc[~np.isnan(tmp_pts['val']),:]
    tmp_pts.loc[:,['X','Y','Z']] = model.rescale(tmp_pts.loc[:,['X','Y','Z']],inplace=False)
    fault_plane_feature = PlaneFitFeature(tmp_pts,'Plane_Fit_{}'.format(fname),model=model)
    
    dv =data2.loc[np.logical_and(~np.isnan(data2['val']),mask),'val'].to_numpy()
    estimated_normal = fault_plane_feature.evaluate_gradient(data2.loc[np.logical_and(np.isnan(data2['val']),mask),
                            ['X','Y','Z']].to_numpy())
    if view is not None:
        view.add_value_data(model.rescale(data2.loc[np.logical_and(~np.isnan(data2['val']),mask),
                            ['X','Y','Z']].to_numpy(),inplace=False),
                            dv,
                            name='points_after',
                            vmin=np.nanmin(value_data),
                            vmax=np.nanmax(value_data),
                            cmap='rainbow')
        view.add_orientation_disks(model.rescale(data2.loc[np.logical_and(np.isnan(data2['val']),mask),
                            ['X','Y','Z']].to_numpy(),inplace=False),
                                   estimated_normal,
                                   name='estimated_normal')
        view.add_isosurface(fault_plane_feature,
                            slices=np.unique(dv),
                            paint_with=fault_plane_feature,
                            vmin=np.nanmin(value_data),
                            vmax=np.nanmax(value_data),
                            cmap='rainbow')
    fv = fault_plane_feature.evaluate_value(data2.loc[np.logical_and(~np.isnan(data2['val']),mask),['X','Y','Z']].to_numpy())
    return fv, dv, fault_plane_feature, estimated_normal, normals
