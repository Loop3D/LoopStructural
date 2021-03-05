import numpy as np
from ._plane_fit import PlaneFitFeature
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
                            view=None):
    ## create the fault in the model
    if view:
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
    mask = model[fname].inside_volume(model.data[['X','Y','Z']].to_numpy())
    # mask only data associated with the faulted strati
    mask = np.logical_and(mask,model.data['feature_name']=='supergroup_0')
    if view:
        view.add_isosurface(model[fname],value=0)
    data2 = model.data.copy()
    value_data = data2.loc[mask,'val'].to_numpy()
    # gx = model[fname][0].evaluate_value(data2.loc[np.logical_and(~np.isnan(data2['val']),mask),['X','Y','Z']].to_numpy(),inplace=False))
    if view:
        view.add_vector_field(model[fname],
                            locations=model.rescale(data2.loc[mask,['X','Y','Z']].to_numpy(),
                            inplace=False))
        view.add_value_data(model.rescale(data2.loc[mask,['X','Y','Z']].to_numpy(),inplace=False),
                            data2.loc[mask,['val']].to_numpy(),name='pts_before',
                            vmin=np.nanmin(value_data),
                            vmax=np.nanmax(value_data),cmap='rainbow')
    data2.loc[mask,['X','Y','Z']] = model[fname].apply_to_points(data2.loc[mask,['X','Y','Z']].to_numpy())
    tmp_pts = data2.loc[mask,['X','Y','Z','val','nx','ny','nz']]
    tmp_pts = tmp_pts.loc[~np.isnan(tmp_pts['val']),:]
    fault_plane_feature = PlaneFitFeature(tmp_pts,'Plane_Fit_{}'.format(fname),model=model)
    
    dv =data2.loc[np.logical_and(~np.isnan(data2['val']),mask),'val'].to_numpy()
    if view:
        view.add_value_data(data2.loc[np.logical_and(~np.isnan(data2['val']),mask),
                            ['X','Y','Z']].to_numpy(),
                            dv,
                            vmin=np.nanmin(value_data),
                            vmax=np.nanmax(value_data),
                            cmap='rainbow')
        view.add_isosurface(fault_plane_feature,
                            values=np.unique(dv),
                            paint_with=fault_plane_feature,
                            vmin=np.nanmin(value_data),
                            vmax=np.nanmax(value_data),
                            cmap='rainbow')
    fv = fault_plane_feature.evaluate_value(data2.loc[np.logical_and(~np.isnan(data2['val']),mask),['X','Y','Z']].to_numpy())
    return fv, dv, fault_plane_feature
