import logging

import matplotlib.pyplot as plt
import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class MapView:
    """

    """
    def __init__(self, model = None, bounding_box=np.zeros((2,2)), nsteps=None, ax = None, **kwargs):
        """

        Parameters
        ----------
        origin - lower left
        maximum - upper right
        nsteps - number of cells
        kwargs
        """
        
        self.xx = None
        self.yy = None
        
        self._bounding_box = bounding_box
        self._nsteps = nsteps
        if self._nsteps is not None and self._bounding_box is not None:
            self._update_grid()
        if model is not None:
            #make sure self._nsteps is 2d
            self.model = model
        self.ax = ax
        if self.ax is None:
            fig, self.ax = plt.subplots(1, figsize=(10, 10))
        self.ax.set_aspect('equal', adjustable='box')

        # set plot limits to model bounding box
        self._xmin = self.bounding_box[0,0]
        self._xmax = self.bounding_box[1,0]
        self._ymin = self.bounding_box[0,1]
        self._ymax = self.bounding_box[1,1]
        self._update_plot_limits()
    @property
    def model(self):
        return self._model

    @model.setter
    def model(self,model):
        if model is not None:
            bb = np.array([model.origin[:2],model.maximum[:2]])
            self.bounding_box = bb#model.bounding_box
            self.nsteps = model.nsteps[:2]
            self._model = model
            self._update_grid()

    @property
    def nsteps(self):
        return self._nsteps
    
    @nsteps.setter
    def nsteps(self,nsteps):
        if len(nsteps) != 2:
            logger.error("Can't update nsteps, needs to be 2D")
            return
        self._nsteps = nsteps
        self._update_grid()
        
    @property
    def bounding_box(self):
        return self._bounding_box
    @bounding_box.setter
    def bounding_box(self,bounding_box):
        self._bounding_box = bounding_box
        self._update_grid()

    @property
    def xmin(self):
        return self._xmin
    
    @xmin.setter
    def xmin(self,xmin):
        self._xmin = xmin
        self._update_plot_limits()
    
    @property
    def xmax(self):
        return self._xmax
    
    @xmax.setter
    def xmax(self,xmax):
        self._xmax = xmax
        self._update_plot_limits()

    @property
    def ymin(self):
        return self._ymin
    
    @ymin.setter
    def ymin(self,ymin):
        self._ymin = ymin
        self._update_plot_limits()

    @property
    def ymax(self):
        return self._ymax
    
    @ymax.setter
    def ymax(self,ymax):
        self._ymax = ymax
        self._update_plot_limits()

    def _update_plot_limits(self):
        self.ax.set_xlim([self._xmin,self._xmax])
        self.ax.set_ylim([self._ymin,self._ymax])

    
    def _update_grid(self):
        """Internal function to update the current grid when the bounding box
        or number of steps changes
        """
        if self.nsteps is None or self.bounding_box is None:
            return        
        x = np.linspace(self.bounding_box[0,0], self.bounding_box[1,0], self.nsteps[0])
        y = np.linspace(self.bounding_box[0,1], self.bounding_box[1,1], self.nsteps[1])
        self.xx, self.yy = np.meshgrid(x, y, indexing='ij')
        self.xx = self.xx.flatten()
        self.yy = self.yy.flatten()

    def add_data(self, feature, val=True, grad=True, unfault=False, dip=True, **kwargs):
        """
        Adds the data associated to the feature to the plot
        Parameters
        ----------
        feature : GeologicalFeature
            the feature whose data you want to add
        val : bool
            whether to add value data
        grad : bool
            whether to add gradient data
        unfault : bool
            plot points in their restored location
        dip : bool
            whether to annotate the dip, default False

        kwargs are passed to matplotlib functions and draw strike
        Returns
        -------

        """
        # logger.warning("Plotting restored data locations")
        ori_data = []
        gradient_data = feature.builder.get_gradient_constraints()
        if unfault:
            gradient_data = feature.interpolator.get_gradient_constraints()

        if gradient_data.shape[0] > 0:
            ori_data.append(gradient_data)
        norm_data = feature.builder.get_norm_constraints()
        if unfault:
            norm_data = feature.interpolator.get_norm_constraints()

        if norm_data.shape[0] > 0:
            ori_data.append(norm_data)
        cmap = kwargs.pop('cmap','rainbow')
        # if single colour then specify kwarg, otherwise use point value
        if val:
            value_data = np.copy(feature.builder.get_value_constraints())
            if unfault:
                value_data = np.copy(feature.interpolator.get_value_constraints())
              
            value_data[:,:3] = self.model.rescale(value_data[:,:3],inplace=False)
            point_colour = kwargs.pop('point_colour',None)
            if point_colour is None:
                self.ax.scatter(value_data[:, 0], value_data[:, 1], c=value_data[:,3],
                            vmin=feature.min(), vmax=feature.max(),cmap=cmap)
            if point_colour is not None:
                self.ax.scatter(value_data[:, 0], value_data[:, 1], c=point_colour)
        if grad:
            symb_colour = kwargs.pop('symb_colour','black')
            symb_scale=kwargs.pop('symb_scale',1.)
            gradient_data = np.hstack(ori_data)
            gradient_data[:,:3] = self.model.rescale(gradient_data[:,:3],inplace=False)
            gradient_data[:, 3:5] /= np.linalg.norm(gradient_data[:, 3:5], axis=1)[:, None]
            t = gradient_data[:, [4, 3]] * np.array([1, -1]).T
            n = gradient_data[:, 3:5]
            t *= symb_scale
            n *= 0.5*symb_scale
            p1 = gradient_data[:, [0, 1]] - t
            p2 = gradient_data[:, [0, 1]] + t
            # plt.scatter(val[:,0],val[:,1],c='black')
            self.ax.plot([p1[:, 0], p2[:, 0]], [p1[:, 1], p2[:, 1]], symb_colour)
            p1 = gradient_data[:, [0, 1]]
            p2 = gradient_data[:, [0, 1]] + n
            self.ax.plot([p1[:, 0], p2[:, 0]], [p1[:, 1], p2[:, 1]], symb_colour)
            if dip:
                dip_v = np.rad2deg(np.arccos(gradient_data[:,5])).astype(int)
                for d, xy, v in zip(dip_v,gradient_data[:,:2],gradient_data[:,3:6]):
                    self.ax.annotate(d,xy,xytext=xy+v[:2]*symb_scale*.1,fontsize='small')
                    
    def add_fault_ellipse(self, faults=None, **kwargs):
        from matplotlib.patches import Ellipse
        for k, f in self.model.stratigraphic_column['faults'].items():
            center = self.model.rescale(f['FaultCenter'])
            e = Ellipse((center[0],center[1]),
                            f['HorizontalRadius']*2,
                            f['InfluenceDistance']*2,
                            360-f['FaultDipDirection'],
                                facecolor='None',edgecolor='k')

        self.ax.add_patch(e)

    def add_scalar_field(self, feature, z=0, **kwargs):
        """
        Plot the scalar field value on a map

        Parameters
        ----------
        feature : GeologicalFeature
            which feature to plot on the map
        z : double/np.array
            height
        kwargs

        Returns
        -------

        """
        zz = np.zeros(self.xx.shape)
        zz[:] = z
        v = feature.evaluate_value(self.model.scale(np.array([self.xx, self.yy, zz]).T,inplace=False))
        return self.ax.imshow(v.reshape(self.nsteps).T,
                       extent=[self.bounding_box[0,0], self.bounding_box[1,0], self.bounding_box[0,1],
                               self.bounding_box[1,1]],
                       vmin=feature.min(), vmax=feature.max(),
                       origin='lower',
                       **kwargs)

    def add_contour(self, feature, values, z=0, mask = None, **kwargs):
        """Add an isoline of a scalar field to the map

        Parameters
        ----------
        feature : GeologicalFeature
            the feature to isosurface
        values : list 
            list of values to contour
        z : double/np.array, optional
            elevation of map, by default 0
        """        
        zz = np.zeros(self.xx.shape)
        zz[:] = z
        v = feature.evaluate_value(self.model.scale(np.array([self.xx, self.yy, zz]).T,inplace=False))
        if mask:
            maskv = mask(self.model.scale(np.array([self.xx, self.yy, zz]).T,inplace=False))
            v[~maskv] = np.nan
        return self.ax.contour(v.reshape(self.nsteps).T,extent=[self.bounding_box[0,0], self.bounding_box[1,0], self.bounding_box[0,1],
                                    self.bounding_box[1,1]],origin='lower',levels=values,**kwargs
                       )

        
    
    def add_model(self, z = 0,cmap=None):
        """Plot the model onto a map

        Parameters
        ----------
        z : int/numpy array, optional
            height of the map surface (could also be a dem), by default 0
        cmap : str/matplotlib colourmap, optional
            specify a colour map, by default 'tab20'
        """  
        if cmap is None:
            import matplotlib.colors as colors
            colours = []
            boundaries = []
            data = []
            for g in self.model.stratigraphic_column.keys():
                if g == 'faults':
                    continue
                for u, v  in self.model.stratigraphic_column[g].items():
                    data.append((v['id'],v['colour']))
                    colours.append(v['colour'])
                    boundaries.append(v['id'])#print(u,v)
            cmap = colors.ListedColormap(colours)

        zz = np.zeros_like(self.xx)
        zz[:] = z#self.bounding_box[1,2]
        pts = np.vstack([self.xx.flatten(),self.yy.flatten(),zz.flatten()])
        if self.model is None:
            logger.error("Mapview needs a model assigned to plot model on map")
            return 
        vals = self.model.evaluate_model(pts.T,scale=True)
        return self.ax.imshow(vals.reshape(self.nsteps).T,extent=[self.bounding_box[0,0], self.bounding_box[1,0], self.bounding_box[0,1],
                                    self.bounding_box[1,1]],origin='lower',cmap=cmap)
                                
    def add_fault_displacements(self,z = 0, cmap='rainbow'):
        
        zz = np.zeros_like(self.xx)
        zz[:] = z#self.bounding_box[1,2]
        pts = np.vstack([self.xx.flatten(),self.yy.flatten(),zz.flatten()])
        if self.model is None:
            logger.error("Mapview needs a model assigned to plot model on map")
            return 
        vals = self.model.evaluate_fault_displacements(pts.T,scale=True)
        return self.ax.imshow(vals.reshape(self.nsteps).T,extent=[self.bounding_box[0,0], self.bounding_box[1,0], self.bounding_box[0,1],
                                    self.bounding_box[1,1]],origin='lower',cmap=cmap)

    def add_faults(self,**kwargs):
        for f in self.model.features:
            if f.type=='fault':
                # create a function to return true if displacement > 0
                def mask(x):
                    val = f.displacementfeature.evaluate_value(x)
                    val[np.isnan(val)] = 0
                    maskv = np.zeros(val.shape).astype(bool)
                    maskv[np.abs(val) > 0.001] = 1
                    return maskv
                
                self.add_contour(f,0,mask=mask,**kwargs)