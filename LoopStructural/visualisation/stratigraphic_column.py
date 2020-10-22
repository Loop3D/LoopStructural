import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

class StratigraphicColumnView:
    def __init__(self,model,ax=None,cmap=None, labels = None):
        n_units = 0 #count how many discrete colours
        xmin = 0
        ymin = 0
        ymax = 1
        xmax = 1
        if ax is None:
            fig, ax = plt.subplots(figsize=(2,10))
        patches = []
        for g in model.stratigraphic_column.keys():
            if g == 'faults':
                continue
            for u in model.stratigraphic_column[g].keys():
                n_units+=1
                ymin = -model.stratigraphic_column[g][u]['min']
                if np.isinf(model.stratigraphic_column[g][u]['min']):
                    ymin = 0
                ymax = -model.stratigraphic_column[g][u]['max']
                if np.isinf(ymax):
                    ymin = ymax + (ymax-ymin)*(1+np.random.rand())
                polygon_points = np.array([[xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]])
                patches.append(Polygon(polygon_points))
                xy = (0,ymin+(ymax-ymin)/2)
                if labels:
                    ax.annotate(labels[u],xy)
                else:
                    ax.annotate(u,xy)
        if cmap is None:
            import matplotlib.colors as colors
            colours = []
            boundaries = []
            data = []
            for g in model.stratigraphic_column.keys():
                if g == 'faults':
                    continue
                for u, v  in model.stratigraphic_column[g].items():
                    data.append((v['id'],v['colour']))
                    colours.append(v['colour'])
                    boundaries.append(v['id'])#print(u,v)
            cmap = colors.ListedColormap(colours)
        else:
            cmap = cm.get_cmap(cmap,n_units-1)
        ci = 0
        p = PatchCollection(patches, cmap=cmap)

        colors = np.arange(len(patches))
        p.set_array(np.array(colors))

        ax.add_collection(p)

        ax.set_ylim(ymax+(ymax-ymin)*-2,0)#ax.set_ylim(0,ymax)
        ax.axis('off')
        