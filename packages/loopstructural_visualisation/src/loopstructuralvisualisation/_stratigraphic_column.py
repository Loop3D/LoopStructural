import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from LoopStructural.utils import rng


class StratigraphicColumnView:
    def __init__(self, model, ax=None, cmap=None, labels=None):
        self.model = model
        self.ax = ax
        self.cmap = cmap
        self.labels = labels

    def plot(self):
        n_units = 0  # count how many discrete colours (number of stratigraphic units)
        xmin = 0
        ymin = 0
        ymax = 1
        xmax = 1
        fig = None
        if self.ax is None:
            fig, self.ax = plt.subplots(figsize=(2, 10))
        patches = []  # stores the individual stratigraphic unit polygons

        total_height = 0
        prev_coords = [0, 0]

        # iterate through groups, skipping faults

        for g in reversed(self.model.stratigraphic_column.get_groups()):
            for u in g.units:
                n_units += 1

                ymax = total_height
                ymin = ymax - (u.thickness)

                if not np.isfinite(ymin):
                    ymin = prev_coords[1] - (prev_coords[1] - prev_coords[0]) * (1 + rng.random())

                total_height = ymin

                prev_coords = (ymin, ymax)

                polygon_points = np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])
                patches.append(Polygon(polygon_points))
                xy = (0, ymin + (ymax - ymin) / 2)
                if self.labels:
                    self.ax.annotate(self.labels[u], xy)
                else:
                    self.ax.annotate(u.name, xy)

        if self.cmap is None:
            import matplotlib.colors as colors

            colours = []
            boundaries = []
            data = []
            for g in self.model.stratigraphic_column.get_groups():
                if g == "faults":
                    continue
                for u in g.units:
                    data.append((u.id, u.colour))
                    colours.append(u.colour)
                    boundaries.append(u.id)  # print(u,v)
            cmap = colors.ListedColormap(colours)
        else:
            cmap = cm.get_cmap(self.cmap, n_units - 1)
        p = PatchCollection(patches, cmap=cmap)

        colors = np.arange(len(patches))
        p.set_array(np.array(colors))

        self.ax.add_collection(p)

        self.ax.set_ylim(total_height - (total_height - prev_coords[0]) * 0.1, 0)

        self.ax.axis("off")

        return fig
