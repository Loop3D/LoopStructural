import enum
from typing import Dict, Optional, List, Tuple
import numpy as np
from LoopStructural.utils import rng, getLogger, Observable
logger = getLogger(__name__)
logger.info("Imported LoopStructural Stratigraphic Column module")
class UnconformityType(enum.Enum):
    """
    An enumeration for different types of unconformities in a stratigraphic column.
    """

    ERODE = 'erode'
    ONLAP = 'onlap'


class StratigraphicColumnElementType(enum.Enum):
    """
    An enumeration for different types of elements in a stratigraphic column.
    """

    UNIT = 'unit'
    UNCONFORMITY = 'unconformity'


class StratigraphicColumnElement:
    """
    A class to represent an element in a stratigraphic column, which can be a unit or a topological object
    for example unconformity.
    """

    def __init__(self, uuid=None):
        """
        Initializes the StratigraphicColumnElement with a uuid.
        """
        if uuid is None:
            import uuid as uuid_module

            uuid = str(uuid_module.uuid4())
        self.uuid = uuid


class StratigraphicUnit(StratigraphicColumnElement):
    """
    A class to represent a stratigraphic unit.
    """

    def __init__(self, *, uuid=None, name=None, colour=None, thickness=None, data=None):
        """
        Initializes the StratigraphicUnit with a name and an optional description.
        """
        super().__init__(uuid)
        self.name = name
        if colour is None:
            colour = rng.random(3)
        self.colour = colour
        self.thickness = thickness
        self.data = data
        self.element_type = StratigraphicColumnElementType.UNIT

    def to_dict(self):
        """
        Converts the stratigraphic unit to a dictionary representation.
        """
        colour = self.colour
        if isinstance(colour, np.ndarray):
            colour = colour.astype(float).tolist()
        return {"name": self.name, "colour": colour, "thickness": self.thickness, 'uuid': self.uuid}

    @classmethod
    def from_dict(cls, data):
        """
        Creates a StratigraphicUnit from a dictionary representation.
        """
        if not isinstance(data, dict):
            raise TypeError("Data must be a dictionary")
        name = data.get("name")
        colour = data.get("colour")
        thickness = data.get("thickness", None)
        uuid = data.get("uuid", None)
        return cls(uuid=uuid, name=name, colour=colour, thickness=thickness)

    def __str__(self):
        """
        Returns a string representation of the stratigraphic unit.
        """
        return (
            f"StratigraphicUnit(name={self.name}, colour={self.colour}, thickness={self.thickness})"
        )


class StratigraphicUnconformity(StratigraphicColumnElement):
    """
    A class to represent a stratigraphic unconformity, which is a surface of discontinuity in the stratigraphic record.
    """

    def __init__(
        self, *, uuid=None, name=None, unconformity_type: UnconformityType = UnconformityType.ERODE
    ):
        """
        Initializes the StratigraphicUnconformity with a name and an optional description.
        """
        super().__init__(uuid)

        self.name = name
        if unconformity_type not in [UnconformityType.ERODE, UnconformityType.ONLAP]:
            raise ValueError("Invalid unconformity type")
        self.unconformity_type = unconformity_type
        self.element_type = StratigraphicColumnElementType.UNCONFORMITY

    def to_dict(self):
        """
        Converts the stratigraphic unconformity to a dictionary representation.
        """
        return {
            "uuid": self.uuid,
            "name": self.name,
            "unconformity_type": self.unconformity_type.value,
        }

    def __str__(self):
        """
        Returns a string representation of the stratigraphic unconformity.
        """
        return (
            f"StratigraphicUnconformity(name={self.name}, "
            f"unconformity_type={self.unconformity_type.value})"
        )

    @classmethod
    def from_dict(cls, data):
        """
        Creates a StratigraphicUnconformity from a dictionary representation.
        """
        if not isinstance(data, dict):
            raise TypeError("Data must be a dictionary")
        name = data.get("name")
        unconformity_type = UnconformityType(
            data.get("unconformity_type", UnconformityType.ERODE.value)
        )
        uuid = data.get("uuid", None)
        return cls(uuid=uuid, name=name, unconformity_type=unconformity_type)
class StratigraphicGroup:
    """
    A class to represent a group of stratigraphic units.
    This class is not fully implemented and serves as a placeholder for future development.
    """

    def __init__(self, name=None, units=None):
        """
        Initializes the StratigraphicGroup with a name and an optional list of units.
        """
        self.name = name
        self.units = units if units is not None else []


class StratigraphicColumn(Observable['StratigraphicColumn']):
    """
    A class to represent a stratigraphic column, which is a vertical section of the Earth's crust
    showing the sequence of rock layers and their relationships.
    """

    def __init__(self):
        """
        Initializes the StratigraphicColumn with a name and a list of layers.
        """
        super().__init__()
        self.order = [StratigraphicUnit(name='Basement', colour='grey', thickness=np.inf),StratigraphicUnconformity(name='Base Unconformity', unconformity_type=UnconformityType.ERODE)]
        self.group_mapping = {}
    def clear(self,basement=True):
        """
        Clears the stratigraphic column, removing all elements.
        """
        if basement:
            self.order = [StratigraphicUnit(name='Basement', colour='grey', thickness=np.inf),StratigraphicUnconformity(name='Base Unconformity', unconformity_type=UnconformityType.ERODE)]
        else:
            self.order = []
        self.group_mapping = {}
        self.notify('column_cleared')
    def add_unit(self, name,*, colour=None, thickness=None, where='top'):
        unit = StratigraphicUnit(name=name, colour=colour, thickness=thickness)

        if where == 'top':
            self.order.append(unit)
        elif where == 'bottom':
            self.order.insert(0, unit)
        else:
            raise ValueError("Invalid 'where' argument. Use 'top' or 'bottom'.")
        self.notify('unit_added', unit=unit)
        return unit

    def remove_unit(self, uuid):
        """
        Removes a unit or unconformity from the stratigraphic column by its uuid.
        """
        for i, element in enumerate(self.order):
            if element.uuid == uuid:
                del self.order[i]
                self.notify('unit_removed', uuid=uuid)
                return True
        
        return False

    def add_unconformity(self, name, *, unconformity_type=UnconformityType.ERODE, where='top' ):
        unconformity = StratigraphicUnconformity(
            uuid=None, name=name, unconformity_type=unconformity_type
        )

        if where == 'top':
            self.order.append(unconformity)
        elif where == 'bottom':
            self.order.insert(0, unconformity)
        else:
            raise ValueError("Invalid 'where' argument. Use 'top' or 'bottom'.")
        self.notify('unconformity_added', unconformity=unconformity)
        return unconformity

    def get_element_by_index(self, index):
        """
        Retrieves an element by its index from the stratigraphic column.
        """
        if index < 0 or index >= len(self.order):
            raise IndexError("Index out of range")
        return self.order[index]

    def get_unit_by_name(self, name):
        """
        Retrieves a unit by its name from the stratigraphic column.
        """
        for unit in self.order:
            if isinstance(unit, StratigraphicUnit) and unit.name == name:
                return unit

        return None

    def get_unconformity_by_name(self, name):
        """
        Retrieves an unconformity by its name from the stratigraphic column.
        """
        for unconformity in self.order:
            if isinstance(unconformity, StratigraphicUnconformity) and unconformity.name == name:
                return unconformity

        return None
    def get_element_by_uuid(self, uuid):
        """
        Retrieves an element by its uuid from the stratigraphic column.
        """
        for element in self.order:
            if element.uuid == uuid:
                return element
        raise KeyError(f"No element found with uuid: {uuid}")
    
    def get_group_for_unit_name(self, unit_name:str) -> Optional[StratigraphicGroup]:
        """
        Retrieves the group for a given unit name.
        """
        for group in self.get_groups():
            if any(unit.name == unit_name for unit in group.units):
                return group
        return None
    def add_element(self, element):
        """
        Adds a StratigraphicColumnElement to the stratigraphic column.
        """
        if isinstance(element, StratigraphicColumnElement):
            self.order.append(element)
        else:
            raise TypeError("Element must be an instance of StratigraphicColumnElement")

    def get_elements(self):
        """
        Returns a list of all elements in the stratigraphic column.
        """
        return self.order

    def get_groups(self):
        groups = []
        i=0
        group = StratigraphicGroup(
            name=(
                f'Group_{i}'
                if f'Group_{i}' not in self.group_mapping
                else self.group_mapping[f'Group_{i}']
            )
        )
        for e in reversed(self.order):
            if isinstance(e, StratigraphicUnit):
                group.units.append(e)
            else:
                if group.units:
                    groups.append(group)
                    i+=1
                    group = StratigraphicGroup(
                        name=(
                            f'Group_{i}'
                            if f'Group_{i}' not in self.group_mapping
                            else self.group_mapping[f'Group_{i}']
                        )
                    )
        if group:
            groups.append(group)
        return groups

    def get_unitname_groups(self):
        groups = self.get_groups()
        groups_list = []
        group = []
        for g in groups:
            group = [u.name for u in g.units if isinstance(u, StratigraphicUnit)]
            groups_list.append(group)
        return groups_list
    
    def get_group_unit_pairs(self) -> List[Tuple[str,str]]:
        """
        Returns a list of tuples containing group names and unit names.
        """
        groups = self.get_groups()
        group_unit_pairs = []
        for g in groups:
            for u in g.units:
                if isinstance(u, StratigraphicUnit):
                    group_unit_pairs.append((g.name, u.name))
        return group_unit_pairs

    def __getitem__(self, uuid):
        """
        Retrieves an element by its uuid from the stratigraphic column.
        """
        for element in self.order:
            if element.uuid == uuid:
                return element
        raise KeyError(f"No element found with uuid: {uuid}")

    def update_order(self, new_order):
        """
        Updates the order of elements in the stratigraphic column based on a new order list.
        """
        if not isinstance(new_order, list):
            raise TypeError("New order must be a list")
        self.order = [
            self.__getitem__(uuid) for uuid in new_order if self.__getitem__(uuid) is not None
        ]
        self.notify('order_updated', new_order=self.order)

    def update_element(self, unit_data: Dict):
        """
        Updates an existing element in the stratigraphic column with new data.
        :param unit_data: A dictionary containing the updated data for the element.
        """
        if not isinstance(unit_data, dict):
            raise TypeError("unit_data must be a dictionary")
        element = self.__getitem__(unit_data['uuid'])
        if isinstance(element, StratigraphicUnit):
            element.name = unit_data.get('name', element.name)
            element.colour = unit_data.get('colour', element.colour)
            element.thickness = unit_data.get('thickness', element.thickness)
        elif isinstance(element, StratigraphicUnconformity):
            element.name = unit_data.get('name', element.name)
            element.unconformity_type = UnconformityType(
                unit_data.get('unconformity_type', element.unconformity_type.value)
            )
        self.notify('element_updated', element=element)

    def __str__(self):
        """
        Returns a string representation of the stratigraphic column, listing all elements.
        """
        return "\n".join([f"{i+1}. {element}" for i, element in enumerate(self.order)])

    def to_dict(self):
        """
        Converts the stratigraphic column to a dictionary representation.
        """
        return {
            "elements": [element.to_dict() for element in self.order],
        }
    def update_from_dict(self, data):
        """
        Updates the stratigraphic column from a dictionary representation.
        """
        if not isinstance(data, dict):
            raise TypeError("Data must be a dictionary")
        with self.freeze_notifications():
            self.clear(basement=False)
            elements_data = data.get("elements", [])
            for element_data in elements_data:
                if "unconformity_type" in element_data:
                    element = StratigraphicUnconformity.from_dict(element_data)
                else:
                    element = StratigraphicUnit.from_dict(element_data)
                self.add_element(element)
    @classmethod
    def from_dict(cls, data):
        """
        Creates a StratigraphicColumn from a dictionary representation.
        """
        if not isinstance(data, dict):
            raise TypeError("Data must be a dictionary")
        column = cls()
        column.clear(basement=False)
        elements_data = data.get("elements", [])
        for element_data in elements_data:
            if "unconformity_type" in element_data:
                element = StratigraphicUnconformity.from_dict(element_data)
            else:
                element = StratigraphicUnit.from_dict(element_data)
            column.add_element(element)
        return column

    def get_isovalues(self) -> Dict[str, float]:
        """
        Returns a dictionary of isovalues for the stratigraphic units in the column.
        """
        surface_values = {}
        for g in reversed(self.get_groups()):
            v = 0
            for u in g.units:
                surface_values[u.name] = {'value':v,'group':g.name,'colour':u.colour}
                v += u.thickness
        return surface_values

    def plot(self,*, ax=None, **kwargs):
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        n_units = 0  # count how many discrete colours (number of stratigraphic units)
        xmin = 0
        ymin = 0
        ymax = 1
        xmax = 1
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=(2, 10))
        patches = []  # stores the individual stratigraphic unit polygons

        total_height = 0
        prev_coords = [0, 0]

        # iterate through groups, skipping faults
        for u in reversed(self.order):
            if u.element_type == StratigraphicColumnElementType.UNCONFORMITY:
                logger.info(f"Plotting unconformity {u.name} of type {u.unconformity_type.value}")
                ax.axhline(y=total_height, linestyle='--', color='black')
                ax.annotate(
                    getattr(u, 'name', 'Unconformity'),
                    xy=(xmin, total_height),
                    fontsize=8,
                    ha='left',
                )

                total_height -= 0.05  # Adjust height slightly for visual separation
                continue

            if u.element_type == StratigraphicColumnElementType.UNIT:
                logger.info(f"Plotting unit {u.name} of type {u.element_type}")

                n_units += 1

                ymax = total_height
                ymin = ymax - (getattr(u, 'thickness', np.nan) if not np.isinf(getattr(u, 'thickness', np.nan)) else np.nanmean([getattr(e, 'thickness', np.nan) for e in self.order if not np.isinf(getattr(e, 'thickness', np.nan))]))

                if not np.isfinite(ymin):
                    ymin = prev_coords[1] - (prev_coords[1] - prev_coords[0]) * (1 + rng.random())

                total_height = ymin

                prev_coords = (ymin, ymax)

                polygon_points = np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])
                patches.append(Polygon(polygon_points))
                ax.annotate(getattr(u, 'name', 'Unknown'), xy=(xmin+(xmax-xmin)/2, (ymax-ymin)/2+ymin), fontsize=8, ha='left')

        if 'cmap' not in kwargs:
            import matplotlib.colors as colors

            colours = []
            boundaries = []
            data = []
            for i, u in enumerate(self.order):
                if u.element_type != StratigraphicColumnElementType.UNIT:
                    continue
                data.append((i, u.colour))
                colours.append(u.colour)
                boundaries.append(i)  # print(u,v)
            cmap = colors.ListedColormap(colours)
        else:
            cmap = cm.get_cmap(kwargs['cmap'], n_units - 1)
        p = PatchCollection(patches, cmap=cmap)

        colors = np.arange(len(patches))
        p.set_array(np.array(colors))

        ax.add_collection(p)

        ax.set_ylim(total_height - (total_height - prev_coords[0]) * 0.1, 0)

        ax.axis("off")

        return fig
