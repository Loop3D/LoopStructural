import geopandas
import pandas
import shapely
from .logging import getLogger
from typing import Any, Optional

logger = getLogger(__name__)

class ContactExtractor:
    def __init__(self, geology: geopandas.GeoDataFrame, faults: Optional[geopandas.GeoDataFrame] = None):
        self.geology = geology
        self.faults = faults
        self.contacts = None
        self.basal_contacts = None
        self.all_basal_contacts = None
    def __call__(self, **kwds: Any) -> Any:
        if 'stratigraphic_column' in kwds:
            return self.extract_basal_contacts(kwds['stratigraphic_column'], save_contacts=kwds.get('save_contacts', True))
        else:
            return self.extract_all_contacts(save_contacts=kwds.get('save_contacts', True))
    def extract_all_contacts(self, save_contacts: bool = True) -> geopandas.GeoDataFrame:
        logger.info("Extracting contacts")
        geology = self.geology.copy()
        geology = geology.dissolve(by="UNITNAME", as_index=False)
        if 'INTRUSIVE' in geology.columns:
            geology = geology[~geology["INTRUSIVE"]]
        if 'SILL' in geology.columns:
            geology = geology[~geology["SILL"]]
        if self.faults is not None:
            faults = self.faults.copy()
            faults["geometry"] = faults.buffer(50)
            geology = geopandas.overlay(geology, faults, how="difference", keep_geom_type=False)
        units = geology["UNITNAME"].unique().tolist()
        column_names = ["UNITNAME_1", "UNITNAME_2", "geometry"]
        contacts = geopandas.GeoDataFrame(crs=geology.crs, columns=column_names, data=None)
        while len(units) > 1:
            unit1 = units[0]
            units = units[1:]
            for unit2 in units:
                if unit1 != unit2:
                    join = geopandas.overlay(
                        geology[geology["UNITNAME"] == unit1],
                        geology[geology["UNITNAME"] == unit2],
                        keep_geom_type=False,
                    )[column_names]
                    join["geometry"] = join.buffer(1)
                    buffered = geology[geology["UNITNAME"] == unit2][["geometry"]].copy()
                    buffered["geometry"] = buffered.boundary
                    end = geopandas.overlay(buffered, join, keep_geom_type=False)
                    if len(end):
                        contacts = pandas.concat([contacts, end], ignore_index=True)
        contacts["length"] = [row.length for row in contacts["geometry"]]
        if save_contacts:
            self.contacts = contacts
        return contacts

    def extract_basal_contacts(self, 
                               stratigraphic_column: list, 
                               save_contacts: bool = True) -> geopandas.GeoDataFrame:
        
        logger.info("Extracting basal contacts")
        units = stratigraphic_column
        
        if self.contacts is None:
            self.extract_all_contacts(save_contacts=True)
            basal_contacts = self.contacts.copy()
        else:
            basal_contacts = self.contacts.copy()
        if any(unit not in units for unit in basal_contacts["UNITNAME_1"].unique()):
            missing_units = (
                basal_contacts[~basal_contacts["UNITNAME_1"].isin(units)]["UNITNAME_1"]
                .unique()
                .tolist()
            )
            logger.error(
                "There are units in the Geology dataset, but not in the stratigraphic column: "
                + ", ".join(missing_units)
                + ". Please readjust the stratigraphic column if this is a user defined column."
            )
            raise ValueError(
                "There are units in stratigraphic column, but not in the Geology dataset: "
                + ", ".join(missing_units)
                + ". Please readjust the stratigraphic column if this is a user defined column."
            )
        basal_contacts["ID"] = basal_contacts.apply(
            lambda row: min(units.index(row["UNITNAME_1"]), units.index(row["UNITNAME_2"])), axis=1
        )
        basal_contacts["basal_unit"] = basal_contacts.apply(lambda row: units[row["ID"]], axis=1)
        basal_contacts["stratigraphic_distance"] = basal_contacts.apply(
            lambda row: abs(units.index(row["UNITNAME_1"]) - units.index(row["UNITNAME_2"])), axis=1
        )
        basal_contacts["type"] = basal_contacts.apply(
            lambda row: "ABNORMAL" if abs(row["stratigraphic_distance"]) > 1 else "BASAL", axis=1
        )
        basal_contacts = basal_contacts[["ID", "basal_unit", "type", "geometry"]]
        basal_contacts["geometry"] = [
            shapely.line_merge(shapely.snap(geo, geo, 1)) for geo in basal_contacts["geometry"]
        ]
        if save_contacts:
            self.all_basal_contacts = basal_contacts
            self.basal_contacts = basal_contacts[basal_contacts["type"] == "BASAL"]
        return basal_contacts
