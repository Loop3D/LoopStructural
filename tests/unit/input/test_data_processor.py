import numpy as np
import pandas as pd

from LoopStructural.modelling import ProcessInputData
from LoopStructural.utils import rng


def test_create_processor():
    df = pd.DataFrame(rng.random(size=(10, 3)), columns=["X", "Y", "Z"])
    df["name"] = [f"unit_{name % 2}" for name in range(10)]
    stratigraphic_order = [("sg", ["unit_0", "unit_1", "basement"])]
    thicknesses = {"unit_0": 1.0, "unit_1": 0.5}
    processor = ProcessInputData(
        contacts=df, stratigraphic_order=stratigraphic_order, thicknesses=thicknesses
    )
    assert (processor.data["val"].unique() == np.array([0.5, 0])).all()


def test_map2loop_processor_gives_each_supergroup_its_own_name(tmp_path, monkeypatch):
    import networkx

    from LoopStructural.modelling.input import map2loop_processor as m2lp

    (tmp_path / "tmp").mkdir()
    (tmp_path / "output").mkdir()
    files = {
        "tmp/all_sorts_clean.csv": ",group,code,colour\n0,g1,u1,red\n1,g2,u2,blue\n",
        "tmp/super_groups.csv": "g1,\ng2,\n",
        "tmp/bbox.csv": "minx,miny,maxx,maxy,lower,upper\n0,0,1,1,0,1\n",
        "output/orientations_clean.csv": "formation,polarity\nu1,1\n",
        "output/formation_summary_thicknesses.csv": "formation,thickness median\nu1,1\n",
        "output/contacts_clean.csv": "formation\nu1\n",
        "output/fault_displacements3.csv": "fname,vertical_displacement,downthrow_dir\n",
        "output/fault_orientations.csv": "formation,DipDirection\n",
        "output/faults.csv": "formation\n",
        "output/supergroup-fault-relationships.csv": "supergroup\n",
        "output/fault_dimensions.csv": "Fault,InfluenceDistance,VerticalRadius,HorizontalRadius\n",
    }
    for path, text in files.items():
        (tmp_path / path).write_text(text)
    networkx.write_gml(networkx.Graph(), tmp_path / "tmp" / "fault_network.gml")

    captured = {}

    def fake_init(self, contacts, orientations, stratigraphic_order, **kwargs):
        captured["order"] = stratigraphic_order

    monkeypatch.setattr(ProcessInputData, "__init__", fake_init)
    m2lp.Map2LoopProcessor(str(tmp_path))
    assert captured["order"] == [("supergroup_0", ["u1"]), ("supergroup_1", ["u2"])]
