from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_laurent2016

def average_axis():
    data, bb = load_laurent2016()

    model = GeologicalModel(bb[0,:],bb[1,:])
    model.set_model_data(data)
    s2 = model.create_and_add_fold_frame('s2',
                                         nelements=10000)

    s1 = model.create_and_add_folded_fold_frame('s1',
                                                limb_wl=.4,
                                                av_fold_axis=True,
                                                nelements=50000
                                               )

    s0 = model.create_and_add_folded_fold_frame('s0',
                                                limb_wl=1.,
                                                av_fold_axis=True,
                                                nelements=50000
                                               )