from .geological_model import GeologicalModel


class ProbabilisticGeologicalModel(GeologicalModel):
    def __init__(self, origin, maximum, rescale=True, nsteps=(40, 40, 40)):
        GeologicalModel.__init__(origin,maximum,rescale,nsteps)
