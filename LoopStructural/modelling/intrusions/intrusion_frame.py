from ..features import StructuralFrame


class IntrusionFrame(StructuralFrame):
    """A StructuralFrame built specifically to parameterise an intrusion's
    curvilinear coordinate system, so that IntrusionFrameBuilder produces a
    type distinguishable from a generic StructuralFrame (mirroring how
    FaultBuilder produces a FaultSegment).
    """

