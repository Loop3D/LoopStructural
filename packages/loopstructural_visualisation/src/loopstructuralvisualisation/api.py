from ._3d_viewer import Loop3DView


def plot_block_model(model, filename=None, **kwargs):
    """
    Plot the model using pyvista

    Parameters
    ----------
    model : LoopStructuralModel
        The model to plot
    kwargs : dict
        Keyword arguments to pass to the plot
    """

    p = Loop3DView(model)
    p.plot_block_model(**kwargs)
    if filename is not None:
        p.screenshot(filename)
    p.show()


def plot_surface(model, geological_feature, **kwargs):
    """
    Plot a surface using pyvista

    Parameters
    ----------
    model : LoopStructuralModel
        The model to plot
    geological_feature : BaseFeature
        The feature to plot
    kwargs : dict
        Keyword arguments to pass to the plot
    """

    p = Loop3DView(model)
    p.plot_surface(geological_feature, **kwargs)
    return p
