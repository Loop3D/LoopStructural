from LoopStructural.utils import rng


def random_colour(n: int, cmap='tab20'):
    """
    Generate a list of random colours

    Parameters
    ----------
    n : int
        Number of colours to generate
    cmap : str, optional
        Name of the matplotlib colour map to use, by default 'tab20'

    Returns
    -------
    list
        List of colours in the form of (r,g,b,a) tuples
    """
    import matplotlib.cm as cm

    colours = []
    for _i in range(n):
        colours.append(cm.get_cmap(cmap)(rng.random()))

    return colours
