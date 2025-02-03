from LoopStructural.utils import rng


def random_colour(n: int = 1, cmap='tab20'):
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
    from matplotlib import colormaps as cm

    colours = []
    for _i in range(n):
        colours.append(cm.get_cmap(cmap)(rng.random()))

    return colours


def random_hex_colour(n: int = 1, cmap='tab20'):
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
        List of colours in the form of hex strings
    """
    from matplotlib import colormaps as cm

    colours = []
    for _i in range(n):
        colours.append(cm.get_cmap(cmap)(rng.random()))

    return [f'#{int(c[0]*255):02x}{int(c[1]*255):02x}{int(c[2]*255):02x}' for c in colours]
