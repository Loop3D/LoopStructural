from LoopStructural.visualisation import model_visualisation as model_visualisation
def _get_loop_visualisation_scraper():
    return Scraper()

class Scraper:
    """
    Save ``pyvista.Plotter`` objects.

    Used by sphinx-gallery to generate the plots from the code in the examples.

    Pass an instance of this class to ``sphinx_gallery_conf`` in your
    ``conf.py`` as the ``"image_scrapers"`` argument.
    """

    def __call__(self, block, block_vars, gallery_conf):
        """Save the figures generated after running example code.

        Called by sphinx-gallery.

        """
        try:
            from sphinx_gallery.scrapers import figure_rst
        except ImportError:
            raise ImportError('You must install `sphinx_gallery`')
        image_names = list()
        image_path_iterator = block_vars["image_path_iterator"]
        figures = model_visualisation._OPEN_VIEWERS
        for address, plotter in figures.items():
            plotter.lv['background'] = 'white'
            fname = next(image_path_iterator)
            plotter.save(fname)
            image_names.append(fname)
        model_visualisation.close_all() # close and clear all plotters
        return figure_rst(image_names, gallery_conf["src_dir"])