def create_dtm_with_rasterio(dtm_path):
    try:

        from map2loop.map import MapUtil
    except ImportError:
        print("map2loop not installed. Please install it and try again")
    try:
        import rasterio
    except ImportError:
        print("rasterio not installed. Please install it and try again.")
        return
    dtm_map = MapUtil(None, dtm=rasterio.open(dtm_path))
    return lambda xyz: dtm_map.evaluate_dtm_at_points(xyz[:, :2])
