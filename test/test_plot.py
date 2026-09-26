import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
from cgeniepy.array import GriddedData
from cgeniepy.table import ScatterData
import cgeniepy

def create_sample_data():
    model = cgeniepy.sample_model()
    return model.get_var('ocn_sur_temp').isel(time=-1)

def test_map_creation():
    """Test that map plotting doesn't crash and creates expected elements"""
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Mollweide()})
    data = create_sample_data()
    
    # Test that plotting works without error
    data.plot(ax=ax, cmap="viridis")
    
    # Test that the plot has expected properties
    assert len(ax.collections) > 0  # Has plot elements
    
    plt.close(fig)

def test_line_plot():
    """Test line plotting functionality"""
    data = create_sample_data()
    fig, ax = plt.subplots()
    
    line_data = data.mean(dim='lon')
    line_data.plot(ax=ax)
    
    # Test plot properties
    assert len(ax.lines) > 0  # Has line elements
    assert ax.get_xlabel() != ""  # Has x-label
    assert ax.get_ylabel() != ""  # Has y-label
    
    plt.close(fig)


def test_contourf_uses_requested_colormap():
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Mollweide()})
    data = create_sample_data()

    contour_set = data.plot(
        ax=ax, pcolormesh=False, contourf=True, cmap="plasma"
    )

    assert contour_set.cmap.name == "plasma"
    plt.close(fig)


def test_cell_edges_follow_genie_grid():
    """Cells sit on the edges written by GENIE, for tracer and MOC grids."""
    model = cgeniepy.sample_model()
    temp = model.get_var("ocn_temp").isel(time=-1)
    opsi = model.get_var("phys_opsi").isel(time=-1)
    lat_edges = model.grid_lat_edges().data.values
    zt_edges = model.grid_zt_edges().data.values
    lat, zt = temp.data.lat.values, temp.data.zt.values

    *_, lat_cells, depth_cells = (
        temp.mean(dim="lon").to_GriddedDataVis()._transect_coordinates()
    )
    np.testing.assert_allclose(lat_cells, lat_edges, atol=1e-9)
    np.testing.assert_allclose(depth_cells, zt_edges, atol=2)

    # MOC values are nodes: their cells split at tracer centres and stop at the
    # poles, the surface and the seafloor
    *_, lat_cells, depth_cells = opsi.to_GriddedDataVis()._transect_coordinates()
    np.testing.assert_allclose(lat_cells, np.r_[-90, lat, 90], atol=1e-9)
    np.testing.assert_allclose(depth_cells, np.r_[0, zt, zt_edges[-1]], atol=2)


def test_cell_edges_for_subsets_and_zonal_sections():
    model = cgeniepy.sample_model()
    temp = model.get_var("ocn_temp").isel(time=-1)
    lon_edges = model.grid_lon_edges().data.values

    tropics = temp.isel(zt=0).sel(lat=slice(-30, 30))
    *_, lat_cells = tropics.to_GriddedDataVis()._map_coordinates()
    np.testing.assert_allclose(lat_cells[[0, -1]], [-30, 30], atol=1e-9)

    section = temp.sel(lat=0, method="nearest")
    name, *_, lon_cells, _ = section.to_GriddedDataVis()._transect_coordinates()
    assert name == "lon"
    np.testing.assert_allclose(lon_cells, lon_edges)


def test_contour_labels_can_be_chosen_or_turned_off():
    data = create_sample_data()
    fig, axes = plt.subplots(1, 2, subplot_kw={"projection": ccrs.PlateCarree()})

    vis = data.to_GriddedDataVis()
    vis.aes_dict["contour_kwargs"]["levels"] = [0, 10, 20]
    vis.aes_dict["contour_label_kwargs"]["levels"] = [10]
    labelled = vis.plot(ax=axes[0], pcolormesh=False, contour=True)
    unlabelled = data.plot(
        ax=axes[1], pcolormesh=False, contour=True, contour_label=False
    )

    assert list(labelled.labelLevelList) == [10]
    assert not unlabelled.labelTexts
    plt.close(fig)


def test_filled_contours_are_returned_with_line_contours():
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    contours = create_sample_data().plot(
        ax=ax, pcolormesh=False, contourf=True, contour=True, colorbar=True
    )
    assert contours.filled
    assert contours.colorbar is not None
    plt.close(fig)


def test_map_contours_close_across_the_seam():
    vis = create_sample_data().to_GriddedDataVis()
    longitude, _, longitude_edges, _ = vis._map_coordinates()
    closed_longitude, closed = vis._close_longitude(longitude, longitude_edges)

    assert closed_longitude[-1] == longitude[0] + 360
    np.testing.assert_array_equal(closed[:, -1], closed[:, 0])


def test_cells_straddling_the_map_edge_are_split_in_place():
    """Issue #4: a 0-359° grid on EckertIV used to crash or lose its facecolor."""
    lat = np.linspace(-89.5, 89.5, 180)
    lon = np.arange(360.0)
    values = np.random.default_rng(4).random((lat.size, lon.size))
    values[50:100, 50:100] = np.nan
    array = xr.DataArray(values, coords=[("lat", lat), ("lon", lon)])
    vis = GriddedData(array).to_GriddedDataVis()
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EckertIV()})

    _, _, longitude_edges, _ = vis._map_coordinates()
    edges, split = vis._split_at_map_edge(ax, longitude_edges)
    vis.plot(ax=ax)
    fig.canvas.draw()

    assert longitude_edges[0] == -0.5  # not nudged half a cell east
    assert np.setdiff1d(edges, longitude_edges).tolist() == [180.0]
    np.testing.assert_array_equal(split[:, 180], split[:, 181])
    plt.close(fig)


def test_transect_grid_font_and_colorbar_options():
    data = create_sample_data()
    model = cgeniepy.sample_model()
    zonal = model.get_var("ocn_temp").isel(time=-1).mean(dim="lon")
    fig, axes = plt.subplots(1, 2)

    default = zonal.plot(ax=axes[0])
    vis = zonal.to_GriddedDataVis()
    vis.aes_dict["general_kwargs"]["font"] = "sans-serif"
    vis.aes_dict["colorbar_kwargs"]["orientation"] = "horizontal"
    custom = vis.plot(ax=axes[1], gridline=False)

    assert default.colorbar.orientation == "vertical"
    assert any(line.get_visible() for line in axes[0].get_xgridlines())
    assert custom.colorbar.orientation == "horizontal"
    assert not any(line.get_visible() for line in axes[1].get_xgridlines())
    assert axes[1].yaxis.label.get_family() == ["sans-serif"]
    assert data.to_GriddedDataVis().aes_dict["general_kwargs"]["font"] is None
    plt.close(fig)

def test_scatterdata_plot():
    """Test scatter data visualization"""
    
    from importlib.resources import files
    ## access the sample data file
    ## but convert to str for the ScatterData class
    file_path = str(files('cgeniepy.data').joinpath('EDC_CO2.tab'))

    data = ScatterData(file_path, sep='\t')
    data.set_index('Age [ka BP]')
    
    fig, ax = plt.subplots()
    data.plot(var='CO2 [ppmv]', ax=ax)
    
    # Test that data was plotted
    assert len(ax.lines) > 0 or len(ax.collections) > 0
    
    plt.close(fig)
