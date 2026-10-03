from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from cgeniepy import EcoModel, GenieModel


@pytest.fixture
def remote_datasets(tmp_path, monkeypatch):
    """Serve real lazy NetCDF arrays through a mocked remote connection."""
    dataset = xr.Dataset(
        {
            "ocn_sur_temp": (
                ("time", "lat", "lon"),
                np.arange(12, dtype=float).reshape(2, 2, 3),
                {"long_name": "Surface temperature", "units": "degC"},
            )
        },
        coords={"time": [0, 1], "lat": [-45, 45], "lon": [0, 120, 240]},
    )
    dataset["ocn_sur_temp"][0, 0, 0] = np.nan
    path = tmp_path / "fields.nc"
    dataset.to_netcdf(path)
    open_dataset = xr.open_dataset
    opened = []
    closed = []

    def open_remote(url, **kwargs):
        ds = open_dataset(path, **kwargs)
        # Preserve the actual backend close callback before installing our spy.
        backend_close = ds._close

        def close_remote():
            closed.append(url)
            backend_close()

        ds.set_close(close_remote)
        opened.append((url, kwargs, ds))
        return ds

    monkeypatch.setattr(xr, "open_dataset", open_remote)
    return opened, closed


@pytest.mark.parametrize(
    "url",
    [
        "https://example.org/dodsC/run/fields_biogem_2d.nc",
        "http://example.org/dap/dataset?ocn_sur_temp[0:1:1][0:1:1][0:1:2]",
    ],
)
def test_remote_selection_and_close(url, remote_datasets):
    opened, closed = remote_datasets
    with GenieModel(url, decode_times=False) as model:
        assert model.is_remote
        assert not model.is_ensemble
        assert model.tsvar_list == []
        assert "OPeNDAP" in repr(model)
        assert "ocn_sur_temp" in model.ncvar_dict[url]
        assert not opened[0][2]["ocn_sur_temp"].variable._in_memory

        field = model.get_var("ocn_sur_temp")
        assert not field.data.variable._in_memory
        assert field.data.dims == ("time", "lat", "lon")
        assert field.attrs["long_name"] == "Surface temperature"
        assert field.attrs["units"] == "°C"
        assert field.isel(time=1).sel(lat=45, lon=240).data.item() == 11
        assert np.isnan(field.isel(time=0, lat=0, lon=0).data.item())
        assert len(opened) == 1
        assert opened[0][:2] == (url, {"engine": "netcdf4", "decode_times": False})
        with pytest.raises(ValueError, match="not found"):
            model.get_var("missing")

    assert closed == [url]
    assert model._datasets == {}
    model.close()
    assert closed == [url]


def test_remote_ensemble(remote_datasets):
    opened, closed = remote_datasets
    urls = ["https://example.org/dap/control", "https://example.org/dap/experiment"]
    with GenieModel(urls) as model:
        assert model.is_ensemble
        assert "ensemble of 2" in repr(model)
        field = model.get_var("ocn_sur_temp")
        assert field.data.dims == ("model", "time", "lat", "lon")
        assert field.data.model.values.tolist() == urls
        np.testing.assert_array_equal(
            field.isel(time=1, lat=1, lon=2).data.values, [11, 11]
        )
        model.get_var("ocn_sur_temp")
        assert len(opened) == 2
    assert closed == urls


def test_remote_reader_options(monkeypatch):
    calls = []
    session = object()

    def open_remote(url, **kwargs):
        calls.append((url, kwargs))
        return xr.Dataset({"temperature": ("time", [10.0])})

    monkeypatch.setattr(xr, "open_dataset", open_remote)
    url = "https://example.org/dap/dataset"
    with GenieModel(url, engine="pydap", backend_kwargs={"session": session}) as model:
        assert model.get_var("temperature").data.item() == 10
    assert calls == [(url, {"engine": "pydap", "backend_kwargs": {"session": session}})]


def test_remote_text_files_unavailable(remote_datasets):
    url = "https://example.org/dap/dataset"
    with GenieModel(url) as model:
        for operation, args in [
            (model.get_ts, ("ocn_temp",)),
            (model.get_diag_avg, (9999,)),
            (model.get_config, ()),
        ]:
            with pytest.raises(NotImplementedError, match="requires local cGENIE text"):
                operation(*args)
    with pytest.raises(NotImplementedError, match="Plankton_params.txt"):
        EcoModel(url)
    assert remote_datasets[1] == [url, url]


@pytest.mark.parametrize(
    "paths, message",
    [
        ([], "at least one"),
        (["https://example.org/dap/dataset", "."], "cannot mix"),
        ("https:///dataset", "server name"),
    ],
)
def test_invalid_remote_sources(paths, message):
    with pytest.raises(ValueError, match=message):
        GenieModel(paths)


def test_connection_error_propagates(monkeypatch):
    def fail(*args, **kwargs):
        raise OSError("OPeNDAP server unavailable")

    monkeypatch.setattr(xr, "open_dataset", fail)
    with pytest.raises(OSError, match="OPeNDAP server unavailable"):
        GenieModel("https://example.org/dap/dataset")


def test_local_paths_and_initialisation_cleanup(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    run = tmp_path / "run"
    (run / "biogem").mkdir(parents=True)
    for dim in ("2d", "3d"):
        xr.Dataset({f"field_{dim}": ("lat", [1.0, 2.0])}).to_netcdf(
            run / "biogem" / f"fields_biogem_{dim}.nc"
        )

    with GenieModel("~/run", gemflag="biogem") as model:
        assert model.model_path == str(run)
        assert model.get_var("field_3d").mean().data.item() == 1.5
    with GenieModel([run, run], gemflag="biogem") as model:
        field = model.get_var("field_2d")
        assert field.data.sizes["model"] == 2
        assert field.data.model.values.tolist() == [str(run), str(run)]

    opened = xr.open_dataset(run / "biogem" / "fields_biogem_2d.nc")
    closed = []
    backend_close = opened._close

    def close():
        closed.append(True)
        backend_close()

    opened.set_close(close)

    def fail_second(path, **kwargs):
        if Path(path).name.endswith("3d.nc"):
            raise OSError("cannot open second dataset")
        return opened

    monkeypatch.setattr(xr, "open_dataset", fail_second)
    with pytest.raises(OSError, match="second dataset"):
        GenieModel(run, gemflag="biogem")
    assert closed == [True]
