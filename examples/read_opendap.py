"""
Read a field over OPeNDAP
========================

Read the last stored time slice without downloading the entire dataset.
Try NOAA's public monthly sea-surface temperature dataset::

    python examples/read_opendap.py \
        https://psl.noaa.gov/thredds/dodsC/Datasets/noaa.ersst.v5/sst.mnmean.nc sst

You can substitute your own cGENIE dataset URL and variable name. The URL must
be an OPeNDAP data endpoint. This example accesses the network only when run
as a script.
"""

import argparse

from cgeniepy import GenieModel

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", help="OPeNDAP dataset URL")
    parser.add_argument("variable", help="NetCDF variable with a time dimension")
    args = parser.parse_args()

    with GenieModel(args.url) as model:
        field = model.get_var(args.variable).isel(time=-1)
        field.data.load()
        print(field)
