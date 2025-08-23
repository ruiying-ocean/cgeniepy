"""
=========================================
Plot the 2D ScatterData
=========================================

This example shows how to plot the 2D ScatterData object. I use a LGM d13C data from Peterson et al. 2014 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013PA002552) as an example.
"""

from cgeniepy.table import ScatterData
from cgeniepy.plot import CommunityPalette
import matplotlib.pyplot as plt
import requests
import os

def download_zenodo_file(record_id, filename, download_path="~/.cgeniepy/"):
    """
    Downloads a specific file from a Zenodo record.

    Args:
        record_id (str): The Zenodo record ID (the numeric part of the DOI).
        filename (str): The name of the file to download.
        download_path (str): The directory to save the file in.
    """
    api_url = f"https://zenodo.org/api/records/{record_id}"
    response = requests.get(api_url)
    response.raise_for_status()  # Raise an exception for bad status codes

    record_data = response.json()
    file_to_download = None
    for f in record_data.get('files', []):
        if f['key'] == filename:
            file_to_download = f
            break

    if not file_to_download:
        raise FileNotFoundError(f"File '{filename}' not found in Zenodo record '{record_id}'.")

    download_url = file_to_download['links']['self']
    file_path = os.path.join(os.path.expanduser(download_path), filename)

    print(f"Downloading {filename} to {file_path}...")
    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        with open(file_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print("Download complete.")
    return file_path

record_id = "13786013"
filename = "LGM_d13c_CLEAN.txt"

local_file_path = download_zenodo_file(record_id, filename)

proxy_d13C = ScatterData(local_file_path, delimiter="\t", header=None)
proxy_d13C.data.columns = ["Longitude", "Latitude", "Depth", "d13C","Event"]
proxy_d13C.set_index(["Latitude", "Longitude"])
cmap = CommunityPalette("BuDaRd18").colormap
proxy_d13C.plot(var='d13C', edgecolor='k', cmap=cmap)

plt.show()
