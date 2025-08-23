import os
from pathlib import Path
import numpy as np
import requests


def check_rm(path):
    """
    Check and delete file

    :param path: path string of target file
    """
    if os.path.isfile(path):
        print(f"removing {path}")
        os.remove(path)


def file_exists(path):
   """
   Check if file exists. If not, raise FileNotFoundError.

   :param path: path string of target file
   :returns: True if file exists
   """
   file_path = Path(path).expanduser()
   if file_path.is_file():
       return True
   else:
       raise FileNotFoundError(f"{path} not exist")


def is_empty(path):
    """
    Check if file is empty. If the file does not exist, then create one.

    :param path: path string of target file
    :returns: boolean operator
    """
    if os.path.isfile(path):
        return os.stat(path).st_size == 0
    else:
        Path(path).touch()

def efficient_log(data, replace_zero=10):
    """A shortcut to efficiently calculate log10 of data.

    :param data: data to calculate log10
    :param replace_zero: value to replace zero in data

    :returns: log10 of data
    """
    return np.where(data == 0, replace_zero, np.log10(data))


def download_zenodo_file(record_id, filename, download_path="~/.cgeniepy/"):
    """
    Downloads a specific file from a Zenodo record.
    
    Args:
        record_id (str): The Zenodo record ID (the numeric part of the DOI).
        filename (str): The name of the file to download.
        download_path (str): The directory to save the file in.
    """
    import requests
    from pathlib import Path
    
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
    
    # Create Path object and ensure directory exists
    download_dir = Path(download_path).expanduser()
    download_dir.mkdir(parents=True, exist_ok=True)
    
    file_path = download_dir / filename
    
    print(f"Downloading {filename} to {file_path}...")
    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        with file_path.open('wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    
    print("Download complete.")
    return str(file_path)
