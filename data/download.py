"""
Simple script to download all .trees files within the `data` directory on GitHub,
saving to a local `data` directory
"""

import os
import json
from urllib.request import urlretrieve, urlopen

if not os.path.isdir("data"):
    os.mkdir("data")  # Make a "data" directory within the current folder
print(f"Downloading data files into {os.path.join(os.getcwd(), 'data')}")
# Save the data files to the data directory
response = urlopen("https://tskit.dev/tutorials/examples/files.txt")
for fn in response:
    fn = fn.decode(response.headers.get_content_charset()).strip()
    if fn.endswith(".trees"):
        urlretrieve("https://tskit.dev/tutorials/examples/" + fn, os.path.join("data", fn))
        print(".", end="")
print(" finished downloading")