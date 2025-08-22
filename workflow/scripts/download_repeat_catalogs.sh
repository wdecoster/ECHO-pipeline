#!/bin/bash

# Download repeat catalog zip file from Zenodo to current directory
wget -O repeat_catalogs.zip "https://zenodo.org/record/16925640/files/repeat_catalogs.zip?download=1"

# Unzip - this will create the 'repeat_catalogs/' folder in the current directory
unzip -o repeat_catalogs.zip

echo "Repeat catalogs downloaded and extracted successfully"
