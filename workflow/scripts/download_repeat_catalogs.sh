#!/bin/bash

# Download repeat catalog zip file from Zenodo to current directory
wget -O repeat-catalogs.zip "https://zenodo.org/record/16925640/files/repeat-catalogs.zip?download=1"

# Unzip - this will create the 'repeat_catalogs/' folder in the current directory
unzip -o repeat-catalogs.zip

# Remove zip file after extraction
rm repeat-catalogs.zip

echo "Repeat catalogs downloaded and extracted successfully"
