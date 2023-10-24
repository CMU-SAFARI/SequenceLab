#!/bin/bash

# Get the directory of the current script
script_dir=$(dirname "$0")
# Change to the directory of the script using pushd
pushd "$script_dir" || exit 1

wget https://zenodo.org/records/10028978/files/ont.tar.gz
tar xf ont.tar.gz
rm ont.tar.gz

popd
