#!/bin/bash

# Get the directory of the current script
script_dir=$(dirname "$0")
# Change to the directory of the script using pushd
pushd "$script_dir" || exit 1

wget https://zenodo.org/records/10028978/files/illumina.tar.gz
tar xf illumina.tar.gz
rm illumina.tar.gz

popd
