#!/usr/bin/env bash

IS_CONDA="$(which conda)"

if [ ${#IS_CONDA} -eq "0" ]
    then
        echo "It seems that you do not have conda installed. Download it from https://www.continuum.io/downloads"
        exit 1
fi

echo "Removing preexisting dada2 environment if exists"
source deactivate
conda remove --name dada2v1.10 --all --yes;

echo "Creating the new dada2 environment from environment.yaml"
conda env create -f environment.yaml;

if [ $? -ne 0 ]
    then
        echo "Problem creating the conda environment. I don't know what to tell you."
        exit 1
fi

################################################################################
echo "------------------------------------------------------------------------------------------------------"
echo "INSTALLATION COMPLETE"
echo "------------------------------------------------------------------------------------------------------"
echo "NOTE: You must specify the required information in the config.yaml file to run the workflow properly."
echo "NOTE: Once configured, enter:"
echo ""
echo "> source activate dada2v1.6"
echo "> snakemake"
echo ""
echo "In this directory to run the workflow."
echo "------------------------------------------------------------------------------------------------------"
