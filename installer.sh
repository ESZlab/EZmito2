#!/bin/bash

echo -e "EZmito2 installer. This script allows you to create the EZmito2 conda environment and install all the softwares...\n"
sleep 2

# Check if Conda is installed
conda=$(which conda)

# Make EZmito executable:
chmod 775 ezmito_env.yml
chmod 775 ezmito.py

# Check if Conda is installed
if [[ -z $conda ]]; then
    echo -e "Conda program not found. Please install conda to proceed\n"
    echo -e "Do you want us to install Conda for you? [y/n]\n"
    read yesorno
    if [[ "${yesorno,,}" == "n" ]]; then
        echo -e "\nPlease install Conda to proceed"
        exit 1
    elif [[ "${yesorno,,}" == "y" ]]; then
        echo -e "\nConda installation...\n"
        # Download the installer
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -P /home/$(whoami)
        # Make it executable
        chmod 775 /home/$(whoami)/Miniconda3-latest-Linux-x86_64.sh
        # Run the installer
        /home/$(whoami)/Miniconda3-latest-Linux-x86_64.sh
        echo -e "Creating EZmito conda environment\n"
        sleep 3
        # Create conda environment
        conda env create -f ezmito_env.yml
        condaenv=$(conda env list | grep "ezmito" | wc -l)
        if [[ $condaenv -gt 0 ]]; then
            echo -e "\nEZmito2 environment installed successfully\n"
        else 
            echo "Error in the environment creation"
            exit 1
        fi
    else
        echo -e "\nOnly 'y' or 'n' are accepted as answers."
        exit 1
    fi
    
WHICH_CONDA=$(conda info | grep -i 'active env location' | cut -d':' -f 2 | sed 's/ //g')
source ${WHICH_CONDA}"/etc/profile.d/conda.sh"


elif [[ -n $conda ]]; then
    echo -e "Conda found"
    condaenv=$(conda env list | grep "ezmito_env" | wc -l)
    if [[ $condaenv -gt 0 ]]; then
    	conda env remove -n ezmito_env
    fi
    # Create conda environment
    conda env create -f ezmito_env.yml
    condaenv=$(conda env list | grep "ezmito" | wc -l)
    if [[ $condaenv -gt 0 ]]; then
    	echo -e "\nDone\n"
    else
    	echo "Error in the environment creation"
    exit 1
       fi
fi

