#!/bin/bash

# Mamba/conda path
_mamba_path=$(which mamba)
_conda_path=$(which conda)
_env_grep=$(conda env list | grep "soprano-dev")

yml_path="src/SOPRANO/local.yml"

function deactivate_all() {
  for _ in $(seq "${CONDA_SHLVL}")
  do
    conda deactivate
  done
}

if [ -n "$_mamba_path" ]
then
	echo "-- mamba detected"
  deactivate_all

	if [ -z "$_env_grep" ]
	then
		echo "-- building environment"
		mamba env create -f $yml_path
    echo "-- activating soprano"
    mamba activate soprano-dev
	else
		echo "-- soprano environment detected"
    echo "-- activating soprano"
    mamba activate soprano-dev
		echo "-- checking for updates to environment"
		mamba env update --file $yml_path --prune
	fi

elif [ -n "$_conda_path" ]
then
	echo "-- conda detected"
  deactivate_all
	if [ -z "$_env_grep" ]
	then
		echo "-- building conda environment"
		conda env create -f src/SOPRANO/local.yml
    echo "-- activating soprano"
    conda activate soprano-dev
	else
		echo "-- soprano environment detected"
    echo "-- activating soprano"
    conda activate soprano-dev
		echo "-- checking for updates to environment"
		conda env update --file $yml_path --prune
	fi

else
	echo "No conda variant available to build environment!"
fi

eval "$_PIP_CMD"