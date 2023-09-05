#!/bin/bash

# Keep track of project root
PROJECT_ROOT="$( cd -- "$(dirname -- "$0")" >/dev/null 2>&1 || exit ; pwd -P )"

# Installation helpers dir
INSTALLERS_DIR_PATH="src/SOPRANO/bash_installers"

# Data directory
DATA_DIR_PATH="src/SOPRANO/data"

if [ "$#" -eq 0 ]
then
  echo "-- installing dependencies for dev"
  DEPS="dev"
elif [ "$1" == "dev" ] || [ "$1" == "ci" ]
then
  echo "-- installing dependencies for $1"
  DEPS=$1
else
  echo "-- unrecognized dependencies argument. (defaulting to dev)"
  DEPS="dev"
fi
_PIP_CMD="pip install -e .[$DEPS]"

# Configure conda environment and install repository
source "$INSTALLERS_DIR_PATH/.setup_env.sh"
source "$INSTALLERS_DIR_PATH/.setup_data.sh"