#!/bin/bash

# Keep track of project root
PROJECT_ROOT="$( cd -- "$(dirname -- "$0")" >/dev/null 2>&1 || exit ; pwd -P )"

# Installation helpers dir
SHELL_UTILS_DIR="src/SOPRANO/shell_utils"

# Data directory
DATA_DIR_PATH="src/SOPRANO/data"

if [ "$#" -eq 0 ]
then
  echo "-- installing standard dependencies"
  DEPS=""
elif [ "$1" == "dev" ] || [ "$1" == "ci" ]
then
  echo "-- installing dependencies for $1"
  DEPS="[$1]"
else
  echo "-- unrecognized dependencies argument. (defaulting to dev)"
  DEPS="[dev]"
fi
_PIP_CMD="pip install -e .$DEPS"

# Configure conda environment and install repository
source "$SHELL_UTILS_DIR/setup_env.sh"
source "$SHELL_UTILS_DIR/setup_data.sh"