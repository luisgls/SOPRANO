#!/bin/bash

# Keep track of project root
PROJECT_ROOT="$( cd -- "$(dirname -- "$0")" >/dev/null 2>&1 || exit ; pwd -P )"

# Installation helpers dir
INSTALLERS_DIR_PATH="src/SOPRANO/bash_installers"

# Configure conda environment and install repository
source "$INSTALLERS_DIR_PATH/.setup_env.sh"