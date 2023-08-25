#!/bin/bash

VEP_DIR="$TOOLS_DIR_PATH/ensembl-vep"
VEP_CHECK_FILE="$TOOLS_DIR_PATH/vep_configured"
VEP_VERSION=110

if [ ! -f "$VEP_CHECK_FILE" ]
then
    # Initialise submodule
    if [ ! -d "$VEP_DIR" ]
    then
      git submodule init
    fi

    # Fetch pinned version
    cd "$VEP_DIR" || exit

    # Checkout pinned version of vep and install
    echo "-- configuring for vep version: $VEP_VERSION"
    git fetch --all --tags --prune
    git checkout release/$VEP_VERSION
    perl INSTALL.pl
    vep_install_status=$?

    if [ $vep_install_status -eq 0 ]
    then
      touch "$VEP_CHECK_FILE"
      echo "-- vep configured successfully"
    else
      echo "-- vep installation had non-zero exit status: $vep_install_status"
    fi

    cd "$PROJECT_ROOT" || exit
else
    echo "-- vep version $VEP_VERSION already configured"
fi