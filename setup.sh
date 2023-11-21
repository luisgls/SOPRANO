#!/usr/bin/env bash

ENV_NAME="soprano"
DEPS_FILE="src/SOPRANO/local.yml"
DEPS_FATAL=false

# Parse installation method for pip
if [ "$#" -eq 0 ]
then
  echo "Will install standard PyPI dependencies"
  INSTALL_METHOD=""
elif [ "$1" == "dev" ] || [ "$1" == "ci" ]
then
  echo "Will install $1 PyPI dependencies"
  INSTALL_METHOD="[$1]"
else
  echo "Will install dev PyPI dependencies"
  INSTALL_METHOD="[dev]"
fi
_PIP_CMD="pip install -e .$INSTALL_METHOD"

if [ ! -f "$DEPS_FILE" ]
then
  echo "Fatal: Conda dependencies file not found"
  DEPS_FATAL=true
else
  echo "Conda dependencies file detected: $DEPS_FILE"
fi

# Get os type
case "$OSTYPE" in
  darwin*)
    OS="osx" ;;
  linux*)
    OS="linux" ;;
  *)
    echo "OS $OSTYPE not implemented."
    OS="other" ;;
esac

# Flag whether to run update
built_from_scratch=false

function _create_for_osx() {
  echo "Creating env for osx"
  built_from_scratch=true
  conda create --name $ENV_NAME
  conda activate $ENV_NAME
  conda config --env --set subdir osx-64
  if command -v mamba &> /dev/null
  then
    echo "  ... building with mamba"
    mamba env update --file $DEPS_FILE
  else
    conda env update --file $DEPS_FILE
  fi
  conda activate $ENV_NAME
}

function _create_for_linux() {
  echo "Creating env for linux"
  built_from_scratch=true
  if command -v mamba &> /dev/null
  then
    echo "  ... building with mamba"
    mamba env create -f $DEPS_FILE --name $ENV_NAME
  else
    conda env create -f $DEPS_FILE --name $ENV_NAME
  fi
  conda activate $ENV_NAME
}

function _update() {
  if command -v mamba &> /dev/null
  then
    echo "  ... updating with mamba"
    mamba env update --file $DEPS_FILE --prune
  else
    conda env update --file $DEPS_FILE --prune
  fi
}

function _update_for_osx() {
  echo "Updating env for osx"
  conda config --env --set subdir osx-64
  _update
}

function _update_for_linux() {
  echo "Updating env for linux"
  _update
}

function activate_env() {
  echo "Attempting activation"
  if [ "$OS" == "osx" ]
  then
    echo "... osx"
    conda activate $ENV_NAME || _create_for_osx

    if [ $built_from_scratch == false ]
    then
      _update_for_osx
    fi
  elif [ "$OS" == "linux" ]
  then
    echo "... linux"
    conda activate $ENV_NAME || _create_for_linux

    if [ $built_from_scratch == false ]
    then
      _update_for_linux
    fi
  else
    exit 1
  fi
  if [ $? == 0 ]
  then
    echo "Installing SOPRANO Python packages"
    eval "$_PIP_CMD"
  fi
  if [ $? == 0 ]
  then
    echo "Installing SOPRANO R packages"
    Rscript "src/SOPRANO/R/pkgs.R"
  fi
}

if command -v conda &> /dev/null
then
  echo "Conda detected: Running environment setup"
  echo "OS system $OS: $OSTYPE"

  if [ $DEPS_FATAL == true ]
  then
    echo "    ... but unable to build due to missing dependencies file."
  else
    activate_env
    if [ $? == 0 ]
    then
      echo ""
      echo "Installation procedure complete."
      echo "After pulling an updated version of the repository"
      echo "run this setup script. Otherwise, call"
      echo "    conda activate $ENV_NAME"
      echo "to load the environment."
    else
      echo ""
      echo "Installation failed."
    fi
  fi

else
  echo "Fatal: Conda not detected."
fi

fasta_file="src/SOPRANO/data/ensemble_transcriptID.fasta"
zipped_fasta_file="src/SOPRANO/data/ensemble_transcriptID.fasta.gz"

if [ ! -f "$fasta_file" ]
then
  echo "Decompressing transcript IDs"
  gunzip -k "$zipped_fasta_file"
  echo "Decompressed: $fasta_file"
else
  echo "Decompressed transcript IDs detected: $fasta_file"
fi
