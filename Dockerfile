FROM condaforge/mambaforge:23.3.1-1 as conda

WORKDIR /app

# Construct conda environment and clean up cache files
COPY src/SOPRANO/local.yml ./
RUN mamba env create -f local.yml --name soprano &&  \
    mamba clean -afy &&  \
    rm local.yml

# Install R dependencies pulled from GH
COPY install_R_pkgs.R ./
RUN mamba run --no-capture-output -n soprano Rscript install_R_pkgs.R &&  \
    rm install_R_pkgs.R

# Prepare Python dependencies (this could be better...)
COPY pyproject.toml ./
RUN sed -i '71,84d' pyproject.toml &&  \
    sed -i '38d' pyproject.toml &&  \
    sed -i '/name = "SOPRANO"/a version="0.0.1"' pyproject.toml && \
    touch README.md

# Run pip install
COPY src ./src
RUN mamba run --no-capture-output -n soprano pip install -e . &&  \
    mamba run --no-capture-output -n soprano pip cache purge

# Decompress transcript ID
RUN if [[ ! -f "./src/SOPRANO/data/ensemble_transcriptID.fasta" ]] ; then gunzip -v ./src/SOPRANO/data/ensemble_transcriptID.fasta.gz ; fi

# Clean up additional stuff ...
RUN rm -rf ./src/SOPRANO/immunopeptidomes/mouse && \
    find ./src/SOPRANO/scripts -name "*.R" -delete && \
    find -name "*.yml" -delete &&  \
    find -name "*.gz" -delete &&  \
    find -name '__pycache__' -type d -exec rm -rf {} + && \
    cd /opt/conda &&  \
    rm -rf conda-meta && \
    find -name '__pycache__' -type d -exec rm -rf {} + && \
    cd /opt/conda/lib && \
    find -name 'tests' -type d -exec rm -rf {} +

# Expose port for streamlit interface
EXPOSE 8501

# Define streamlit run as entry point
ENTRYPOINT ["mamba", "run", "-n", "soprano", "streamlit", "run", "./src/SOPRANO/app.py", "--server.port=8501", "--server.address=0.0.0.0"]

# Bind mount to downloads with:
# docker run -d -p 8501:8501 --name devtest -v "$(pwd)"/ensembl_downloads/homo_sapiens:/app/ensembl_downloads/homo_sapiens soprano
# docker run -d -p 8501:8501 --name devtest -e SOPRANO_DISABLE_DOWNLOADS=True -v "$(pwd)"/ensembl_downloads/homo_sapiens:/app/ensembl_downloads/homo_sapiens soprano
# Should be visible on e.g. http://localhost:8501/