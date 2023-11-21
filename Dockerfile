FROM condaforge/mambaforge:23.3.1-1 as conda

WORKDIR /app

COPY pyproject.toml ./
COPY install_R_pkgs.R ./
COPY src ./src
RUN touch README.md

# Construct conda environment and clean up cache files
RUN mamba env create -f src/SOPRANO/local.yml --name soprano
RUN mamba clean -afy

# Install R dependencies pulled from GH
RUN ["mamba", "run", "--no-capture-output", "-n", "soprano", "Rscript", "install_R_pkgs.R"]

# Tweak TOML to support pip installation (basically omit hatch version bits)
RUN sed -i '71,84d' pyproject.toml
RUN sed -i '38d' pyproject.toml
RUN sed -i '/name = "SOPRANO"/a version="0.0.1"' pyproject.toml

# Run pip install
RUN ["mamba", "run", "--no-capture-output", "-n", "soprano", "pip", "install", "-e", "."]

# Expose port for streamlit interface
EXPOSE 8501

# Run health check
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

# Define streamlit run as entry point
ENTRYPOINT ["mamba", "run", "-n", "soprano", "streamlit", "run", "./src/SOPRANO/app.py", "--server.port=8501", "--server.address=0.0.0.0"] # , ">", "dev/null"]