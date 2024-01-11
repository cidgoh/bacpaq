FROM --platform=linux/amd64 mambaorg/micromamba:0.14.0

# metadata
LABEL image.authors="CIDGOH"
LABEL image.version="0.1"
LABEL image.peppan="1.0.5"
LABEL organization="CIDGOH"

RUN apt-get update && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/* && \
    CONDA_DIR="/opt/conda" && \
    micromamba install -n base -y -c bioconda -c conda-forge mmseqs2 blast rapidnj diamond fasttree python=3.9.1 && \
    pip install --no-cache peppan && \
    micromamba clean --all --yes && \
    rm -rf $CONDA_DIR/conda-meta && \
    rm -rf $CONDA_DIR/include && \
    rm -rf $CONDA_DIR/lib/python3.*/site-packages/pip && \
    find $CONDA_DIR -name '__pycache__' -type d -exec rm -rf '{}' '+'

# ARG MAMBA_DOCKERFILE_ACTIVATE=1
# RUN python3 -m pip install --no-cache-dir peppan
# RUN pip install --no-cache peppan