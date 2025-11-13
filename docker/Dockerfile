# Usage: 
# git clone https://github.com/RosettaCommons/RFdiffusion.git
# cd RFdiffusion
# docker build -f docker/Dockerfile -t rfdiffusion .
# mkdir $HOME/inputs $HOME/outputs $HOME/models
# bash scripts/download_models.sh $HOME/models
# wget -P $HOME/inputs https://files.rcsb.org/view/5TPN.pdb

# docker run -it --rm --gpus all \
#   -v $HOME/models:$HOME/models \
#   -v $HOME/inputs:$HOME/inputs \
#   -v $HOME/outputs:$HOME/outputs \
#   rfdiffusion \
#   inference.output_prefix=$HOME/outputs/motifscaffolding \
#   inference.model_directory_path=$HOME/models \
#   inference.input_pdb=$HOME/inputs/5TPN.pdb \
#   inference.num_designs=3 \
#   'contigmap.contigs=[10-40/A163-181/10-40]'

FROM nvcr.io/nvidia/cuda:11.6.2-cudnn8-runtime-ubuntu20.04

# Install base utilities
RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# install mamba
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh -O ~/mamba.sh  && \
    /bin/bash  ~/mamba.sh -b
# Put conda in path so we can use conda activate
ENV PATH $CONDA_DIR/bin:$PATH

COPY env/SE3nv.yml /tmp/SE3nv.yml 

RUN  /root/mambaforge/bin/mamba  env create -f /tmp/SE3nv.yml  --prefix /opt/conda/envs/SE3nv

ENV PATH /opt/conda/envs/SE3nv/bin:$PATH

COPY . /app/RFdiffusion/

WORKDIR /app/RFdiffusion/env/SE3Transformer
RUN pip install --no-cache-dir -r requirements.txt
RUN python setup.py install 
RUN pip install /app/RFdiffusion

ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/conda/envs/SE3nv/lib/
CMD /app/RFdiffusion/scripts/run_inference.py 
