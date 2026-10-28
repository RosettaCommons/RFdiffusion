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

# Install packages with uv rather than pip: no pip release that supports Python 3.9
# includes the fix for the October 2026 emSign CA incident. uv 0.13.4+ does.
COPY --from=ghcr.io/astral-sh/uv:0.13.6 /uv /uvx /bin/

COPY . /app/RFdiffusion/

RUN apt-get -q update \
  && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
  ca-certificates \
  git \
  python3.9 \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get autoremove -y \
  && apt-get clean \
  && uv pip install --system --python python3.9 -q --no-cache \
  'numpy<2' \
  'certifi>=2026.10.19' \
  dgl==1.0.2+cu116 -f https://data.dgl.ai/wheels/cu116/repo.html \
  torch==1.12.1+cu116 --extra-index-url https://download.pytorch.org/whl/cu116 \
  e3nn==0.3.3 \
  wandb==0.12.0 \
  pynvml==11.0.0 \
  git+https://github.com/NVIDIA/dllogger#egg=dllogger \
  decorator==5.1.0 \
  hydra-core==1.3.2 \
  pyrsistent==0.19.3 \
  /app/RFdiffusion/env/SE3Transformer \
  && uv pip install --system --python python3.9 --no-cache --no-deps /app/RFdiffusion
  
WORKDIR /app/RFdiffusion

ENV DGLBACKEND="pytorch"

ENTRYPOINT ["python3.9", "scripts/run_inference.py"]
