# Use an NVIDIA CUDA base image that includes CUDA 11.8 and is compatible with Ubuntu 20.04
FROM nvidia/cuda:11.8.0-devel-ubuntu20.04
# Set the working directory
WORKDIR /usr/src/app
# Install Python 3.8 and essential packages
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        software-properties-common \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        python3.8 python3.8-distutils python3.8-venv \
        curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install pip for Python 3.8
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.8
# Remove existing python links if they exist and create new ones
RUN if [ -f /usr/bin/python ]; then rm /usr/bin/python; fi && \
    if [ -f /usr/bin/python3 ]; then rm /usr/bin/python3; fi && \
    ln -s /usr/bin/python3.8 /usr/bin/python3 && \
    ln -s /usr/bin/python3.8 /usr/bin/python

# Install compatible versions of torch, pandas, scipy, numpy, tqdm
# These versions are selected based on compatibility with Python 3.8 and PyTorch 2.1.2
RUN pip install torch --index-url https://download.pytorch.org/whl/cu118
RUN pip install pandas scipy numpy tqdm
# Optional: Clean up pip cache to save space
RUN rm -rf /root/.cache/pip/*
