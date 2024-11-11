FROM ubuntu:18.04

LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at"

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# Install basic dependencies
RUN apt-get update && \
    apt-get install -y wget && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    mkdir /root/.conda && \
    bash Miniconda3-latest-Linux-x86_64.sh -b && \
    rm -f Miniconda3-latest-Linux-x86_64.sh

# Set up Miniconda environment
COPY environment.yml .
RUN conda env create -f environment.yml && \
    conda clean -a

# Activate the conda environment
ENV PATH="/root/miniconda3/envs/haplotyping_env/bin:${PATH}"
ARG PATH="/root/miniconda3/envs/haplotyping_env/bin:${PATH}"

# Set the default command to activate the conda environment
CMD ["/bin/bash", "-c", "source activate haplotyping_env && exec bash"]
