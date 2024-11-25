FROM ubuntu:20.04

LABEL org.opencontainers.image.version="20.04"
LABEL org.opencontainers.image.ref.name="ubuntu"

ARG LAUNCHPAD_BUILD_ARCH
ARG RELEASE

RUN dpkg --add-architecture i386 && \
    apt-get update && \
    apt-get install -y \
        libc6:i386 \
        libncurses5:i386 \
        libstdc++6:i386 \
        zlib1g:i386 \
        gcc-multilib \
        g++-multilib \
        binutils \
        git \
        curl \
        libxrender1 \
        csh && \
    apt-get clean

RUN curl -L -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /opt/conda/bin/conda init

RUN /opt/conda/bin/conda create -n ligpargen -y python=3.7 && \
    /opt/conda/bin/conda install -n ligpargen -y -c rdkit rdkit && \
    /opt/conda/bin/conda install -n ligpargen -y -c conda-forge openbabel

ENV PATH="/opt/conda/envs/ligpargen/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

RUN git clone https://github.com/Isra3l/ligpargen.git /opt/ligpargen && \
    cd /opt/ligpargen && \
    /opt/conda/envs/ligpargen/bin/pip install -e .

COPY ./boss /opt/BOSSdir

RUN chmod +x /opt/BOSSdir/*

ENV BOSSdir="/opt/BOSSdir"

WORKDIR /opt/output

RUN echo "source activate ligpargen" > ~/.bashrc

SHELL ["/bin/bash", "-c"]

CMD ["/bin/bash"]

