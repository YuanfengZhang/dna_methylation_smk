# To run this Dockerfile, use the following command:
#   docker built \
#        -t ${pkg}:{ver} \
#        -f Rpkgs.Dockerfile \
#        --build-arg CORES=8 \
#        --build-arg PACKAGE_NAME=${pkg} \
#        .
FROM ubuntu:24.04

ENV DEBIAN_FORNTED=noninteractive

SHELL ["/bin/bash", "-c"]

RUN apt-get clean -y && rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential gfortran libblas-dev \
    liblapack-dev libcurl4-openssl-dev \
    libssl-dev libxml2-dev libcairo2-dev \
    libpango1.0-dev libreadline-dev \
    libtcl8.6 libtk8.6 tcl-dev \
    tk-dev xorg-dev zlib1g-dev \
    libfontconfig1-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev \
    libpng-dev libtiff5-dev libjpeg-dev \
    libbz2-dev liblzma-dev libdeflate-dev \
    autoconf pkg-config \
    software-properties-common dirmngr \
    curl gnupg pigz zstd

ARG CORES=4

RUN mkdir -p /opt/htslib
RUN curl \
    -L -o /opt/htslib/htslib-1.21.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
RUN cd /opt/htslib/ && \
    tar -xvf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    make -j "${CORES}" && \
    make install

# ! If R=4.4.3 is needed, uncomment line 45 to 52 and comment line 54 to 62.
# RUN gpg --keyserver keyserver.ubuntu.com --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9
# RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 |\
#     apt-key add -

# RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/"

# RUN apt-get update && \
#     apt-get install -y --no-install-recommends r-base=4.4.3* r-base-dev r-base-html

RUN curl \
    -L -o /opt/R-alpha.tar.xz \
    https://cran.r-project.org/src/base-prerelease/R-alpha_2025-03-18_r87998.tar.xz
RUN cd /opt && \
    tar -xvf R-alpha.tar.xz && \
    cd R-alpha && \
    ./configure && \
    make -j "${CORES}" && \
    make install

COPY install_rpkg.R /opt/install_rpkg.R

ARG PACKAGE_NAME

RUN if [ -n "$PACKAGE_NAME" ]; then \
    echo "Installing ${PACKAGE_NAME}"; \
    Rscript /opt/install_rpkg.R --cores ${CORES} --package ${PACKAGE_NAME}; \
    else \
    echo "No package to install"; \
    fi

RUN apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

RUN rm -rf \
    /opt/htslib \
    /opt/htslib/htslib-1.21.tar.bz2 \
    R-alpha \
    /opt/R-alpha.tar.xz
