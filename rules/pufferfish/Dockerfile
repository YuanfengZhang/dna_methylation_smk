FROM ubuntu:24.04
# 安装基础依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates build-essential \
    aria2 git curl perl-base wget ssh gcc g++ cmake make \
    zlib1g-dev liblzma-dev libbz2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*\
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN aria2c -c https://www.github.com/Kitware/CMake/releases/download/v3.22.1/cmake-3.22.1-linux-x86_64.sh
RUN chmod +x cmake-3.22.1-linux-x86_64.sh && mkdir -p /opt/cmake && \
    ./cmake-3.22.1-linux-x86_64.sh --skip-license --prefix=/opt/cmake
RUN aria2c -c https://ftp.gnu.org/gnu/make/make-4.3.tar.gz
RUN tar -xvf make-4.3.tar.gz && \
    cd make-4.3 && \
    ./configure --prefix=/opt/make && \
    make -j 4 && \
    make install
ENV PATH="/opt/cmake/bin:/opt/make/bin:$PATH"

RUN git clone https://github.com/COMBINE-lab/pufferfish -b develop 

RUN mkdir -p pufferfish/build && cd pufferfish/build && \
    { /opt/cmake/bin/cmake ../  || cmake ../; } && \
    { /opt/make/bin/make VERBOSE=1 -j 24 || \
      make VERBOSE=1 -j 24 || \
      /opt/make/bin/make VERBOSE=1 -j 24 || \
      make VERBOSE=1 -j 24 || \
      /opt/make/bin/make VERBOSE=1 -j 24 || \
      make VERBOSE=1 -j 24; }

RUN rm cmake-3.22.1-linux-x86_64.sh

SHELL ["/bin/bash", "--login", "-c"]
CMD ["bash"]
