FROM ubuntu:24.04

ENV DEBIAN_FORNTED=noninteractive

SHELL ["/bin/bash", "-c"]

RUN mkdir -p /opt/dna_methylation_smk && \
    chmod 777 -R /opt
WORKDIR /opt/dna_methylation_smk

RUN apt-get update -y && \
    apt-get install -y build-essential \
        aria2 git curl perl \
        zlib1g-dev libbz2-dev libz-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#! create a normal user and install other packages
#! as the normal user to prevent permission issues
RUN printf 'CREATE_MAIL_SPOOL=no' >> /etc/default/useradd \
    && mkdir -p /home/snake \
    && groupadd snake \
    && useradd snake -g snake -d /home/snake \
    && chown snake:snake /home/snake \
    && usermod -aG sudo snake

USER snake

RUN aria2c -c \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -d /tmp -o miniforge.sh && \
    bash /tmp/miniforge.sh -b -f -p /opt/miniforge && \
    rm /tmp/miniforge.sh

ENV PATH=/opt/miniforge/bin:$PATH

RUN mamba init && \
    source ~/.bashrc

RUN mamba create -n snakemake \
        -c conda-forge -c bioconda snakemake htslib -y && \
    mamba create -n cpython \
        -c conda-forge python=3.12 pybind11 -y

SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

COPY ./utils/download.sh utils/download.sh
RUN bash utils/download.sh

COPY ./fq2bedgraph.smk fq2bedgraph.smk
COPY ./config config
COPY ./rules rules
COPY ./utils/sample_sheet_parser.py utils/sample_sheet_parser.py

RUN mkdir -p input && \
    touch input/BS1_HF_1.R1.fq.gz input/BS1_HF_1.R2.fq.gz \
        input/EM1_HF_1.R1.fq.gz input/EM1_HF_1.R2.fq.gz \
        input/PS1_HF_1.R1.fq.gz input/PS1_HF_1.R2.fq.gz

COPY ./test.csv test.csv
RUN snakemake --snakefile fq2bedgraph.smk \
        --config sample_sheet=test.csv \
        --cores 1 --use-conda --conda-create-envs-only
