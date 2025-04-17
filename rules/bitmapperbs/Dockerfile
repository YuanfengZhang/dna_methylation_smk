FROM ubuntu:24.04
# 安装基础依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 安装micromamba
# RUN curl -L https://github.com/mamba-org/micromamba-releases/releases/download/2.1.0-0/micromamba-linux-64 \
#     -o /usr/local/bin/micromamba \
#     && chmod +x /usr/local/bin/micromamba

COPY micromamba-linux-64 /usr/local/bin/micromamba
RUN chmod 777 /usr/local/bin/micromamba

# 设置环境变量
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH="/opt/conda/bin:$PATH"

# 创建conda环境
RUN micromamba create -n bitmapperbs -y \
    -c bioconda -c conda-forge -c defaults \
    samtools>=1.21 \
    bitmapperbs=1.0.2.3=hc088bd4_2 \
    mbuffer

# 配置shell环境
RUN /usr/local/bin/micromamba \
     shell init --shell bash --root-prefix="$MAMBA_ROOT_PREFIX" \
    && echo "micromamba activate bitmapperbs" >> ~/.bashrc

# 清理缓存
RUN /usr/local/bin/micromamba clean -a -y

# 设置默认启动命令
SHELL ["/bin/bash", "--login", "-c"]
CMD ["bash"]