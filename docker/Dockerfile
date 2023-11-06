FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y build-essential

RUN mkdir /tools
WORKDIR /tools
COPY environment.yaml . 
RUN . /opt/conda/bin/activate && \
    conda env create -n strainline -f environment.yaml && \
    conda clean -a && \
    wget https://github.com/gt1/daccord/releases/download/0.0.10-release-20170526170720/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz && \
    tar -zvxf daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz && \
    ln -fs /tools/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu/bin/daccord /opt/conda/envs/strainline/bin/daccord && \
    git clone https://github.com/HaploKit/Strainline.git && \
    rm daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz environment.yaml
ENV PATH=/opt/conda/envs/strainline/bin/:$PATH
ENV PATH=/tools/Strainline/src/:$PATH
