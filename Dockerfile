FROM mambaorg/micromamba:2.0.5

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.12 plink=1.90b7.7 plink2=2.0.0a.6.9 polars=1.32 pyyaml=6.0 tomli=2.2 awscli=1.40 \
    && micromamba clean --all --yes

WORKDIR /opt/pgs_pipeline
COPY config.py combine_scores.py PLINK_PIPELINE_*.py ./
COPY sumstats/ ./sumstats/

ENTRYPOINT []
CMD ["python", "--version"]
