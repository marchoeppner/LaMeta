FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for LaMeta pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/LaMeta-1.0/bin:$PATH
