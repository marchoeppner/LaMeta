FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for LaMeta pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/LaMeta-1.0/bin:$PATH

RUN mkdir -p /db/checkm_data && \
	cd /db/checkm_data && \
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
	tar xzf checkm_data_2015_01_16.tar.gz && \
	rm *.tar.gz && \
	printf "/db/checkm_data\n/db/checkm_data\n" | checkm data setRoot

