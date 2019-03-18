Bootstrap:docker
From:continuumio/anaconda

%labels
    MAINTAINER Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for the LaMeta pipeline
    VERSION 1.0

%environment
    PATH=/opt/conda/envs/LaMeta-1.0/bin:$PATH
    export PATH
    GTDBTK_DATA_PATH=/refdata/
    export GTDBTK_DATA_PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a

    mkdir -p /ifs
    mkdir -p /refdata
    apt-get -y install procps

