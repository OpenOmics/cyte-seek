# Dockerfile for chicyte cite pipeline
# Using Ubuntu Focal (20.04 LTS) as base image
FROM ubuntu:20.04
LABEL maintainer=kuhnsa@nih.gov

# Overview of Dependencies
#  • python/3.8    (apt-get python3)
#    - pandas      (pip install pandas)
#    - xlsxwriter  (pip install XlsxWriter)
#    - argparse    (standard libary >= 3.5)
#    - glob        (standard libary >= 3.5)
#    - csv         (standard libary >= 3.5)
#    - ntpath      (standard libary >= 3.5)
#    - shutil      (standard libary >= 3.5)
#  • R/4.1
#    - pandoc        (agt-get pandoc)
#    - argparse      (CRAN)
#    - Seurat/4.0.1  (CRAN)
#    - ggplot2       (CRAN)
#    - dplyr         (CRAN)
#    - gridExtra	 (CRAN)
#    - RColorBrewer	 (CRAN)
#    - future	     (CRAN)
#    - future.apply	 (CRAN)
#    - dsb           (CRAN)
#    - htmlwidgets   (CRAN)
#    - knitr         (CRAN)
#    - stringr       (CRAN)
#    - flextable     (CRAN)
#    - scales        (CRAN)
#    - rmarkdown     (CRAN)
#    - remotes       (CRAN)
#    - DT            (CRAN)
#    - plotly        (CRAN)
#    - devtools      (CRAN)
#    - BiocManager   (CRAN)
#    - MAST          (Bioconductor)
#    - Azimuth/0.4.5 (github satijalab/azimuth)
#  • demuxlet        (github)
#    - htslib/1.10   (github)
#  • vcftools/0.1.16 (apt-get)
#  • bcftools/1.10   (apt-get)
#  • samtools/1.10   (apt-get)

# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with host filesyetem 
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apt-transport-https \
        autotools-dev \
        automake \
        autoconf \
        bcftools \
        build-essential \
        bzip2 \
        ca-certificates \
        cmake \
        curl \
        dirmngr \
        g++ \
        gcc \
        gfortran \
        git \
        gnupg \
        gzip \
        libfontconfig1-dev \
        libgeos-dev \
        libhdf5-dev \
        libtool \
        libz-dev \
        locales \
        make \
        pandoc \
        python3 \
        python3-pip \
        samtools \
        software-properties-common \
        unzip \
        vcftools \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install python packages
RUN ln -s /usr/bin/python3.8 /usr/bin/python
RUN pip3 install --upgrade pip \
    && pip3 install argparse \
    && pip3 install pandas \
    && pip3 install XlsxWriter

# Install R (4.X) and required packages
# For more information, check out: 
# https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev
# Install required R packages
# Seurat (needs libgeos-dev), flextable (needs libfontconfig1-dev), Azimuth (libhdf5-dev)
# CRAN packages
RUN Rscript -e 'install.packages(c("DT", "RColorBrewer", "Seurat"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("argparse", "devtools", "dplyr"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("flextable", "future", "future.apply"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("ggplot2", "gridExtra", "htmlwidgets"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("knitr", "plotly", "rmarkdown"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("scales", "stringr", "BiocManager"), repos="http://cran.r-project.org")'
# Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("MAST", "limma"))'
# dsb needs limma installed
RUN Rscript -e 'install.packages(c("dsb"), repos="http://cran.r-project.org")'
# Github packages
RUN Rscript -e 'update.packages(oldPkgs = c("withr", "rlang"))'
RUN Rscript -e 'remotes::install_github("satijalab/azimuth", ref = "master")'

# Install and build HTSlib 1.10 (required by demuxlet)
# Needs to be a sibling directory to demuxlet install,
# the name of the htslib directory must be called htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 \
    && tar -vxjf htslib-1.10.tar.bz2 \
    && rm htslib-1.10.tar.bz2 \
    && cd htslib-1.10 \
    && autoreconf -i \
    && ./configure \
    && make \
    && make install \
    && ln -s /opt2/htslib-1.10 /opt2/htslib
WORKDIR /opt2 

# Install and build demuxlet, 
# configure command needs libtool installed first
RUN git clone https://github.com/statgen/demuxlet.git \
    && cd demuxlet \
    && autoreconf -vfi \
    && ./configure \
    && make \
    && make install
WORKDIR /opt2

# Add Dockerfile and export env variables
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2