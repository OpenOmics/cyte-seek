# Dockerfile for ACTINN
# https://github.com/mafeiyang/ACTINN
# Using Ubuntu Focal (20.04 LTS) as base image
FROM ubuntu:20.04
LABEL maintainer=kuhnsa@nih.gov

# Overview of Dependencies
#  • python>=3.6         (make python3.6)
#    - tensorflow>=1.10  (pip install tensorflow)
#    - numpy>=1.14       (pip install numpy)
#    - pandas>=0.23      (pip install pandas)
#    - scipy>=1.1        (pip install scipy)
#    - argparse          (standard libary >= 3.5)

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
        build-essential \
        ca-certificates \
        curl \
        git \
        gzip \
        libreadline-gplv2-dev \
        libncursesw5-dev \
        libssl-dev \
        libsqlite3-dev \
        tk-dev \
        libgdbm-dev \
        libc6-dev \
        libbz2-dev \
        liblzma-dev \
        lzma \
        locales \
        make \
        unzip \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3.6 and set as default interpreter
# and install python packages
RUN wget https://www.python.org/ftp/python/3.6.10/Python-3.6.10.tgz \
    && tar -xvzf Python-3.6.10.tgz \
    && rm Python-3.6.10.tgz \
    && cd Python-3.6.10 \
    && ./configure --enable-optimizations \
    && make altinstall
RUN ln -sf /usr/local/bin/python3.6 /usr/bin/python \
    && ln -sf /usr/local/bin/python3.6 /usr/bin/python3 \
    && ln -sf /usr/local/bin/pip3.6 /usr/local/bin/pip3
RUN pip3 install --upgrade pip \
    && pip3 install tensorflow==1.10.1 \
    && pip3 install numpy==1.14 \
    && pip3 install pandas==0.23 \
    && pip3 install scipy==1.1 \
    && pip3 install tables

# Install ACTINN from Github
# https://github.com/mafeiyang/ACTINN
# from commit/a32686697aa9aaadbd96159523b0dfad04e93edd
WORKDIR /opt2
RUN git clone https://github.com/mafeiyang/ACTINN.git \
    && sed -i '1 i\#!/usr/bin/env python3' /opt2/ACTINN/*.py \
    && chmod +x /opt2/ACTINN/*.py

# Add Dockerfile and export env variables
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2/ACTINN:$PATH"
WORKDIR /data2
