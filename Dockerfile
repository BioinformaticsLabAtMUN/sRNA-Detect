# syntax=docker/dockerfile:1
# Download ubuntu as base image
FROM ubuntu

# install Python3 and its libraries
RUN apt-get update --fix-missing -qq && apt-get -y install -y -q \
    software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3 && \
    apt-get install -y python3-pip && \
    apt-get install -y build-essential python-numpy python-pysam python-htseq

# install bedtools
RUN apt-get update && apt-get install -y bedtools
