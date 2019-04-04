FROM ubuntu:14.04

RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y build-essential
RUN apt-get install -y gcc
RUN apt-get install -y wget
RUN apt-get install -y cmake
RUN apt-get install -y libz-dev

RUN wget https://github.com/fritzsedlazeck/Sniffles/archive/master.tar.gz -O Sniffles.tar.gz
RUN tar xzvf Sniffles.tar.gz
WORKDIR Sniffles-master/
RUN mkdir -p build/
WORKDIR build
RUN cmake ..
RUN make

