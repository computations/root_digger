FROM ubuntu:22.04
RUN apt-get update && \
apt-get upgrade --assume-yes && \
DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC \
apt-get install --assume-yes build-essential git cmake libopenblas-dev liblapacke-dev pkgconf && mkdir root_digger
RUN git clone https://github.com/computations/root_digger/ --depth=1 --recursive
RUN cd root_digger && make && cp bin/rd /usr/local/bin/
CMD ["--help"]
ENTRYPOINT ["rd"]
