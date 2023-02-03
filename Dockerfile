FROM ubuntu:22.04
RUN apt-get update && \
apt-get upgrade --assume-yes && \
DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC \
apt-get install --assume-yes build-essential git cmake libopenblas-dev liblapacke-dev && mkdir root_digger
COPY CMakeLists.txt makefile root_digger/
COPY lib/ root_digger/lib/
COPY src/ root_digger/src/
COPY .git/ root_digger/.git/
RUN cd root_digger && make && cp bin/rd /usr/local/bin/
CMD ["rd","--help"]
