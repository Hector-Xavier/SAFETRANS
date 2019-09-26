FROM r-base:latest
USER root
RUN apt-get update && apt-get -y upgrade && apt-get -y install libcurl4-openssl-dev

WORKDIR /tmp
RUN R -e "install.packages('plotrix', dependencies = TRUE)"
RUN cd .. && rm -R tmp

RUN mkdir /data
RUN mkdir /home/SAFETRANS
ADD SAFETRANS_visibility_module.R /home/SAFETRANS/SAFETRANS_visibility_module.R
ADD scripta.R /home/SAFETRANS/scripta.R
WORKDIR /home/SAFETRANS

CMD Rscript SAFETRANS_visibility_module.R 1 FALSE azimuth null FALSE null null TRUE TRUE
