FROM r-base:latest
USER root
RUN apt-get update && apt-get -y upgrade && apt-get -y install libcurl4-openssl-dev

WORKDIR /tmp
RUN R -e "install.packages('plotrix', dependencies = TRUE)"
RUN cd .. && rm -R tmp

RUN mkdir /data
RUN mkdir /home/SAFETRANS
ADD Visibility.R /home/SAFETRANS/Visibility.R
ADD scripta.R /home/SAFETRANS/scripta.R
ADD AERONET_data.txt /home/SAFETRANS/AERONET_data.txt.R
WORKDIR /home/SAFETRANS

CMD Rscript Visibility.R 1 FALSE azimuth null 355 FALSE null null TRUE TRUE
