FROM ubuntu:16.04

# R setup

ENV R_BASE_VERSION 3.4.1

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get -qq update
RUN apt-get upgrade -y


RUN apt-get install -y --no-install-recommends \
    littler \
    r-base-core=${R_BASE_VERSION}* \
    r-base-dev=${R_BASE_DEV}* \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfftw3-dev \
    git \
    wget \
    zip
    

RUN echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site
RUN echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r

# Install Normalyzer dependencies

RUN Rscript -e 'install.packages("devtools", dependencies=TRUE)'

RUN Rscript -e 'install.packages("Rcmdr", dependencies=TRUE)'
RUN Rscript -e 'install.packages("ape", dependencies=TRUE)'
RUN Rscript -e 'install.packages("raster", dependencies=TRUE)'
RUN Rscript -e 'install.packages("gridExtra", dependencies=TRUE)'
RUN Rscript -e 'install.packages("ggplot2", dependencies=TRUE)'
RUN Rscript -e 'install.packages("grid", dependencies=TRUE)'
RUN Rscript -e 'install.packages("hexbin", dependencies=TRUE)'

## Bioconductor

RUN > rscript.R \
	&& echo 'source("https://bioconductor.org/biocLite.R")' >> rscript.R \
	&& echo 'biocLite(ask=FALSE)' >> rscript.R \
	#&&echo 'biocLite("BiocUpgrade")' >> rscript.R \
	&& echo 'biocLite("vsn", ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("preprocessCore", ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("limma", ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("MASS", ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("car", ask=FALSE)' >> rscript.R \
	&& Rscript rscript.R

# Normalyzer setup

WORKDIR /root
RUN wget -O $(pwd)/normalyzer_source.zip https://lu.box.com/shared/static/aiy9zrpj6aeej48b6elwhckk1a20ue7a.zip
RUN unzip $(pwd)/normalyzer_source.zip
RUN wget -O $(pwd)/normalyzer_test_data.zip https://lu.box.com/shared/static/72h0xqcvycyz1qck9f7c1913n6adwbz7.zip
RUN unzip $(pwd)/normalyzer_test_data.zip
RUN rm *.zip

