BootStrap: docker
From: ubuntu:18.04

%post
    R_BASE_VERSION="3.5.1"

    # Setup system packages
    apt-get -qq update
    apt-get upgrade -y
    apt-get install -y \
        apt-transport-https \
        gnupg \
        ca-certificates \
        libc6 \
        libcurl4-openssl-dev \
        libxml2-dev \
        libfftw3-dev \
        git \
        wget \
        zip \
        libssl-dev \
        vim-tiny \
        libglu1-mesa-dev \
        locales
 
    locale-gen en_US
    update-locale LANG=en_US

    # Setup R repository
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
    apt-get -qq update
    apt-get upgrade -y

    # Install R
    apt-get install -y --no-install-recommends \
    littler \
    r-base-core=${R_BASE_DEV}* \
    r-base-dev=${R_BASE_DEV}*
 
    # Prepare R repositories
    echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site
    echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r

# Install NormalyzerDE dependencies

    Rscript -e 'install.packages("BiocManager", dependencies=TRUE)'

    # CRAN
    Rscript -e 'install.packages("devtools", dependencies=TRUE)'
    Rscript -e 'install.packages("Rcmdr", dependencies=TRUE)'
    Rscript -e 'install.packages("ape", dependencies=TRUE)'
    Rscript -e 'install.packages("raster", dependencies=TRUE)'
    Rscript -e 'install.packages("RcmdrMisc", dependencies=TRUE)'
    Rscript -e 'install.packages("gridExtra", dependencies=TRUE)'
    Rscript -e 'install.packages("ggplot2", dependencies=TRUE)'
    # Rscript -e 'BiocManager::install("grid", dependencies=TRUE)'
    Rscript -e 'install.packages("hexbin", dependencies=TRUE)'
    Rscript -e 'install.packages("htmlTable", dependencies=TRUE)'
    Rscript -e 'install.packages("stringi", dependencies=TRUE)'

    # Bioconductor
    Rscript -e 'BiocManager::install("vsn", dependencies=TRUE)'
    Rscript -e 'BiocManager::install("preprocessCore", dependencies=TRUE)'
    Rscript -e 'BiocManager::install("limma", dependencies=TRUE)'
    Rscript -e 'BiocManager::install("MASS", dependencies=TRUE)'
    Rscript -e 'BiocManager::install("car", dependencies=TRUE)'
    Rscript -e 'BiocManager::install("SummarizedExperiment", dependencies=TRUE)'

# NormalyzerDE setup

    # NORMALYZERDE_VERSION="0.99.5"
    # wget https://github.com/ComputationalProteomics/NormalyzerDE/releases/download/${NORMALYZERDE_VERSION}/NormalyzerDE_${NORMALYZERDE_VERSION}.tar.gz
    
    Rscript -e 'devtools::install_github("ComputationalProteomics/NormalyzerDE")'

    # R CMD INSTALL NormalyzerDE_${NORMALYZERDE_VERSION}.tar.gz

    

