BootStrap: docker
From: bioconductor/bioconductor_docker

# Last build date: 202117
# Test build locally with singularity build test_out.simg Singularity

%post

    Rscript -e 'install.packages("devtools")'
    Rscript -e 'devtools::install_github("ComputationalProteomics/NormalyzerDE", dependencies=TRUE)'
    

