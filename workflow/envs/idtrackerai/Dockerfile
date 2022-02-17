FROM continuumio/miniconda3:4.10.3p1
ENV DEBIAN_FRONTEND=noninteractive
RUN conda install python=3.6 &&\
    /opt/conda/bin/pip3 install "idtrackerai[gui]"
RUN apt update -y &&\
    apt install -y locales &&\
    echo "LC_ALL=en_US.UTF-8" >> /etc/environment &&\
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen &&\
    echo "LANG=en_US.UTF-8" > /etc/locale.conf &&\
    locale-gen en_US.UTF-8 &&\
    apt install -y qtbase5-dev
# If run with singularity the guest user is the same as the host user
# If run from docker the user is root
# this is required for idtrackerai_test run on singularity since idtrackerai_test writes to its own directory
RUN chmod -R 777 /opt/conda/lib/python3.6/site-packages/idtrackerai*
ENTRYPOINT idtrackerai
LABEL note="\
This is not an official container for idtrackerai. \
This countainer can only be run with singularity (required for the GUI). \
It was built from continuumio/miniconda3:4.10.3p1 and the pip package idtrackerai. \
Refer to the respective maintainers for problems with those. \
Absolutely no warranty implied, use at your own risk.\
"
LABEL how.to.test="\
Open a shell in the container using singularity \'run --writable-tmpfs\' and run \
\'idtrackerai_test -o [output_dir]\'. Note that the output directory must be specified \
otherwise idtrackerai_test tries to write on its own directory in /opt and fails since \
there is not enough space. The \'--writable-tmpfs\' flag is also required since \
idtrackerai_test also edits some small files in its own directory in /opt even when the \
output directory is set.\
"
LABEL maintainer="saul@ebi.ac.uk"
