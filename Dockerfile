FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update

# -----------------------------------------

RUN apt-get install libffi-dev libssl-dev
RUN pip install --upgrade requests[security]

# Install Mummer
RUN \
    wget http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz && \
    tar xf MUMmer3.23.tar.gz && \
    mv MUMmer3.23 /kb/deployment/mummer && \
    cd /kb/deployment/mummer && make all && \
    echo 'export PATH=$PATH:/kb/deployment/mummer' >> /kb/deployment/user-env.sh

# RUN \
#     cp ../tools/dnadiff_genomes.pl /kb/deployment && \
#     mv /kb/deployment/dnadiff_genomes.pl /kb/deployment/dnadiff_genomes && \
#     chmod a+x /kb/deployment/dnadiff_genomes

# Copy local wrapper files, and build

COPY ./ /kb/module
RUN mkdir -p /kb/module/work

WORKDIR /kb/module

RUN make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
