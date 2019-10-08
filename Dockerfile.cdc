# build singularity

FROM golang:1.12.5-alpine3.9 as singularity_builder 

RUN apk update \
    && apk add --no-cache \
        ca-certificates \
        build-base \
        linux-headers \
        git \
        openssl-dev \
        util-linux-dev \
        libseccomp-dev \
        gawk \
    && mkdir -p $GOPATH/src/github.com/sylabs \
    && cd $GOPATH/src/github.com/sylabs \
    && git clone -b v3.3.0 https://github.com/sylabs/singularity.git \
    && cd $GOPATH/src/github.com/sylabs/singularity \
    && ./mconfig -s -S -p /opt/singularity \
    && cd $GOPATH/src/github.com/sylabs/singularity/builddir \
    && make \
    && make install


# build final container

FROM alpine:3.9

COPY --from=singularity_builder /opt/singularity /opt/singularity

COPY pipeline.sh /opt/pipeline/
COPY supportFiles/Phila_NC_002942.fna /opt/pipeline/reference/
COPY supportFiles/NC_002942.gff /opt/pipeline/reference/
COPY entrypoint.sh /

# custom certs required for containers to run in CDC network
COPY certs/* /etc/ssl/certs/

RUN apk update \
    && apk add --no-cache \
        ca-certificates \
        squashfs-tools \
        libseccomp \
        su-exec \
        bash \
    && chmod u+s /opt/singularity/libexec/singularity/bin/starter-suid \
    && chmod -R 755 /opt/pipeline \
    && chmod 644 /etc/ssl/certs/*

ENV PATH /opt/singularity/bin:$PATH

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]

