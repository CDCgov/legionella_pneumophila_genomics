#!/bin/bash

# get user and group IDs based on bind mounted output folder
DOCKER_UID=$(stat -c %u /output)
DOCKER_GID=$(stat -c %g /output)

# add user and group to container
addgroup -g $DOCKER_GID dockergroup
adduser -u $DOCKER_UID -g $DOCKER_GID -D -s /bin/sh dockeruser dockergroup

# make sure the pipeline directory is owned by the user/group
chown -R $DOCKER_UID:$DOCKER_GID /opt/pipeline

# run the pipeline
su-exec $DOCKER_UID:$DOCKER_GID /opt/pipeline/pipeline.sh --reference=/opt/pipeline/reference/Phila_NC_002942.fna --gff=/opt/pipeline/reference/NC_002942.gff --output=/output "$@"

