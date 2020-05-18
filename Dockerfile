# Dockerfile to build the OSS DBS environment
#
# Specifically the following stages are contained in that order extending the
# previous stage:
# 1. python_dep
# 2. app
#
# Author:
# Max Schroeder <max.schroeder@uni-rostock.de>


###############################################################################
# Stage 1
FROM sfbelaine/oss_dbs:fenics19_latest as python_dep

ARG OSS_UID="1000"
ARG OSS_GID="1000"

ENV OSS_UID=$OSS_UID \
    OSS_GID=$OSS_GID

WORKDIR /opt/oss_platform

# Install python packages
COPY requirements.txt /opt/oss_platform/requirements.txt
RUN pip3 install -r requirements.txt

# Create a new user and group in order to run processes as non-root
RUN groupadd -g $OSS_GID oss_platform
RUN useradd -m -s /bin/bash -g $OSS_GID -u $OSS_UID oss_platform

# switch user
USER $OSS_UID

###############################################################################
# Stage 2
FROM python_dep as app

# Copy all from the current source code repository
COPY --chown=oss_platform:oss_platform . /opt/oss_platform

WORKDIR /opt/oss_platform/OSS_platform

CMD [ "/bin/bash" ]
