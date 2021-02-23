#!/bin/bash
# My example bash script
#DIRECTORY=`dirname $0`
#echo $DIRECTORY

docker run --name OSS_container --volume $1:/opt/OSS-DBS -it --rm custom_oss_platform python3 Launcher_OSS_lite.py
