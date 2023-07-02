#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

HELP(){
	echo "

	Installing align_annotate pipeline.

	This pipeline requires two coding environments to run properly:

	1) conda: a package manager that we will use to install all dependencies.
   	   to install conda, use the instructions found at 

	   https://conda.io/projects/conda/en/latest/user-guide/install/index.html

	2) docker: a sandbox for running software with many dependencies.
	   we need to use docker because some of the dependencies for running
	   smoove, the program we use to call indels, do not have conda packages 
	   available for OSX or Windows.

		to install docker, follow the instructions found at:
		
		https://docs.docker.com/get-docker/"
}

#check if docker and conda are installed
if ! command -v docker &> /dev/null
then 
		echo "docker could not be found"
		HELP
else
		echo "found docker"
fi


if ! command -v conda &> /dev/null
then 
		echo "conda could not be found"
		HELP
else
		echo "found conda"
fi

##install conda dependencies
echo "CREATING CONDA ENVIRONMENT"
conda env create -f align_annotate.yml
conda activate align_annotate

##download smoove docker image
echo "DOWNLOADING SMOOVE" 
docker pull brentp/smoove
#docker run -it brentp/smoove smoove -h


