# --------------------------------------------------------------------------- #
#    Dockerfile for CI/CD                                                     #
#                                                                             #
#    This file contains the commands to build a Docker image containing       #
#    all the packages needed to build and test the project.                   #
#    Once built and uploaded in the repository's container registry           #
#    (See: https://gitlab.com/smspp/milpsolver/container_registry),           #
#    the image can be fetched and used by the GitLab Runner.                  #
#                                                                             #
#    Note: Once built, this image will contain parts of the                   #
#          IBM ILOG CPLEX Optimization Studio suite. For this reason,         #
#          the image cannot be publicly distributed and can be used only      #
#          by team members that own a valid license for that software.        #
#    ---------------------------------------------------------------------    #
#    DISCLAIMER: The author of this file is not affiliated, associated,       #
#    authorized, endorsed by, or in any way officially connected with IBM,    #
#    or any of its subsidiaries or its affiliates. The names IBM, ILOG and    #
#    CPLEX as well as related names, marks, emblems and images are            #
#    registered trademarks of their respective owners.                        #
#    ---------------------------------------------------------------------    #
#                                                                             #
#    Copy the ILOG directory from a Linux package in the build directory.     #
#    Note that the directory must be in the Docker context.                   #
#    (See: https://docs.docker.com/engine/reference/builder/#copy).           #
#                                                                             #
#    Build the image with:                                                    #
#                                                                             #
#        $ docker build --build-arg ILOG_DIR=<ILOG-DIR> \                     #
#                       -t registry.gitlab.com/smspp/milpsolver .             #
#                                                                             #
#    Upload with:                                                             #
#                                                                             #
#        $ docker push registry.gitlab.com/smspp/milpsolver                   #
#                                                                             #
#    Run (locally) with:                                                      #
#                                                                             #
#        $ docker run --rm -it registry.gitlab.com/smspp/milpsolver:latest    #
#                                                                             #
#    Note: you need to rebuild and upload the image only when this file       #
#          changes, not when MILPSolver changes.                              #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# An image containing everything needed to build SMS++
FROM registry.gitlab.com/smspp/smspp

# Install required packages
RUN set -ex; \
		apt-get update; \
		apt-get install -y --no-install-recommends \
		liblapack-dev libgfortran-7-dev libtbb-dev; \
		rm -rf /var/lib/apt/lists/*;

# Install SCIP
RUN set -ex; \
	mkdir -p "scip" && cd "scip"; \
	curl -SL -o scip.sh "https://www.scipopt.org/download/release/SCIPOptSuite-7.0.1-Linux.sh"; \
	chmod u+x scip.sh; \
	./scip.sh --prefix=/usr/local --exclude-subdir --skip-license; \
	cd .. && rm -r "scip";

# Install CPLEX
ARG ILOG_DIR=ILOG
COPY $ILOG_DIR /opt/ibm/ILOG
