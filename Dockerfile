# Build this with:
# docker build -t registry.gitlab.com/smspp/milpsolver .
# Run with:
# docker run --rm -it registry.gitlab.com/smspp/milpsolver:latest
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
