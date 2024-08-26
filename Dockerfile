# Base OS
FROM python:3.10-slim-bookworm

ARG VERSION="1.5.0"
ARG SIMULATOR_VERSION="2.73"

# metadata
LABEL \
    org.opencontainers.image.title="Smoldyn" \
    org.opencontainers.image.version="${SIMULATOR_VERSION}" \
    org.opencontainers.image.description="BioSimulators-compliant command-line interface to the Smoldyn simulation program" \
    org.opencontainers.image.url="https://github.com/ssandrews/Smoldyn" \
    org.opencontainers.image.documentation="https://Smoldyn.readthedocs.io/" \
    org.opencontainers.image.source="https://github.com/biosimulators/Biosimulators_Smoldyn" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="BSD-3-Clause" \
    \
    base_image="python:3.9-slim-buster" \
    version="${VERSION}" \
    software="Smoldyn" \
    software.version="${SIMULATOR_VERSION}" \
    about.summary="BioSimulators-compliant command-line interface to the Smoldyn simulation program" \
    about.home="https://github.com/ssandrews/Smoldyn" \
    about.documentation="https://Smoldyn.readthedocs.io/" \
    about.license_file="https://github.com/ssandrews/Smoldyn/blob/master/LICENSE.md" \
    about.license="SPDX:BSD-3-Clause" \
    about.tags="BioSimulators,mathematical model,kinetic model,simulation,systems biology,computational biology,stochastic,spatial,particle-simulation,SED-ML,COMBINE,OMEX" \
    extra.identifiers.biotools="Smoldyn" \
    maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        g++ \
        libatlas-base-dev \
        swig \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# fonts for matplotlib
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends libfreetype6 \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_Smoldyn
RUN pip install pip==23.0.1
RUN pip install sympy /root/Biosimulators_Smoldyn \
    && rm -rf /root/Biosimulators_Smoldyn
#RUN pip install sympy /root/Biosimulators_Smoldyn Smoldyn==${SIMULATOR_VERSION} \
#    && rm -rf /root/Biosimulators_Smoldyn
ENV VERBOSE=0 \
    MPLBACKEND=PDF

# Entrypoint
ENTRYPOINT ["biosimulators-smoldyn"]
CMD []