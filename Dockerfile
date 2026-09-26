# Build image for the obdeect runtime consumed by simtools.
#
# The image deliberately contains the current reference executables.  The
# production `obdeect-simtools-raytrace` executable referenced by
# simtools will be added once the production scene/arrival contract is complete.
#
# Build from the obdeect repository root:
#   docker build -t ghcr.io/gammasim/obdeect:v0.1.0 .

ARG BUILD_BASE_IMAGE=docker.io/library/almalinux:9.8
ARG RUNTIME_IMAGE=docker.io/library/almalinux:9.8-minimal
ARG OBDEECT_VERSION=0.1.0
ARG OBDEECT_REVISION=unknown

FROM ${BUILD_BASE_IMAGE} AS build

ARG OBDEECT_VERSION

RUN dnf install -y cmake gcc-c++ make && \
    dnf clean all && \
    rm -rf /var/cache/dnf/* /tmp/* /var/tmp/*

WORKDIR /src/obdeect
COPY CMakeLists.txt ./
COPY cmake/ cmake/
COPY cpp/ cpp/

RUN cmake -S . -B build -DOBDEECT_BUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build --parallel && \
    cmake --install build --prefix /opt/obdeect

RUN printf '%s\n' \
      "obdeect_version: \"${OBDEECT_VERSION}\"" \
      'obdeect_build_type: "Release"' \
      'obdeect_build_testing: false' \
      'reference_executables:' \
      '  - obdeect_toy' \
      '  - obdeect_ctao' > /opt/obdeect/build_opts.yml

FROM ${RUNTIME_IMAGE}

ARG OBDEECT_VERSION
ARG OBDEECT_REVISION

RUN microdnf install -y libstdc++ && \
    microdnf clean all && \
    rm -rf /var/cache/dnf/* /tmp/* /var/tmp/*

COPY --from=build /opt/obdeect /opt/obdeect

ENV PATH="/opt/obdeect/bin:${PATH}"
WORKDIR /opt/obdeect

LABEL org.opencontainers.image.source="https://github.com/gammasim/obdeect" \
      org.opencontainers.image.revision="${OBDEECT_REVISION}" \
      org.opencontainers.image.version="${OBDEECT_VERSION}"
