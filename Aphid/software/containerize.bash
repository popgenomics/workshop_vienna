#!/bin/bash

# Use this script to generate consistent singularity / docker definition files.

NAME="Aphid" # Project name.
REPO="aphid" # Source folder name.
WORKDIR="aphid" # Workding directory name.
BIN="aphid" # Executable name.

# Collect in this temporary variable.
DEPS=""
deps() {
    DEPS="${DEPS} $@"
}

# "Build" dependencies are removed after compilation.
BUILD_DEPS=""
bdeps() {
    BUILD_DEPS="${BUILD_DEPS} $@"
}

#==== Dependencies =============================================================
# Build-time only dependencies.
bdeps git               # To get source code.
bdeps gcc               # Compiler.
bdeps pacman-contrib    # For paccache.

# Runtime dependencies.
# deps <NONE>

#=== Construction layers =======================================================
read -r -d '' DEPENDENCIES_LAYER <<-EOF
    #---- Dependencies ---------------------------------------------------------
    echo 'Server = https://mirrors.kernel.org/archlinux/\$repo/os/\$arch' \\
         > /etc/pacman.d/mirrorlist
    pacman -Syu --noconfirm
    pacman -Sy --noconfirm ${BUILD_DEPS} ${DEPS}
    #---------------------------------------------------------------------------
EOF

read -r -d '' COMPILATION_LAYER <<-EOF
    #---- Compilation ----------------------------------------------------------
    # Get source code.
    cd /opt
    git clone --recursive https://gitlab.mbb.univ-montp2.fr/ibonnici/${REPO}
    cd ${REPO}

    # Compile
    mkdir bin
    gcc -lm *.c -o bin/${BIN}
    #---------------------------------------------------------------------------
EOF

read -r -d '' INSTALLATION_LAYER <<-EOF
    #---- Installation ---------------------------------------------------------
    ln -s /opt/${REPO}/bin/${BIN} /usr/bin/${BIN}
    #---------------------------------------------------------------------------
EOF

read -r -d '' CLEANUP_LAYER <<-EOF
    #---- Cleanup --------------------------------------------------------------
    # Remove the packages downloaded to image's Pacman cache dir.
    paccache -r -k0

    # Uninstall build dependencies.
    pacman -Rns --noconfirm ${BUILD_DEPS}
    #---------------------------------------------------------------------------
EOF

#=== Generate Singularity file =================================================

FILENAME="singularity.def"
cat <<EOF > $FILENAME
BootStrap: docker
From: archlinux

# This file was automatically generated from ./containerize.bash.

%post
    echo "Building container.."

    ${DEPENDENCIES_LAYER}

    ${COMPILATION_LAYER}

    ${INSTALLATION_LAYER}

    ${CLEANUP_LAYER}

    echo "export CONTAINER_BUILD_TIME=\\"\$(date)\\"" \\
        >> \${SINGULARITY_ENVIRONMENT}

%runscript
    echo "Running ${NAME} container (created on \${CONTAINER_BUILD_TIME})"
    ${BIN} \$@
EOF
echo "Generated $FILENAME"

#=== Generate Docker file ======================================================

FILENAME="Dockerfile"
cat <<DEOF > $FILENAME
# syntax=docker/dockerfile:1.3-labs
FROM archlinux

# This file was automatically generated from ./containerize.bash.

RUN <<EOF
    ${DEPENDENCIES_LAYER}
EOF

RUN <<EOF
    ${COMPILATION_LAYER}
EOF

RUN <<EOF
    ${INSTALLATION_LAYER}
EOF

RUN <<EOF
    ${CLEANUP_LAYER}
EOF

# Now pick a folder to work within.
RUN mkdir -p /home/${WORKDIR}
WORKDIR /home/${WORKDIR}

ENTRYPOINT ["${BIN}"]
DEOF
echo "Generated $FILENAME"

exit 0
