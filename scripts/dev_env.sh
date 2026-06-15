#!/usr/bin/env bash
# Source this file before local GEOSldas development commands:
#   source scripts/dev_env.sh

if [[ "${BASH_SOURCE[0]:-}" == "${0}" ]]; then
  echo "Source this file instead of executing it: source scripts/dev_env.sh" >&2
fi

# Keep system/Homebrew toolchain ahead of conda to avoid old linker/tool shadowing.
SAFE_PATH_PREFIX="/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:/opt/homebrew/sbin"
case "${PATH}" in
  "${SAFE_PATH_PREFIX}:"*) ;;
  *)
    export PATH="${SAFE_PATH_PREFIX}:${PATH}"
    ;;
esac

export FC="${FC:-gfortran}"
export LDFLAGS="${LDFLAGS:--L/opt/homebrew/opt/openmpi/lib}"
export CPPFLAGS="${CPPFLAGS:--I/opt/homebrew/opt/openmpi/include}"

# Some local MPI setups do not populate this CMake variable; keep a safe default.
export GEOSLDAS_MPI_VERSION_STRING="${GEOSLDAS_MPI_VERSION_STRING:-Unknown}"
