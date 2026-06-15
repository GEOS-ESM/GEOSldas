# GEOSldas Mac Build Notes

This note captures the known-good local build environment for this workspace on macOS.

## Component Checkout

Use `mepo` to bring the component tree into the versions requested by `components.yaml`.
Do this before configuring CMake:

```bash
mepo init
mepo fetch --tags
mepo restore-state --parallel
mepo develop GEOSgcm_GridComp GEOSldas_GridComp
```

The final `mepo develop` command is useful when actively working from the `develop`
branches of the grid components after restoring the tagged dependencies.

## Quick Startup Check

Before configuring/building, verify linker precedence:

```bash
which ld
```

Expected: `/usr/bin/ld`

If `which ld` points into Anaconda (for example `/Users/amfox/opt/anaconda3/bin/ld`), builds can fail with:

- `unsupported tapi file type '!tapi-tbd'`

## Known-Good Configure and Build Commands

Use a sanitized PATH. This avoids accidentally picking up linker/compiler tools from
Anaconda or other user environments:

```bash
SANEPATH="/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin"
```

Configure directly with CMake. Do not use `parallel_build.csh` on this local Mac
workspace; it expects environment setup such as `BASEDIR` that is not present here.

```bash
ROOT="/Users/amfox/Desktop/GEOSldas_develop"
SRC="$ROOT/GEOSldas"
PREFIX="$ROOT/aarch64-apple-darwin24.5.0/gfortran/Darwin"
BUILD="$SRC/build-develop-$(date +%Y%m%d)"
INSTALL="$SRC/install-develop-$(date +%Y%m%d)"

env PATH="$SANEPATH" cmake -S "$SRC" -B "$BUILD" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL" \
  -DCMAKE_Fortran_FLAGS:STRING="-I/opt/homebrew/opt/openmpi/include -I/opt/homebrew/include" \
  -DESMFMKFILE:PATH="$PREFIX/lib/esmf.mk" \
  -DCMAKE_PREFIX_PATH:PATH="$PREFIX"

env PATH="$SANEPATH" cmake --build "$BUILD" --target install -j 4
```

Known-good local example from 2026-05-20:

```bash
env PATH="$SANEPATH" cmake -S /Users/amfox/Desktop/GEOSldas_develop/GEOSldas \
  -B /Users/amfox/Desktop/GEOSldas_develop/GEOSldas/build-develop-20260520 \
  -DCMAKE_INSTALL_PREFIX=/Users/amfox/Desktop/GEOSldas_develop/GEOSldas/install-develop-20260520 \
  -DCMAKE_Fortran_FLAGS:STRING="-I/opt/homebrew/opt/openmpi/include -I/opt/homebrew/include" \
  -DESMFMKFILE:PATH=/Users/amfox/Desktop/GEOSldas_develop/aarch64-apple-darwin24.5.0/gfortran/Darwin/lib/esmf.mk \
  -DCMAKE_PREFIX_PATH:PATH=/Users/amfox/Desktop/GEOSldas_develop/aarch64-apple-darwin24.5.0/gfortran/Darwin

env PATH="$SANEPATH" cmake --build /Users/amfox/Desktop/GEOSldas_develop/GEOSldas/build-develop-20260520 \
  --target install -j 4
```

## Errors These Settings Fix

If OpenMPI include is missing from Fortran flags, builds can fail with:

- `Fatal Error: mpif.h: No such file or directory`

The `-I/opt/homebrew/opt/openmpi/include` setting resolves this.

If the Homebrew HDF5 Fortran module include path is missing, builds can fail late in
`GEOSlandassim_GridComp` with:

- `Fatal Error: Cannot open module file 'hdf5.mod' for reading`

The `-I/opt/homebrew/include` setting resolves this for the Homebrew HDF5 install.

If CMake cannot find ESMF or GFTL/GFTL_SHARED, point it at the local baselibs-style
prefix with:

- `-DESMFMKFILE:PATH="$PREFIX/lib/esmf.mk"`
- `-DCMAKE_PREFIX_PATH:PATH="$PREFIX"`

## Successful Build Check

After install, these files should exist in the install tree:

```bash
find "$INSTALL/bin" -maxdepth 1 -type f \( \
  -name 'GEOSldas.x' -o \
  -name 'preprocess_ldas.x' -o \
  -name 'ldas_setup' -o \
  -name 'mk_GEOSldasRestarts.x' -o \
  -name 'mwrtm_bin2nc4.x' \
\) -print
```

Warnings from older Fortran constructs, MPI argument mismatches, missing F2PY3,
missing ImageMagick/LaTeX, and Homebrew HDF5/libaec are expected in this local build.
