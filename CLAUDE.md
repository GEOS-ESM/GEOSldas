# GEOSldas — Claude Code Guide

## Start here

Read `doc/README.Components.FortranDevMemory.md` before touching anything in `src/Components`.
It contains the full component topology, CMake target map, module inventory, and Fortran file index.
Run `scripts/refresh_components_memory.sh` to update it after significant component changes.

## Repo structure

- **Fixture** repo — source lives in external sub-repos managed by `mepo`.
- Sub-repos use an `@` prefix: `@env`, `@cmake`, `src/Components/@GEOSldas_GridComp`, `src/Components/@GEOSgcm_GridComp`, `src/Shared/@GMAO_Shared`, `src/Shared/@MAPL`.
- Each `@` sub-repo is its own git repo. Use `git -C <path>` to run git commands inside them.
- `src/Components/@GEOSldas_GridComp` is where nearly all LDAS Fortran work happens.
- `src/Components/@GEOSgcm_GridComp` is a partial checkout (land/surface subtree only).

## Build

```bash
# Full build (head node)
./parallel_build.csh              # Release → build/ + install/
./parallel_build.csh -debug       # Debug   → build-Debug/ + install-Debug/

# Incremental (after first build)
cmake --build build-develop-20260520 --target install -j4
cmake --build build-develop-20260520 --target GEOSlandassim_GridComp

# Manual CMake (if debugging build system)
mkdir build && cd build
cmake .. -DBASEDIR=$BASEDIR/Linux -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_INSTALL_PREFIX=../install
make -j6 install
```

Never do in-source builds — `CMakeLists.txt` rejects them.

## Running experiments

```bash
install/bin/ldas_setup    # creates experiment dir with exeinp/batinp templates
# then: cd <exp>/run && sbatch lenkf.j
```

## Python

Default env for data work (cartopy, eccodes, xarray, netCDF4, numpy, matplotlib):
```
/Users/amfox/mamba/envs/regrid/bin/python
```
Other envs: `cygnss` env for CYGNSS work; `xr` env (xarray/netCDF4 only, no eccodes).

## Conventions

- Fortran/C++ linkage follows ESMF/MAPL patterns — preserve target names and legacy aliases.
- Don't change `@env` or other mepo-managed sub-repo layouts without a corresponding mepo/ecbuild update.
- Mirror existing flag patterns in `@env/build.csh` and `parallel_build.csh` when changing build flags.
- External deps: ESMF, MAPL, NetCDF, HDF5, GFTL_SHARED, ZLIB.

## Key files

| Purpose | Path |
|---------|------|
| Component map (read first) | `doc/README.Components.FortranDevMemory.md` |
| Build orchestration | `parallel_build.csh`, `@env/build.csh` |
| Top LDAS component | `src/Components/@GEOSldas_GridComp/GEOS_LdasGridComp.F90` |
| Data assimilation | `src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/` |
| ASCAT H SAF reader notes | `doc/README.ASCAT_HSAF_Reader.md` |
| Fortran linter/formatter | gfortran (linter), fprettify (formatter) via VS Code |
| Dev env script | `scripts/dev_env.sh` |

## Current work (as of 2026-06)

Branch `feature/amfox/ascat-hsaf-v8`: integrating H SAF ASCAT SSM products (H121 CDR v8
and H139 ICDR) into `clsm_ensupd_read_obs.F90`. New reader handles 12.5 km Fibonacci grid
netCDF files for Metop-A/B/C. See `doc/README.ASCAT_HSAF_Reader.md` for format details.

## Permissions

Global `~/.claude/settings.json` allows all Bash/Read/Edit/Write. Always confirm before:
- `rm`, `mv`, `cp`, `rsync` — destructive file ops
- `git push`, `git reset --hard`, `git rebase` — remote or history-rewriting
- `sbatch`, `scancel` — HPC job submission
- `pip/conda/mamba install` — environment changes
