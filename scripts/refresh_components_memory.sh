#!/usr/bin/env bash

set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "${repo_root}"

components_root="src/Components"
memory_doc="doc/README.Components.FortranDevMemory.md"
file_inventory="doc/src_components_file_inventory.txt"
fortran_inventory="doc/src_components_fortran_inventory.txt"
cmake_inventory="doc/src_components_cmake_targets_inventory.txt"

mkdir -p doc

# Regenerate inventories.
find "${components_root}" -path '*/.git' -prune -o -path '*/.git/*' -prune -o -type f -print | sort > "${file_inventory}"
find "${components_root}" -path '*/.git' -prune -o -path '*/.git/*' -prune -o -type f -print \
  | grep -Ei '\.(F90|f90|F|f|F95|f95|F03|f03|F08|f08)$' | sort > "${fortran_inventory}"
find "${components_root}" -path '*/.git' -prune -o -path '*/.git/*' -prune -o -name CMakeLists.txt -print0 \
  | xargs -0 awk '/esma_add_library[[:space:]]*\(|ecbuild_add_executable[[:space:]]*\(|add_executable[[:space:]]*\(|add_library[[:space:]]*\(|esma_add_subdirectory[[:space:]]*\(|esma_add_subdirectories[[:space:]]*\(|add_subdirectory[[:space:]]*\(/ {print FILENAME ":" FNR ":" $0}' \
  > "${cmake_inventory}"

captured="$(date '+%Y-%m-%d %H:%M:%S %Z')"
fixture_commit="$(git rev-parse --short HEAD)"

if [[ -d "${components_root}/@GEOSldas_GridComp/.git" ]]; then
  ldas_commit="$(git -C "${components_root}/@GEOSldas_GridComp" rev-parse --short HEAD)"
else
  ldas_commit="unknown"
fi

if [[ -d "${components_root}/@GEOSgcm_GridComp/.git" ]]; then
  gcm_commit="$(git -C "${components_root}/@GEOSgcm_GridComp" rev-parse --short HEAD)"
else
  gcm_commit="unknown"
fi

file_count="$(wc -l < "${file_inventory}" | tr -d ' ')"
fortran_count="$(wc -l < "${fortran_inventory}" | tr -d ' ')"
cmake_count="$(find "${components_root}" -path '*/.git' -prune -o -path '*/.git/*' -prune -o -name CMakeLists.txt -print | wc -l | tr -d ' ')"

if [[ -f "${memory_doc}" ]]; then
  tmp_doc="$(mktemp)"
  awk \
    -v captured="${captured}" \
    -v fixture_commit="${fixture_commit}" \
    -v ldas_commit="${ldas_commit}" \
    -v gcm_commit="${gcm_commit}" \
    -v file_count="${file_count}" \
    -v fortran_count="${fortran_count}" \
    -v cmake_count="${cmake_count}" \
    '
      BEGIN {
        captured_line = "- Captured: `" captured "`"
        fixture_line = "- Fixture commit: `" fixture_commit "`"
        ldas_line = "- `@GEOSldas_GridComp` commit: `" ldas_commit "`"
        gcm_line = "- `@GEOSgcm_GridComp` commit: `" gcm_commit "`"
        file_line = "- Non-`.git` files under `src/Components`: `" file_count "`"
        fortran_line = "- Fortran source files (`.F90/.f90/.F/.f/.F95/.F03/.F08`): `" fortran_count "`"
        cmake_line = "- `CMakeLists.txt` files: `" cmake_count "`"
      }
      /^- Captured:/ { print captured_line; next }
      /^- Fixture commit:/ { print fixture_line; next }
      /^- `@GEOSldas_GridComp` commit:/ { print ldas_line; next }
      /^- `@GEOSgcm_GridComp` commit:/ { print gcm_line; next }
      /^- Non-`\.git` files under `src\/Components`:/ { print file_line; next }
      /^- Fortran source files \(`\.F90\/\.f90\/\.F\/\.f\/\.F95\/\.F03\/\.F08`\):/ { print fortran_line; next }
      /^- `CMakeLists.txt` files:/ { print cmake_line; next }
      { print }
    ' "${memory_doc}" > "${tmp_doc}"
  mv "${tmp_doc}" "${memory_doc}"
fi

echo "Refreshed:"
echo "  ${memory_doc}"
echo "  ${file_inventory}"
echo "  ${fortran_inventory}"
echo "  ${cmake_inventory}"
