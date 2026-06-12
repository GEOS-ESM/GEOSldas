#!/usr/bin/env python3
"""
make_hsaf_flists.py — build per-day flist files for the GEOSldas H SAF ASCAT reader.

Searches metop_{a,b,c}/YYYY/MM/ under DATA_DIR for H121/H139 netCDF files,
extracts the sensing-start day from each filename, and writes one flist file
per satellite per day at:

    FLIST_DIR/Y{YYYY}/M{MM}/D{DD}/{flistname}.txt

Each flist file lists the bare filenames (no path) for that satellite and day,
sorted by sensing-start time.  GEOSldas read_obs_fnames (obs_dir_hier=1) will
prepend Y{YYYY}/M{MM}/ and combine with obs_param%path to reach the actual file.

Usage
-----
    python make_hsaf_flists.py --data-dir /path/to/data --flist-dir /path/to/flists

    # Override flist filename for a satellite (default shown):
    python make_hsaf_flists.py ... --name-a H121_METOPA.txt \
                                   --name-b H121_H139_METOPB.txt \
                                   --name-c H121_H139_METOPC.txt

    # Process only one satellite:
    python make_hsaf_flists.py ... --satellites a

Filename convention (H SAF H121/H139)
--------------------------------------
W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{A|B|C}-12.5km-H{121|139}_C_LIIB_{c14}_{s14}_{e14}____.nc
                                                                          ^
                                         {s14} = sensing start YYYYMMDDHHMMSS
                                         position: fname[-36:-22] (0-indexed from end)
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path


SATELLITES = {
    'a': {'subdir': 'metop_a', 'default_name': 'H121_METOPA.txt'},
    'b': {'subdir': 'metop_b', 'default_name': 'H121_H139_METOPB.txt'},
    'c': {'subdir': 'metop_c', 'default_name': 'H121_H139_METOPC.txt'},
}


def sensing_start(fname: str) -> str:
    """Return the 14-char sensing-start string from a bare H SAF filename."""
    return fname[-36:-22]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-dir', required=True,
                   help='Base directory containing metop_a/, metop_b/, metop_c/ subdirs')
    p.add_argument('--flist-dir', required=True,
                   help='Output directory for Y{YYYY}/M{MM}/D{DD}/*.txt flist files')
    p.add_argument('--satellites', default='abc',
                   help='Which satellites to process (any combination of a, b, c). Default: abc')
    p.add_argument('--name-a', default=SATELLITES['a']['default_name'],
                   help=f"Flist filename for Metop-A (default: {SATELLITES['a']['default_name']})")
    p.add_argument('--name-b', default=SATELLITES['b']['default_name'],
                   help=f"Flist filename for Metop-B (default: {SATELLITES['b']['default_name']})")
    p.add_argument('--name-c', default=SATELLITES['c']['default_name'],
                   help=f"Flist filename for Metop-C (default: {SATELLITES['c']['default_name']})")
    return p.parse_args()


def collect_files(sat_dir: Path) -> dict:
    """
    Walk sat_dir/YYYY/MM/*.nc and return {(YYYY, MM, DD): [sorted bare filenames]}.
    """
    by_day = defaultdict(list)
    if not sat_dir.is_dir():
        print(f'  WARNING: directory not found, skipping: {sat_dir}', file=sys.stderr)
        return by_day

    for nc_file in sorted(sat_dir.rglob('*.nc')):
        name = nc_file.name
        if len(name) < 36:
            print(f'  WARNING: filename too short to parse, skipping: {name}', file=sys.stderr)
            continue
        s14 = sensing_start(name)
        if not s14.isdigit():
            print(f'  WARNING: non-numeric sensing start "{s14}" in {name}, skipping',
                  file=sys.stderr)
            continue
        yyyy, mm, dd = s14[0:4], s14[4:6], s14[6:8]
        by_day[(yyyy, mm, dd)].append(name)

    # sort by sensing-start time within each day
    for key in by_day:
        by_day[key].sort(key=lambda f: sensing_start(f))

    return by_day


def write_flists(by_day: dict, flist_dir: Path, flist_name: str, dry_run: bool = False):
    n_written = 0
    for (yyyy, mm, dd), names in sorted(by_day.items()):
        out_dir = flist_dir / f'Y{yyyy}' / f'M{mm}' / f'D{dd}'
        out_file = out_dir / flist_name
        if not dry_run:
            out_dir.mkdir(parents=True, exist_ok=True)
            with open(out_file, 'w') as f:
                f.write('\n'.join(names) + '\n')
        else:
            print(f'  [dry-run] would write {out_file} ({len(names)} entries)')
        n_written += 1
    return n_written


def main():
    args = parse_args()
    data_dir  = Path(args.data_dir)
    flist_dir = Path(args.flist_dir)

    flist_names = {'a': args.name_a, 'b': args.name_b, 'c': args.name_c}

    for sat in sorted(set(args.satellites.lower()) & set('abc')):
        info     = SATELLITES[sat]
        sat_dir  = data_dir / info['subdir']
        fl_name  = flist_names[sat]

        print(f'Metop-{sat.upper()}: scanning {sat_dir} ...')
        by_day = collect_files(sat_dir)

        if not by_day:
            print(f'  no files found.')
            continue

        n_days  = len(by_day)
        n_files = sum(len(v) for v in by_day.values())
        print(f'  found {n_files} files across {n_days} days')

        n_written = write_flists(by_day, flist_dir, fl_name)
        print(f'  wrote {n_written} flist files -> {flist_dir}/**/{fl_name}')


if __name__ == '__main__':
    main()
