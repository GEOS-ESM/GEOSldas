#!/usr/local/bin/bash
#
# lossless compression of (bit-shaved) nc4 output, operates recursively in <directory>
#
#    usage: compress_bit-shaved_nc4.sh <directory>
#
# ---------------------------------------------------------------

usage(){
	echo -e "\nUsage: $0 <directory>
	<directory> = name of output directory with bit-shaved nc4 files
"
	exit 1;
}

if [[ $# -lt 1 ]]; then
	usage
	exits;
fi

# deflation_level 1-9; higher level (e.g. 9) takes longer; recommended: 3 (based on L4_SM) 
deflate_level=3

# By default, the ncks command is executed sequentially for each nc4 file
# in <directory>.  Modify this script if multi-threading is needed to run
# the script on a compute node. 

for file in $(find $1 -name "*.nc4" -type f); do
        #echo "Processing $file..."
        #                                               # For simple multi-threading:
        ncks -L $deflate_level -O $file $file           #  <-- comment out this line
###     ncks -L $deflate_level -O $file $file &         #  <-- uncomment this line
done;
###wait                                                 #  <-- uncomment this line

# ===================== EOF =======================================================
