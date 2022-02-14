#!/usr/local/bin/bash
#
# compress bit-shaved nc output

usage(){
	echo -e "\nUsage: $0 <directory>
	<directory> = name of output directory with bit-shaved nc files
"
	exit 1;
}

if [[ $# -lt 1 ]]; then
	usage
	exits;
fi

# deflation_level 1-9; higher level (e.g. 9) takes longer; recommended: 3 (based on L4_SM) 
deflate_level=3

for file in $(find $1 -name "*.nc4" -type f); do
	#echo "Processing $file..."
	ncks -L $deflate_level -O $file $file &
done;
wait
