#!/usr/local/bin/bash

usage(){
	echo -e "\nUsage: $0 <directory>
	<directory> - name of the directory
"
	exit 1;
}

if [[ $# -lt 1 ]]; then
	usage
	exits;
fi

# deflation_level 1-9, higher level (e.g. 9) takes longer 
deflate_level=5

for file in $(find $1 -name "*.nc4" -type f); do
	#echo "Processing $file..."
	ncks -L $deflate_level -O $file $file.z
        mv $file.z $file
done;

