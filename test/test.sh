#!/bin/bash

function md5test() {
	file=$1
	md5=$2
	test=$3
	echo "$(md5sum $file | awk '{print $1}')" -eq "$md5"
	if [[ "$(md5sum $file | awk '{print $1}')" == "$md5" ]]
	then
		echo "$file OK"
		return 0
	else
		echo "Test failed: $file's md5 is not $md5"
		#exit -1
	fi
}


# test with non-zipped fastq
$1 $2/bcd.fasta $2/r1.fastq
md5test $2/ 81ac378802a56ca9f43e091934729e4a
$1 $2/bcd.fasta $2/r1.fastq


