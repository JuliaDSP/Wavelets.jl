#!/bin/bash

usage="Usage: `basename $0` pathtojulia"
if [ $# -ne 1 ]; then 		# variable supplied?
	echo $usage 1>&2
	exit 1
fi

results="results_1d.txt"
echo "" >"$results"
# warm up
$1 bm_dwt_filt.jl > /dev/null

$1 bm_dwt_filt.jl >>"$results"
$1 bm_dwt_ls.jl >>"$results"
$1 bm_fft.jl >>"$results"

exit
