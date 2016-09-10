#!/bin/bash

usage="Usage: `basename $0` pathtojulia"
if [ $# -ne 1 ]; then 		# variable supplied?
	echo $usage 1>&2
	exit 1
fi

results="results_1d_inp.txt"
echo "" >"$results"
# warm up
$1 bm_dwt_filt_inp.jl > /dev/null

$1 bm_dwt_filt_inp.jl >>"$results"
$1 bm_dwt_ls_inp.jl >>"$results"
$1 bm_fft_inp.jl >>"$results"

exit
