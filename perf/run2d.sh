#!/bin/bash

usage="Usage: `basename $0` "
if [ $# -ne 0 ]; then 		# variable supplied?
	echo $usage 1>&2
	exit 1 
fi

results="results_2d.txt"
echo "" >"$results"
# warm up
julia bm_dwt2_filt.jl > /dev/null

julia bm_dwt2_filt.jl >>"$results"
julia bm_dwt2_ls.jl >>"$results"
julia bm_fft2.jl >>"$results"

exit


