#!/bin/bash

usage="Usage: `basename $0` "
if [ $# -ne 0 ]; then 		# variable supplied?
	echo $usage 1>&2
	exit 1 
fi

results="results_1d_inp.txt"
echo "" >"$results"
# warm up
julia bm_dwt_filt_inp.jl > /dev/null

julia bm_dwt_filt_inp.jl >>"$results"
julia bm_dwt_ls_inp.jl >>"$results"
julia bm_fft_inp.jl >>"$results"

exit


