#!/bin/bash

usage="Usage: `basename $0` "
if [ $# -ne 0 ]; then 		# variable supplied?
	echo $usage 1>&2
	exit 1 
fi


# warm up
julia bm_fwt.jl > /dev/null

julia bm_fwt.jl
julia bm_lsfwt.jl
julia bm_fft.jl

exit


