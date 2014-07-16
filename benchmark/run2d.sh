#!/bin/bash

usage="Usage: `basename $0` "
if [ $# -ne 0 ]; then 		# variable supplied?
	echo $usage 1>&2
	exit 1 
fi


# warm up
julia bm_fwt2.jl > /dev/null

julia bm_fwt2.jl
julia bm_iwt2.jl
julia bm_fft2.jl

exit


