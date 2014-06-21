#!/bin/bash

usage="Usage: `basename $0` "
if [ $# -ne 0 ]; then 		# variable supplied?
	echo $usage 1>&2
	exit 1 
fi


# warm up
julia benchmarks_fwt.jl > /dev/null

julia benchmarks_fwt.jl
julia benchmarks_iwt.jl
julia benchmarks_fft.jl

exit


