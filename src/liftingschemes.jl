module LiftingSchemes
export WaveletLS, LSstep, GPLS

abstract WaveletLS

immutable LSstep
    stept::Char             # step type: 'p' for predict, 'u' for update
    coef::Vector            # lifting coefficients
    shift::Integer          # + left shift, - right shift
end

# general periodic lifting scheme
immutable GPLS <: WaveletLS
    step::Vector{LSstep}    # steps to be taken
    norm1::Real             # normalization of scaling coefs.
    norm2::Real             # normalization of detail coefs.
    name::String            # name of scheme
end



end


