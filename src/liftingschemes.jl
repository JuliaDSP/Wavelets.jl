module LiftingSchemes
import ..POfilters: WaveletType
export WaveletLS, LSstep, GPLS

abstract WaveletLS <: WaveletType

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

function GPLS(name::String)
    step,norm1,norm2 = get(SCHEMES, name, nothing)
    step == nothing && error("scheme not found")
    return GPLS(step,norm1,norm2,name)
end

# class => namebase
FILTERC2N=(ASCIIString=>ASCIIString)[
"Coiflet" => "coif",
"Daubechies" => "db",
"Symmlet" => "sym",
"Battle" => "batt"
]

# in matlab (d,p) -> (p,u)
SCHEMES=(ASCIIString=>Tuple)[
# cdf 5/3 -> bior 2.2, cdf 9/7 -> bior 4.4
# Cohen-Daubechies-Feauveau [Do Quan & Yo-Sung Ho. Optimized median lifting scheme for lossy image compression.]
"cdf9/7" => ([ LSstep('u',1.5861343420604  *[1.0,1.0],0),
			LSstep('p',0.05298011857291494 *[1.0,1.0],1),
			LSstep('u',-0.882911075531393  *[1.0,1.0],0),
            LSstep('p',-0.44350685204384654*[1.0,1.0],1)],
            1.1496043988603355,
            0.8698644516247099)
,

# Haar, same as db1
"haar" => ([ LSstep('p',[-1.0],0),
            LSstep('u',[0.5],0)],
            0.7071067811865475,
            1.4142135623730951)
,
# Daubechies
"db1" => ([ LSstep('p',[-1.0],0),
            LSstep('u',[0.5],0)],
            0.7071067811865475,
            1.4142135623730951)
,
"db2" => ([ LSstep('p',[-1.7320508075688772],0),
            LSstep('u',[-0.0669872981077807,0.4330127018922193],1),
            LSstep('p',[1.0],-1)],
            0.5176380902050414,
            1.9318516525781364)

]           


end


