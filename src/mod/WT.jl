module WT
export
    DiscreteWavelet,
    FilterWavelet,
    LSWavelet,
    OrthoFilter,
    GLS,
    wavelet
using ..Util
import Base.length
using Compat.LinearAlgebra
using Compat: ComplexF64, undef, rmul!

# TYPE HIERARCHY

abstract type DiscreteWavelet{T} end
#abstract ContinuousWavelet{T}
# discrete transforms via filtering
abstract type FilterWavelet{T} <: DiscreteWavelet{T} end
# discrete transforms via lifting
abstract type LSWavelet{T} <: DiscreteWavelet{T} end
# all wavelet types
#const WaveletTransformType = Union{DiscreteWavelet, ContinuousWavelet}

"""Get wavelet type name."""
function name(::DiscreteWavelet) end
"""Get wavelet filter length."""
function length(::FilterWavelet) end

struct FilterTransform end
struct LiftingTransform end
"""Transform by filtering."""
const Filter = FilterTransform()
"""Transform by lifting."""
const Lifting = LiftingTransform()


# BOUNDARY TYPES

abstract type WaveletBoundary end
# periodic (default)
struct PerBoundary <: WaveletBoundary end
# zero padding
#struct ZPBoundary <: WaveletBoundary end
# constant padding
#struct CPBoundary <: WaveletBoundary end
# and so on...

const Periodic = PerBoundary()
const DEFAULT_BOUNDARY = PerBoundary()


# WAVELET CLASSES

"""
The `WaveletClass` type has subtypes `OrthoWaveletClass`
and `BiOrthoWaveletClass`.

The `WT` module has for convenience constants defined named
as the class short name and optionally amended with a number
specifing the number of vanishing moments.

A class can also be explicitly constructed as e.g. `Daubechies{4}()`.

# Examples
`WT.db2`, `WT.haar`, `WT.cdf97`

**See also:** `WT.class`, `WT.name`, `WT.vanishingmoments`
"""
abstract type WaveletClass end
abstract type OrthoWaveletClass <: WaveletClass end
abstract type BiOrthoWaveletClass <: WaveletClass end

# Single classes
for (TYPE, NAMEBASE, MOMENTS) in (
        (:Haar, "haar", 1),
        (:Beylkin, "beyl", -1), # TODO moments
        (:Vaidyanathan, "vaid",-1), # TODO moments
        )
    @eval begin
        struct $TYPE <: OrthoWaveletClass end
        class(::$TYPE) = $(string(TYPE))
        name(::$TYPE) = string($NAMEBASE)
        vanishingmoments(::$TYPE) = $MOMENTS
    end
    CONSTNAME = Symbol(NAMEBASE)
    @eval begin
        const $CONSTNAME = $TYPE()                  # type shortcut
    end
end

# Parameterized classes
for (TYPE, NAMEBASE, RANGE) in (
        (:Daubechies, "db", 1:10),
        (:Coiflet, "coif", 2:2:8),
        (:Symlet, "sym", 4:10),
        (:Battle, "batt", 2:2:6),
        )
    @eval begin
        struct $TYPE{N} <: OrthoWaveletClass end
        class(::$TYPE) = $(string(TYPE))
        name(::$TYPE{N}) where N = string($NAMEBASE,N)
        vanishingmoments(::$TYPE{N}) where N = N
    end
    for NUM in RANGE
        CONSTNAME = Symbol(string(NAMEBASE, NUM))
        @eval begin
            const $CONSTNAME = $TYPE{$NUM}()        # type shortcut
        end
    end
end

# Parameterized BiOrtho classes
for (TYPE, NAMEBASE, RANGE1, RANGE2) in (
        (:CDF, "cdf", [9], [7]),
        )
    @eval begin
        struct $TYPE{N1, N2} <: BiOrthoWaveletClass end
        class(::$TYPE) = $(string(TYPE))
        name(::$TYPE{N1, N2}) where {N1, N2} = string($NAMEBASE,N1,"/",N2)
        vanishingmoments(::$TYPE{N1, N2}) where {N1, N2} = (N1, N2)
    end
    for i in length(RANGE1)
        CONSTNAME = Symbol(string(NAMEBASE,RANGE1[i],RANGE2[i]))
        @eval begin
            const $CONSTNAME = $TYPE{$RANGE1[$i],$RANGE2[$i]}()        # type shortcut
        end
    end
end



# IMPLEMENTATIONS OF FilterWavelet

"""
Wavelet type for discrete orthogonal transforms by filtering.

**See also:** `GLS`, `wavelet`
"""
struct OrthoFilter{T<:WaveletBoundary} <: FilterWavelet{T}
    qmf     ::Vector{Float64}        # quadrature mirror filter
    name    ::String                 # filter short name
end

function OrthoFilter(w::WC, ::T=DEFAULT_BOUNDARY) where {WC<:OrthoWaveletClass, T<:WaveletBoundary}
    name = WT.name(w)
    if WC <: Daubechies
        qmf = daubechies(vanishingmoments(w))
    else
        qmf = get(FILTERS, name, nothing)
        qmf == nothing && throw(ArgumentError("filter not found"))
    end
    # make sure it is normalized in l2-norm
    return OrthoFilter{T}(qmf./norm(qmf), name)
end

length(f::OrthoFilter) = length(f.qmf)
qmf(f::OrthoFilter) = f.qmf
name(f::OrthoFilter) = f.name

"""Scale filter by scalar."""
function scale(f::OrthoFilter{T}, a::Number) where T<:WaveletBoundary
    return OrthoFilter{T}(f.qmf.*a, f.name)
end

"""Quadrature mirror filter pair."""
function makeqmfpair(f::OrthoFilter, fw::Bool=true, T::Type=eltype(qmf(f)))
    scfilter, dcfilter = makereverseqmfpair(f, fw, T)
    return reverse(scfilter), reverse(dcfilter)
end

"""Reversed quadrature mirror filter pair."""
function makereverseqmfpair(f::OrthoFilter, fw::Bool=true, T::Type=eltype(qmf(f)))
    h = convert(Vector{T}, qmf(f))
    if fw
        scfilter = reverse(h)
        dcfilter = mirror(h)
    else
        scfilter = h
        dcfilter = reverse(mirror(h))
    end
    return scfilter, dcfilter
end

#struct BiOrthoFilter{T<:WaveletBoundary} <: FilterWavelet{T}
#    qmf1::Vector{Float64}       # quadrature mirror filter 1
#    qmf2::Vector{Float64}       # quadrature mirror filter 2
#    name::String                # filter short name
#    BiOrthoFilter(qmf1, qmf2, name) = new(qmf1, qmf2, name)
#end


# IMPLEMENTATIONS OF LSWavelet

abstract type StepType end
struct PredictStep <: StepType end
struct UpdateStep <: StepType end
const Predict = PredictStep()
const Update = UpdateStep()

struct LSStepParam{T<:Number}
    coef    ::Vector{T}        # lifting coefficients
    shift   ::Int              # + left shift, - right shift
end

struct LSStep{T<:Number}
    param::LSStepParam{T}
    steptype::StepType
end

function LSStep(st::StepType, coef::Vector{T}, shift::Int) where T
    return LSStep{T}(LSStepParam{T}(coef, shift), st)
end

length(s::LSStep) = length(s.param)
length(s::LSStepParam) = length(s.coef)

"""
Wavelet type for discrete general (bi)orthogonal transforms
by using a lifting scheme.

**See also:** `OrthoFilter`, `wavelet`
"""
struct GLS{T<:WaveletBoundary} <: LSWavelet{T}
    step    ::Vector{LSStep{Float64}}    # steps to be taken
    norm1   ::Float64           # normalization of scaling coefs.
    norm2   ::Float64           # normalization of detail coefs.
    name    ::String            # name of scheme
end

function GLS(w::WC, ::T=DEFAULT_BOUNDARY) where {WC<:WaveletClass, T<:WaveletBoundary}
    name = WT.name(w)
    schemedef = get(SCHEMES, name, nothing)
    schemedef == nothing && throw(ArgumentError("scheme not found"))
    return GLS{T}(schemedef[1], schemedef[2], schemedef[3], name)
end

name(s::GLS) = s.name

# IMPLEMENTATIONS OF ContinuousWavelet

# ...

# TRANSFORM TYPE CONSTRUCTORS

"""
    wavelet(c[, t=WT.Filter][, boundary=WT.Periodic])

Construct wavelet type where `c` is a wavelet class,
`t` is the transformation type (`WT.Filter` or `WT.Lifting`),
and `boundary` is the type of boundary treatment.

# Examples
```julia
wavelet(WT.coif6)
wavelet(WT.db1, WT.Lifting)
```

**See also:** `WT.WaveletClass`
"""
function wavelet end

wavelet(c::WaveletClass, boundary::WaveletBoundary=DEFAULT_BOUNDARY) = wavelet(c, Filter, boundary)
wavelet(c::OrthoWaveletClass, t::FilterTransform, boundary::WaveletBoundary=DEFAULT_BOUNDARY) = OrthoFilter(c, boundary)
wavelet(c::WaveletClass, t::LiftingTransform, boundary::WaveletBoundary=DEFAULT_BOUNDARY) = GLS(c, boundary)

# ------------------------------------------------------------

# Compute filters from the Daubechies class
# N is the number of zeros at -1
function daubechies(N::Int)
    @assert N > 0
    # Create polynomial
    C = Vector{Int}(undef, N)
    @inbounds for n = 0:N-1
        C[N-n] = binomial(N-1+n, n)
    end

    # Find roots in y domain (truncated binomial series; (1 - y)^{-N})
    Y = roots(C)

    # Find roots in z domain:
    # z + z^{-1} = 2 - 4*y
    # where y is a root from above
    Z = zeros(ComplexF64, 2*N-2)
    @inbounds for i = 1:N-1
        Yi = Y[i]
        d = 2*sqrt( Yi*Yi - Yi )
        y2 = 1 - 2*Yi
        Z[i] = y2 + d
        Z[i+N-1] = y2 -d
    end

    # Retain roots inside unit circle
    nr = 0  # count roots
    @inbounds for i = eachindex(Z)
        if abs(Z[i]) <= 1 + eps()
            nr += 1
        end
    end

    # Find coefficients of the polynomial
    # (1 + z)^N * \prod_i (z - z_i)
    R = Vector{ComplexF64}(undef, N + nr)
    @inbounds for i = 1:N
        R[i] = -1
    end
    k = N
    @inbounds for i = eachindex(Z)
        if abs(Z[i]) <= 1 + eps()
            k += 1
            R[k] = Z[i]
        end
    end
    HH = vieta( R )

    # Normalize coefficients
    rmul!(HH, 1/norm(HH))
    return real(HH)
end

# Compute roots of polynomial
# Input is a coefficient vector with highest powers first
function roots(C::AbstractVector)
    A = compan(C)
    return eigvals(A)
end

# Create companion matrix for a polynomial
# Input is a coefficient vector with highest powers first
function compan(C::AbstractVector)
    n = length(C)
    A = zeros(n-1, n-1)

    if n > 1
        @inbounds A[1,:] .= -C[2:end] ./ C[1]
        @inbounds A[2:n:end] .= 1
    end
    return A
end

# Vieta-like formula for computing polynomial coefficients from roots
# See
# http://www.mathworks.se/help/matlab/ref/poly.html
function vieta(R::AbstractVector)
    n = length( R )
    C = zeros(ComplexF64, n+1)
    C[1] = 1
    Ci::ComplexF64 = 0
    Cig::ComplexF64 = 0

    @inbounds for k = 1:n
        Ci = C[1]
        for i = 1:k
            Cig = C[i+1]
            C[i+1] = Cig - R[k] * Ci
            Ci = Cig
        end
    end
    return C
end


# scaling filters h (low pass)
# the number at end of a filter name is the
# number of vanishing moments of the mother wavelet function
# sources:
# http://statweb.stanford.edu/~wavelab/ (Orthogonal/MakeONFilter.m)
# http://www.mathworks.com/matlabcentral/fileexchange/5502-filter-coefficients-to-popular-wavelets
### https://github.com/nigma/pywt/blob/master/src/wavelets_coeffs.template.h
# name => qmf
const FILTERS = Dict{String, Vector{Float64}}(
# Haar filter
"haar" =>
[0.7071067811865475,0.7071067811865475]
,
# Daubechies filters, see daubechies()

# Coiflet filters
"coif2" =>
[-0.072732619513,0.337897662458,0.852572020212,0.384864846864,-0.072732619513,-0.015655728135]
,
"coif4" =>
[0.0163873364635998,-0.0414649367819558,-0.0673725547222826,0.3861100668229939,0.8127236354493977,0.4170051844236707,-0.0764885990786692,-0.0594344186467388,0.0236801719464464,0.0056114348194211,-0.0018232088707116,-0.0007205494453679]
,
"coif6" =>
[-0.0037935128644910,0.0077825964273254,0.0234526961418362,-0.0657719112818552,-0.0611233900026726,0.4051769024096150,0.7937772226256169,0.4284834763776168,-0.0717998216193117,-0.0823019271068856,0.0345550275730615,0.0158805448636158,-0.0090079761366615,-0.0025745176887502,0.0011175187708906,0.0004662169601129,-0.0000709833031381,-0.0000345997728362]
,
"coif8" =>
[0.0008923136685824,-0.0016294920126020,-0.0073461663276432,0.0160689439647787,0.0266823001560570,-0.0812666996808907,-0.0560773133167630,0.4153084070304910,0.7822389309206135,0.4343860564915321,-0.0666274742634348,-0.0962204420340021,0.0393344271233433,0.0250822618448678,-0.0152117315279485,-0.0056582866866115,0.0037514361572790,0.0012665619292991,-0.0005890207562444,-0.0002599745524878,0.0000623390344610,0.0000312298758654,-0.0000032596802369,-0.0000017849850031]
,
"coif10" =>
[-0.0002120808398259,0.0003585896879330,0.0021782363583355,-0.0041593587818186,-0.0101311175209033,0.0234081567882734,0.0281680289738655,-0.0919200105692549,-0.0520431631816557,0.4215662067346898,0.7742896037334738,0.4379916262173834,-0.0620359639693546,-0.1055742087143175,0.0412892087544753,0.0326835742705106,-0.0197617789446276,-0.0091642311634348,0.0067641854487565,0.0024333732129107,-0.0016628637021860,-0.0006381313431115,0.0003022595818445,0.0001405411497166,-0.0000413404322768,-0.0000213150268122,0.0000037346551755,0.0000020637618516,-0.0000001674428858,-0.0000000951765727]
,
# Symmlet filter
"sym4" =>
[0.0455703458960000,-0.0178247014420000,-0.1403176241790000,0.4212345342040000,1.1366582434079999,0.7037390686560000,-0.0419109651250000,-0.1071489014180000]
,
"sym5" =>
[0.0276321529580000,-0.0298424998690000,-0.2479513626130000,0.0234789231360000,0.8965816483800000,1.0230529668940000,0.2819906968540000,-0.0553441861170000,0.0417468644220000,0.0386547959550000]
,
"sym6" =>
[-0.0110318675090000,0.0024999220930000,0.0632505626600000,-0.0297837512990000,-0.1027249698620000,0.4779043713330000,1.1138927839260000,0.6944579729580000,-0.0683231215870000,-0.1668632154120000,0.0049366123720000,0.0217847003270000]
,
"sym7" =>
[0.0145213947620000,0.0056713426860000,-0.1524638718960000,-0.1980567068070000,0.4081839397250000,1.0857827098140000,0.7581626019640000,0.0246656594890000,-0.0700782912220000,0.0960147679360000,0.0431554525820000,-0.0178704316510000,-0.0014812259150000,0.0037926585340000]
,
"sym8" =>
[-0.0047834585120000,-0.0007666908960000,0.0448236230420000,0.0107586117510000,-0.2026486552860000,-0.0866536154060000,0.6807453471900000,1.0991066305370001,0.5153986703740000,-0.0734625087610000,-0.0384935212630000,0.0694904659110000,0.0053863887540000,-0.0211456865280000,-0.0004283943000000,0.0026727933930000]
,
"sym9" =>
[0.0019811937360000,0.0008765025390000,-0.0187693968360000,-0.0163033512260000,0.0427444336020000,0.0008251409290000,-0.0771721610970000,0.3376589236020000,1.0152597908320000,0.8730484073490000,0.0498828309590000,-0.2708937835030000,-0.0257864459300000,0.0877912515540000,0.0125288962420000,-0.0145155785530000,-0.0006691415090000,0.0015124873090000]
,
"sym10" =>
[-0.0006495898960000,0.0000806612040000,0.0064957283750000,-0.0011375353140000,-0.0287862319260000,0.0081528167990000,0.0707035675500000,-0.0452407722180000,-0.0502565400920000,0.5428130112130000,1.0882515305000000,0.6670713381540000,-0.1002402150310000,-0.2255589722340000,0.0164188694260000,0.0649509245790000,-0.0020723639230000,-0.0122206426300000,0.0001352450200000,0.0010891704470000]
,
# Battle-Lemarie filter
"batt2" =>
[-0.0000867523000000,-0.0001586010000000,0.0003617810000000,0.0006529220000000,-0.0015570100000000,-0.0027458800000000,0.0070644200000000,0.0120030000000000,-0.0367309000000000,-0.0488618000000000,0.2809310000000000,0.5781630000000000,0.2809310000000000,-0.0488618000000000,-0.0367309000000000,0.0120030000000000,0.0070644200000000,-0.0027458800000000,-0.0015570100000000,0.0006529220000000,0.0003617810000000,-0.0001586010000000,-0.0000867523000000]
,
"batt4" =>
[0.0001033070000000,-0.0001642640000000,-0.0002018180000000,0.0003267490000000,0.0003959460000000,-0.0006556200000000,-0.0007804680000000,0.0013308600000000,0.0015462400000000,-0.0027452900000000,-0.0030786300000000,0.0057993200000000,0.0061414300000000,-0.0127154000000000,-0.0121455000000000,0.0297468000000000,0.0226846000000000,-0.0778079000000000,-0.0354980000000000,0.3068300000000000,0.5417360000000000,0.3068300000000000,-0.0354980000000000,-0.0778079000000000,0.0226846000000000,0.0297468000000000,-0.0121455000000000,-0.0127154000000000,0.0061414300000000,0.0057993200000000,-0.0030786300000000,-0.0027452900000000,0.0015462400000000,0.0013308600000000,-0.0007804680000000,-0.0006556200000000,0.0003959460000000,0.0003267490000000,-0.0002018180000000,-0.0001642640000000,0.0001033070000000]
,
"batt6" =>
[0.0001011130000000,0.0001107090000000,-0.0001591680000000,-0.0001726850000000,0.0002514190000000,0.0002698420000000,-0.0003987590000000,-0.0004224850000000,0.0006355630000000,0.0006628360000000,-0.0010191200000000,-0.0010420700000000,0.0016465900000000,0.0016413200000000,-0.0026864600000000,-0.0025881600000000,0.0044400200000000,0.0040788200000000,-0.0074684800000000,-0.0063988600000000,0.0128754000000000,0.0099063500000000,-0.0229951000000000,-0.0148537000000000,0.0433544000000000,0.0208414000000000,-0.0914068000000000,-0.0261771000000000,0.3128690000000000,0.5283740000000000,0.3128690000000000,-0.0261771000000000,-0.0914068000000000,0.0208414000000000,0.0433544000000000,-0.0148537000000000,-0.0229951000000000,0.0099063500000000,0.0128754000000000,-0.0063988600000000,-0.0074684800000000,0.0040788200000000,0.0044400200000000,-0.0025881600000000,-0.0026864600000000,0.0016413200000000,0.0016465900000000,-0.0010420700000000,-0.0010191200000000,0.0006628360000000,0.0006355630000000,-0.0004224850000000,-0.0003987590000000,0.0002698420000000,0.0002514190000000,-0.0001726850000000,-0.0001591680000000,0.0001107090000000,0.0001011130000000]
,
# Beylkin filter
"beyl" =>
[0.0993057653740000,0.4242153608130000,0.6998252140570000,0.4497182511490000,-0.1109275983480000,-0.2644972314460000,0.0269003088040000,0.1555387318770000,-0.0175207462670000,-0.0885436306230000,0.0196798660440000,0.0429163872740000,-0.0174604086960000,-0.0143658079690000,0.0100404118450000,0.0014842347820000,-0.0027360316260000,0.0006404853290000]
,
# Vaidyanathan filter
"vaid" =>
[-0.0000629061180000,0.0003436319050000,-0.0004539566200000,-0.0009448971360000,0.0028438345470000,0.0007081375040000,-0.0088391034090000,0.0031538470560000,0.0196872150100000,-0.0148534480050000,-0.0354703986070000,0.0387426192930000,0.0558925236910000,-0.0777097509020000,-0.0839288843660000,0.1319716614170000,0.1350842271290000,-0.1944504717660000,-0.2634948024880000,0.2016121617750000,0.6356010598720000,0.5727977932110000,0.2501841295050000,0.0457993341110000]


)


# biortho filters
# name => (qmf1,qmf2)
const BIFILTERS = Dict{String, NTuple{2, Vector{Float64}}}(
# test
"test" =>
([0.7071067811865475,0.7071067811865475],
 [0.7071067811865475,0.7071067811865475])

)


# in matlab (d,p) -> (predict, update)
const SCHEMES = Dict{String,NTuple{3, Any}}(
# cdf 5/3 -> bior 2.2, cdf 9/7 -> bior 4.4
# Cohen-Daubechies-Feauveau [Do Quan & Yo-Sung Ho. Optimized median lifting scheme for lossy image compression.]
"cdf9/7" => ([  LSStep(Update,  [1.0,1.0]*1.5861343420604, 0),
                LSStep(Predict, [1.0,1.0]*0.05298011857291494, 1),
                LSStep(Update,  [1.0,1.0]*-0.882911075531393, 0),
                LSStep(Predict, [1.0,1.0]*-0.44350685204384654, 1)],
                1.1496043988603355,
                0.8698644516247099)
,

# Haar, same as db1
"haar" => ([LSStep(Predict, [-1.0], 0),
            LSStep(Update,  [0.5], 0)],
            0.7071067811865475,
            1.4142135623730951)
,
# Daubechies
"db1" => ([ LSStep(Predict, [-1.0], 0),
            LSStep(Update,  [0.5], 0)],
            0.7071067811865475,
            1.4142135623730951)
,
"db2" => ([ LSStep(Predict, [-1.7320508075688772], 0),
            LSStep(Update,  [-0.0669872981077807,0.4330127018922193], 1),
            LSStep(Predict, [1.0], -1)],
            0.5176380902050414,
            1.9318516525781364)

)




end
