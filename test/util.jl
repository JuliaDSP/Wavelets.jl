using Wavelets
using Base.Test


@assert detailn(0) == 1
@assert detailn(1) == 2
@assert detailindex(0,1) == 2

j = 5
@assert detailn(j)+1 == detailrange(j)[1]
@assert detailindex(j,3) == detailrange(j)[3]
J = 13
@assert nscales(2^J) == J
@assert detailn(maxlevel(2^J)) == 2^(J-1)

