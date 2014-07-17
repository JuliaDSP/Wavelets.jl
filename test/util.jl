using Wavelets
using Base.Test


@test detailn(0) == 1
@test detailn(1) == 2
@test detailindex(0,1) == 2

j = 5
@test detailn(j)+1 == detailrange(j)[1]
@test detailindex(j,3) == detailrange(j)[3]
J = 13
@test nscales(2^J) == J
@test detailn(maxlevel(2^J)) == 2^(J-1)


@test mirror([1]) == [1]
@test mirror([2,3]) == [2,-3]
@test mirror([2,3,4]) == [2,-3,4]
@test mirror([4.9,5,6,7]) == [4.9,-5,6,-7]


a0 = [1,2]
exp = [1,2]
n = 2
@test split!(copy(a0)) == exp
@test split!(copy(a0), n, similar(a0)) == exp
@test split!(similar(a0), copy(a0), n) == exp
@test merge!(copy(a0)) == exp
@test merge!(copy(a0), n, similar(a0)) == exp
@test merge!(similar(a0), copy(a0), n) == exp

a0 = [-1,2,3.3,4]
exp = [-1,3.3,2,4]
n = 4
@test split!(copy(a0)) == exp
@test split!(copy(a0), n, similar(a0)) == exp
@test split!(similar(a0), copy(a0), n) == exp
a0 = [-1,3.3,2,4]
exp = [-1,2,3.3,4]
@test merge!(copy(a0)) == exp
@test merge!(copy(a0), n, similar(a0)) == exp
@test merge!(similar(a0), copy(a0), n) == exp

a0 = randn(128)
@test a0 == merge!(split!(a0))





