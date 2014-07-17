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

n = 64
@test wcount(randn(n)) == n
@test wcount(rand(n),0) == n
@test wcount(rand(n),1.01) == 0
@test wcount(rand(n).-2.5,1.01) == n
@test wcount([0,2,6,7,8],6.6) == 2
@test wcount([10,-11,6,7,8.0,-5,-8,0],3,level=2) == 3
@test wcount(randn(n,n)) == n*n
@test wcount(rand(n,n),1.01) == 0
@test wcount([-1 2;3 4],2.5) == 2

n = 32
x = randn(n)
@test circshift(x,5) == circshift!(copy(x),5)
@test circshift(x,-5) == circshift!(copy(x),-5)
@test x == circshift!(circshift!(copy(x),29),-29)
@test x == circshift!(circshift!(copy(x),292),-292)

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





