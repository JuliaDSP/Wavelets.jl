using Wavelets
using Base.Test

print("util ...\n")

@test detailn(0) == 1
@test detailn(1) == 2
@test detailindex(0,1) == 2

j = 5
@test detailn(j)+1 == detailrange(j)[1]
@test detailindex(j,3) == detailrange(j)[3]
J = 7
@test nscales(2^J) == J
@test detailn(maxlevel(2^J)) == 2^(J-1)
@test nscales(rand(2^J)) == J
@test maxlevel(rand(2^J)) == J-1

@test tl2level(512,1) == 8
@test level2tl(512,9) == 0
@test tl2level(rand(2^9),1) == 8
@test level2tl(rand(2^9),9) == 0


# UTILITY FUNCTIONS

@test mirror([1]) == [1]
@test mirror([2,3]) == [2,-3]
@test mirror([2,3,4]) == [2,-3,4]
@test mirror([4.9,5,6,7]) == [4.9,-5,6,-7]

@test upsample([1]) == [1,0]
@test upsample([1],1) == [0,1]
@test upsample([1,2]) == [1,0,2,0]
@test upsample([1,2],1) == [0,1,0,2]
@test downsample([1,2]) == [1]
@test downsample([1,2],1) == [2]
@test downsample([1,2,3,4]) == [1,3]
@test downsample([1,2,3,4],1) == [2,4]

@test iscube(rand(4)) == true
@test iscube(rand(4,4,4,4)) == true
@test iscube(rand(4,5)) == false
@test isdyadic(rand(4)) == true
@test isdyadic(rand(32,32)) == true
@test isdyadic(rand(6)) == false

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
y = randn(n)
for sh in (-43, -5, 0, 1, 28, 41)
    @test circshift(x,sh) == circshift!(copy(x),sh)
    @test circshift(x,-sh) == circshift!(copy(x),-sh)
    @test x == circshift!(circshift!(copy(x),sh),-sh)
    @test circshift(x,sh) == circshift!(similar(x),x,sh)
    @test circshift(x,-sh) == circshift!(similar(x),x,-sh)
end

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

# with strides
a0 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
ia = 3
inca = 2
@test split!(similar(a0), copy(a0), ia, inca, 2)[1:2] == [3,5]
@test split!(similar(a0), copy(a0), ia, inca, 4)[1:4] == [3,7,5,9]
@test split!(similar(a0), copy(a0), 1, 3, 4)[1:4] == [1,7,4,10]
@test split!(similar(a0), copy(a0), 1, 1, 8)[1:8] ==  split!(similar(a0), copy(a0), 8)[1:8]

@test merge!(similar(a0), ia, inca, copy(a0), 2)[ia:inca:ia+(2-1)*inca] == [1,2]
@test merge!(similar(a0), ia, inca, copy(a0), 4)[ia:inca:ia+(4-1)*inca] == [1,3,2,4]
ia = 1
inca = 3
@test merge!(similar(a0), ia, inca, copy(a0), 4)[ia:inca:ia+(4-1)*inca] == [1,3,2,4]
@test merge!(similar(a0), 1, 1, copy(a0), 8)[1:8] ==  merge!(similar(a0), copy(a0), 8)[1:8]




