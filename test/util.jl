using Wavelets
using Base.Test

print("util ...\n")

@test dyadicdetailn(0) == 1
@test dyadicdetailn(1) == 2
@test dyadicdetailindex(0,1) == 2

j = 5
@test dyadicdetailn(j)+1 == dyadicdetailrange(j)[1]
@test dyadicdetailindex(j,3) == dyadicdetailrange(j)[3]
J = 7
@test ndyadicscales(2^J) == J
@test dyadicdetailn(maxdyadiclevel(2^J)) == 2^(J-1)
@test ndyadicscales(rand(2^J)) == J
@test maxdyadiclevel(rand(2^J)) == J-1

@test tl2dyadiclevel(512,1) == 8
@test dyadiclevel2tl(512,9) == 0
@test tl2dyadiclevel(rand(2^9),1) == 8
@test dyadiclevel2tl(rand(2^9),9) == 0

n = 64
L = 2
@test maxtransformlevels(n) == ndyadicscales(n)
@test detailindex(n, L, 3) == dyadicdetailindex(tl2dyadiclevel(n,L), 3)
@test detailrange(n, L) == dyadicdetailrange(tl2dyadiclevel(n,L))
@test detailn(n, L) == dyadicdetailn(tl2dyadiclevel(n,L))


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
    @test circshift(x,sh) == Util.circshift!(copy(x),sh)
    @test circshift(x,-sh) == Util.circshift!(copy(x),-sh)
    @test x == Util.circshift!(Util.circshift!(copy(x),sh),-sh)
    @test circshift(x,sh) == Util.circshift!(similar(x),x,sh)
    @test circshift(x,-sh) == Util.circshift!(similar(x),x,-sh)
end

a0 = [1,2]
exp = [1,2]
n = 2
@test Util.split!(copy(a0)) == exp
@test Util.split!(copy(a0), n, similar(a0)) == exp
@test Util.split!(similar(a0), copy(a0), n) == exp
@test Util.merge!(copy(a0)) == exp
@test Util.merge!(copy(a0), n, similar(a0)) == exp
@test Util.merge!(similar(a0), copy(a0), n) == exp

a0 = [-1,2,3.3,4]
exp = [-1,3.3,2,4]
n = 4
@test Util.split!(copy(a0)) == exp
@test Util.split!(copy(a0), n, similar(a0)) == exp
@test Util.split!(similar(a0), copy(a0), n) == exp
a0 = [-1,3.3,2,4]
exp = [-1,2,3.3,4]
@test Util.merge!(copy(a0)) == exp
@test Util.merge!(copy(a0), n, similar(a0)) == exp
@test Util.merge!(similar(a0), copy(a0), n) == exp

a0 = randn(128)
@test a0 == Util.merge!(Util.split!(a0))

# with strides
a0 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
ia = 3
inca = 2
@test Util.split!(similar(a0), copy(a0), ia, inca, 2)[1:2] == [3,5]
@test Util.split!(similar(a0), copy(a0), ia, inca, 4)[1:4] == [3,7,5,9]
@test Util.split!(similar(a0), copy(a0), 1, 3, 4)[1:4] == [1,7,4,10]
@test Util.split!(similar(a0), copy(a0), 1, 1, 8)[1:8] ==  Util.split!(similar(a0), copy(a0), 8)[1:8]

@test Util.merge!(similar(a0), ia, inca, copy(a0), 2)[ia:inca:ia+(2-1)*inca] == [1,2]
@test Util.merge!(similar(a0), ia, inca, copy(a0), 4)[ia:inca:ia+(4-1)*inca] == [1,3,2,4]
ia = 1
inca = 3
@test Util.merge!(similar(a0), ia, inca, copy(a0), 4)[ia:inca:ia+(4-1)*inca] == [1,3,2,4]
@test Util.merge!(similar(a0), 1, 1, copy(a0), 8)[1:8] ==  Util.merge!(similar(a0), copy(a0), 8)[1:8]

n = 128
x = randn(n)
for L = 0:ndyadicscales(n)
    for st in (:full, :dwt)
        @test isvalidtree(x, maketree(n, L, st))
    end
end
tree = maketree(n, 4)
tree[5] = false
@test !isvalidtree(x, tree)
tree = maketree(n, 4)
tree[18] = true
@test isvalidtree(x, tree)
tree[7] = false
@test !isvalidtree(x, tree)
tree[7] = true
tree[9] = false
@test !isvalidtree(x, tree)

EE = Exception
@test_throws EE maketree(n, 4, :foo)


makewavelet(wavelet(WT.db2, WT.Filter))
for tf in ("Blocks", "HeaviSine")
    testfunction(8, tf)
end
EE = Exception
@test_throws EE testfunction(8, "sksksks")
