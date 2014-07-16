
using Wavelets
using Images  # for imread and imwrite

img = imread("lena.png")
x = permutedims(img.data, [ndims(img.data):-1:1])
L = 2
xts = wplotim(x, L, POfilter("db3"))

imwrite(xts, "transform2d_lena.png")

