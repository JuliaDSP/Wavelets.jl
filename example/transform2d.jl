
using Wavelets
using Images  # for imread and imwrite

img = imread("lena.tiff")
x = permutedims(img.data, [ndims(img.data):-1:1])
L = 2
xts = wplotim(x, L, wavelet(WT.db3, WT.Filter))

imwrite(xts, "lena_2d.png")
# convert with ImageMagick
run(`convert lena_2d.png transform2d_lena.jpg`)


