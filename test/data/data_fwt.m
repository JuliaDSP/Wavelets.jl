
% in octave, using wavelab
% 2D
n=32;
L=0; % J-L in julia
x=1:n;
x(4:6)=[-9.9,-1.0,0];
x(25:29)=-5:-1;
x(15:17)=0.0
x
qmf=MakeONFilter('Daubechies',6);
y=FWT_PO(x,L,qmf)

% 2D
n=16;
L=0; % J-L in julia
x2=zeros(n,n);
x2(1:n)=1:n;
x2(4,5:7)=[2.5,1.1,-9.9];
x2(n,n-2)=4.2;
x2(5,12)=-3;
x2(10,4)=1;
x2(11:14,7:10)=[4 7 -2 -9
9 8 2 1
-5 -2.0 -8 9
1 8 -11 2];
qmf=MakeONFilter('Daubechies',6);
y2=FWT2_PO(x2,L,qmf)
