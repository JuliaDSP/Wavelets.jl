using Wavelets
using PyPlot
supertype(OrthoFilter)
dt=0.01;
t=1e-5:dt:40;
Y=0.1*randn(length(t))+sin(2*π*5*t);

wave,period,scale,coi=Wavelets.cwtft(Y,dt);

f=1./period;
ax=subplot(111);
cf=ax[:contourf](t,f,log(float(abs(wave).^2)),50);
cb=colorbar(cf);
cb[:set_label]("Wavelet Power");
ax[:fill_between](t,1./coi, color="white", alpha=0.5);
yscale("log");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
ylim([minimum(f),maximum(f)]);
