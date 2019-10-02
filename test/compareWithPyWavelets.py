import pywt
import numpy as np
import matplotlib.pyplot as plt
x = np.arange(512)
y = np.sin(2*np.pi*x/32)
coef, freqs = pywt.cwt(y, np.arange(1, 129), 'gaus1')
plt.matshow(coef)
plt.show()

pywt.wavelist(kind='continuous')


wav = pywt.ContinuousWavelet('morl')
wav.upper_bound
