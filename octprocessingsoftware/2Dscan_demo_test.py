"""
2Dscan_demo

Author = JdeWit
Date: 2019-08-15
email: J.deWit-1@tudelft.nl

This is a demo file to load and process 1D or 2D OCT data. 
1D data will still be plotted as a 2D dataset, the horizontal axis being the Ascan number
"""
# load modules
import numpy as np
import matplotlib.pyplot as plt
from time import time
import sys
sys.path.append('modules')
import ReadOCTfile
import DataProcessingOCT

plt.close('all')

# load dechirp to be used in the calculations
dechirp=np.fromfile('Chirp.data',np.float32)
# set data_folder and file name
data_folder="data\\"
file_name="Default_0002_Mode2D.oct"

#%%
t0=time()
# load the data from the file specified above
header,rawdata,spectrum,FOV=ReadOCTfile.OCTgetDataCombined(data_folder+file_name,spectrumindex=0)   
Ascanav=int(header.Ocity.Acquisition.IntensityAveraging.AScans.string)
# process the data into a Bscan image
image=DataProcessingOCT.calc_Bscan(rawdata,spectrum,dechirp,Ascanav,apodization='hanning',filters='none',objective='LK4')
t1=time()
print('It took ',t1-t0,' s to load and process file ',file_name)



#%%
OCTintensity=np.abs(image)**2
dB_im=10*np.log10(OCTintensity)
plt.figure(1)
plt.imshow(dB_im)

plt.figure(2, figsize=[10,5])
plt.subplot(121)
plt.title('Mean OCT intensity')
plt.plot(np.mean(OCTintensity, axis=1), '-r')

N_z = OCTintensity.shape[0] # Number of pixels in depth
N_t = OCTintensity.shape[1] # Number of pixels in time


for i in range(N_t):
    OCTintensity[:,i] = OCTintensity[:,i] - np.mean(OCTintensity, axis=1) # Subtract mean from signal

norm = np.linspace(N_z,1,N_z) # Normalization for unbiased correlation
g = np.zeros((N_z, N_t)) # Preallocate correlation
signal_FFT = np.fft.fft(np.concatenate((OCTintensity, np.zeros((N_z-1, N_t))),axis=0), axis=0) # FFT of signal
G = np.fft.ifft(np.conj(signal_FFT)*signal_FFT, axis=0) # Autocorrelation
G = np.real(G[:N_z, :])/norm[:, None] # Truncate autocorrelation
g = G/G[0,:] # Normalization

# define range over which to calculate average correlation coefficient
start=200
stop=250
ACmean=np.mean(g[:, start:stop], axis=1)

plt.subplot(122)
plt.plot(ACmean)


#%% plot the obtained image
#plt.figure(2)
#fig1=DataProcessingOCT.plot_Bscan_image(image,dBlevel=80,FOV=FOV,title='image (dB) file: '+file_name,plot_no=0)



  
    
