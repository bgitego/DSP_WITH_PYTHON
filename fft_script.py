import csv
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd

#Calculation FFT of Accelerometer Data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

t_list=[]
x_list=[]
y_list=[]
z_list=[]

with open('accel.csv','rt') as csvfile:
	sp = csv.DictReader(csvfile)
	for row in sp:
		t_list.append(float(row['Second']))
		x_list.append(float(row['X']))
		y_list.append(float(row['Y']))
		z_list.append(float(row['Z']))


#Windowing using Hanning Approach

hann = np.hanning(len(x_list))


# Finding the FFT of the Z Axis
f_x = np.fft.fft(x_list*hann)

# Finding the Inverse of the FFT
invf_x = np.fft.ifft(f_x)

N_1_x = int((len(f_x)/2) + 1) # FFT Frequency Sepctrum(Number of Sample/2  + 1 )

# Finding the Frequency Axis

dt_x = float(t_list[2])-float(t_list[1])
fa_x = 1.0/dt_x

print ('dt_x=%.5fs (Sample Time)' % dt_x)
print ('fa_x=%.2fHz (Frequency)' % fa_x)

#scaling of the X axis based on the half sampling frequency
x_x_axis = np.linspace(0,int(fa_x/2),N_1_x,endpoint=True)

#Creating Simple Sine Wave
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

t = np.linspace(0,2*np.pi, 512, endpoint=True)
print(len(t))
f_hz1 = 3.0 #Frequency In Hz
f_hz2 = 15
A = 10.00 # Amplitude
s = (A * np.sin(2*np.pi*f_hz1*t)) + ((A/10) * np.sin(2*np.pi*f_hz2*t))  #Signal

hann_2 = np.hanning(len(s))

Y = np.fft.fft(s)

N_2 = int((len(Y)/2) + 1) # FFT Frequency Sepctrum(Number of Sample/2  + 1 )

#Scaling the FFT 
f_scl = (2/N_1_x)
y_scl = (2/N_2)

plt.figure(2)
plt.subplot(211)
plt.plot(s*hann_2)
plt.xlabel('Time s')
plt.ylabel('Magnitude')

plt.subplot(212)
plt.plot(abs(Y[:N_2]))
plt.xlabel('Frequency')
plt.ylabel('Magnitude')

plt.figure(1)

plt.subplot(212)
plt.plot(x_x_axis,(2*abs(f_x[:N_1_x]))/N_1_x)
plt.xlabel('Frequency')
plt.ylabel('Magnitude')

plt.subplot(211)
plt.plot(x_list*hann)
plt.xlabel('Time (s)')
plt.ylabel('Acceleration(g)')
plt.show()
