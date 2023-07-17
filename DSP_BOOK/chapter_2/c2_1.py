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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

t = np.linspace(0,2*np.pi,512, endpoint=True)
print(len(t))
f_hz1 = 3.0 #Frequency In Hz
f_hz2 = 15
A = 10.00 # Amplitude
s = (A * np.sin(2*np.pi*f_hz1*t)) + ((A/10) * np.sin(2*np.pi*f_hz2*t))  #Signal

hann_2 = np.hanning(len(s))

Y = np.fft.fft(s)

N_2 = int((len(Y)/2) + 1) # FFT Frequency Sepctrum(Number of Sample/2  + 1 )

plt.figure(1)
plt.plot(s*hann_2)
plt.xlabel('Time s')
plt.ylabel('Magnitude')
