import matplotlib.pyplot as plt
import numpy as np
import cmath as cm
print 'DSP Tutorial Implementation of The Discrete fourier transform'

N = 512                        # Number of Sample
N_MINUS_1 =  N-1               # Used to Define Array access 0 Index Max Value
N_OVER_2 = (N/2)               # Used to Define Array Acccess 0 Index Max Value
N_OVER_2_P_1 = N_OVER_2 + 1    # Number of Frequency Samples

PI = np.pi                     # Definition of the PI Constant

XX =  np.zeros(N)              #XX[] Hold The Time Domain Signal
REX = np.zeros(N_OVER_2_P_1)   #REX[] Hold the Real Part of the Frequency Domain
IMX = np.zeros(N_OVER_2_P_1)   #IMX[] Hold the Imaginary part of the frequency domain
MAG = np.zeros(N_OVER_2_P_1)
PHASE = np.zeros(N_OVER_2_P_1)
F_XX =  np.zeros(N)            #XX[] Hold The Time Domain Signal
F_REX = np.zeros(N_OVER_2_P_1) #REX[] Hold the Real Part of the Frequency Domain
F_IMX = np.zeros(N_OVER_2_P_1) #IMX[] Hold the Imaginary part of the frequency domain

#Creating Simple Sine Wave
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

t = np.linspace(0,2*np.pi,N, endpoint=True)
#print(len(t))
#print 'number of samples', t

f_hz1 = 71.0 #Frequency In Hz
f_hz2 = 15
f_hz3 = 1034

A = 22.00    #Amplitude
B = 15.00
C = 45
s = (A * np.sin(f_hz1*t)) + (B *  np.cos(f_hz2*t)) + (C*np.sin(f_hz3*t))  #Signal
#print 'Length of Time Domain Signal: ',len(s)
#Windowing
hamm = np.hamming(len(t))
#s = s*hamm

#Find the discrete fourier transform of the sine wave
Y = np.fft.rfft(s)

#print np.array(Y).real
#print np.array(Y).imag
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Assign the real and Imaginary part of the frequency response to REX and IMX

REX = np.array(Y).real
IMX = np.array(Y).imag

#Assign the time domain signal to F_XX[]
F_XX = s

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#INVERSE FFT Function in Rectangular 

def INV_FFT(REX,IMX,F_RANGE,SAMPLE_SIZE):
	
	for K in range(0,F_RANGE,1):
		REX[K] =  REX[K]/(SAMPLE_SIZE/2)
		IMX[K] = -IMX[K]/(SAMPLE_SIZE/2)
	REX[0] = REX[0]/2
	REX[(F_RANGE-1)] = REX[(F_RANGE-1)]/2
	
	for K in range(0,F_RANGE,1):
		for I in range(0,SAMPLE_SIZE):
			XX[I] = XX[I] + REX[K]*np.cos(2*PI*K*I/N) + IMX[K]*np.sin(2*PI*K*I/N)
	return XX		
		
#Calculating the Inverst FFT

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calculating the FFT of a sign signal

def FFT(F_XX,F_RANGE,SAMPLE_SIZE):
	
	for K in range(0,(F_RANGE),1):
		for I in range(0,SAMPLE_SIZE,1):
			F_REX[K] = F_REX[K] + F_XX[I]*np.cos(2*PI*K*I/N)
			F_IMX[K] = F_IMX[K] - F_XX[I]*np.sin(2*PI*K*I/N)
	return {'F_REX':F_REX,'F_IMX':F_IMX}


def REC_2_POL(REX,IMX,F_RANGE):
	
	for K in range(0,F_RANGE,1):
		MAG[K] = np.sqrt(np.square(REX[K]) + np.square(IMX[K]))
		if (REX[K] == 0):
			REX[K] = 1E-20
		PHASE[K] = np.arctan(IMX[K]/REX[K])
		if (REX[K] < 0 and IMX[K] < 0):
			PHASE[K] = PHASE[K] - PI
		if (REX[K] < 0 and IMX[K] >= 0):
			PHASE[K] = PHASE[K] + PI 
		
	return {'MAG':MAG,'PHASE':PHASE}

XX = INV_FFT(REX,IMX,N_OVER_2_P_1,N)
F_FRQ = FFT(F_XX,N_OVER_2_P_1,N)
POL = REC_2_POL(F_FRQ['F_REX'],F_FRQ['F_IMX'],N_OVER_2_P_1)

print POL

plt.figure(1)
plt.subplot(411)
plt.plot(XX)
plt.xlabel('Time s')
plt.ylabel('Magnitude')
plt.title('Inverse FT Graph')
plt.grid(1)

plt.subplot(412)
plt.plot(s)
plt.title('Original Signal')
plt.xlabel('Time s')
plt.ylabel('Magnitude')
plt.grid(1)

plt.subplot(413)
plt.plot(np.array(Y).real)
plt.title('FFT Python')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.grid(1)

plt.subplot(414)
plt.plot(POL['MAG'])
plt.title('FFT Computed')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.grid(1)


plt.show()
