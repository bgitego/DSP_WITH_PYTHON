import matplotlib.pyplot as plt
import numpy as np
import cmath as cm
print('DSP Tutorial Implementation of The Discrete fourier transform')

N = 128                         # Number of Sample
N_MINUS_1 =  N-1               # Used to Define Array access 0 Index Max Value
N_OVER_2 = int((N/2))              # Used to Define Array Acccess 0 Index Max Value
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

f_hz1 = 1.0 #Frequency In Hz
f_hz2 = 10
f_hz3 = 100

A = 1    #Amplitude
B = 5
C = 2.5
s = (A * np.sin(f_hz1*t)) #+ (B *  np.cos(f_hz2*t)) #+ (C*np.sin(f_hz3*t))  #Signal
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
#INVERSE DFT Function in Rectangular 

def INV_DFT(REX,IMX,F_RANGE,SAMPLE_SIZE):
	
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
#Calculating the DFT of a sign signal

def DFT(F_XX,F_RANGE,SAMPLE_SIZE):
	
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

XX = INV_DFT(REX,IMX,N_OVER_2_P_1,N)
F_FRQ = DFT(F_XX,N_OVER_2_P_1,N)
POL = REC_2_POL(F_FRQ['F_REX'],F_FRQ['F_IMX'],N_OVER_2_P_1)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calculation of the FFT

#FFT Constant
CREX = np.zeros(N)   #REX[] Hold the Complex Real Part of the Frequency Domain
CIMX = np.zeros(N)   #IMX[] Hold the Complex Imaginary part of the frequency domain
CMAG = np.zeros(N)
CPHASE = np.zeros(N)
CF_XX =  np.zeros(N)            #XX[] Derived From Complex FFT: The Time Domain Signal
CF_REX = np.zeros(N_OVER_2_P_1) #REX[] Derived From Complex FFT: Hold Real Part of the Complex Frequency Domain
CF_IMX = np.zeros(N_OVER_2_P_1) #IMX[] Derived From Complex FFT: Hold Imaginary part of the Complex Frequency Domain

#Negative Frequency Generation For Inverse FFT

def NEG_HZ_GEN(REAL_HZ,IMG_HZ,SAMPLE_SIZE):

	COMPLEX_RHZ = np.zeros(SAMPLE_SIZE)	#Create Temporary Complex Frequency Domain 
	COMPLEX_RHZ[0:(int((SAMPLE_SIZE/2)+1))] = REAL_HZ #Copy Real data into Compex Frequency Domain
	
	COMPLEX_IHZ = np.zeros(SAMPLE_SIZE)
	COMPLEX_IHZ[0:(int((SAMPLE_SIZE/2)+1))] = IMG_HZ
	
	for K in range((int((SAMPLE_SIZE/2)+1)),SAMPLE_SIZE,1):
		COMPLEX_RHZ[K] =  COMPLEX_RHZ[SAMPLE_SIZE-K]
		COMPLEX_IHZ[K] = -COMPLEX_IHZ[SAMPLE_SIZE-K]			
	
	return COMPLEX_RHZ,COMPLEX_IHZ

CREX,CIMX = NEG_HZ_GEN(REX,IMX,N)
print (REX)
print (IMX)

print (CREX)
print (CIMX)
plt.figure(1)
plt.subplot(411)
plt.plot(REX)
plt.xlabel('Real Frequency')
plt.ylabel('Frequency')
plt.title('Magnitude')
plt.grid(1)

plt.subplot(412)
plt.plot(IMX)
plt.title('Imaginary Frequency')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.grid(1)

plt.subplot(413)
plt.plot(CREX)
plt.title('Complex Real Hz')
plt.xlabel('Frequency')
plt.ylabel('Phase(Radians)')
plt.grid(1)

plt.subplot(414)
plt.plot(CIMX)
plt.title('Complex Imaginary Hz')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.grid(1)


plt.show()
