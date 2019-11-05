import numpy as np
import matplotlib.pyplot as plt
 
#TASK 1 - 1
#create python FIR filter class in realtime
class FIR_filter:
    def __init__(self,name,taps,fs):
        if taps % 2 ==1:
            taps = taps +1
        self.fs=fs
        self.name=name
        self.taps=taps
         
        self.ifr = np.ones(taps)
        k1 = int(taps/fs*45)
        k2 = int(taps/fs*55)
        self.ifr[k1:k2+1] = 0
        self.ifr[taps-k2:taps-k1+1]=0
        #create impulse response by doing ifft of ideal freq response
        impulse = np.real(np.fft.ifft(self.ifr))
        #swapping +ve and -ve parts
        h2 = np.zeros(taps) #creates empty impulse response
        h2[0:int(taps/2)] = impulse[int(taps/2):taps]
        h2[int(taps/2):taps] = impulse[0:int(taps/2)] #swapping
        #apply hamming window
        self.h2=h2
        self.window=np.hamming(taps)
        self.hn=self.window*self.h2
        self.datalist=np.zeros(taps)
        self.cleandata=[]
        return
     
    def dofilter(self,v):
        #shift 1 foward
        shift1=np.zeros(self.taps)
        shift1[1:self.taps]=self.datalist[0:(self.taps-1)]
        #feed in data
        shift1[0]=v
        self.datalist=shift1
         
        #do convolution
        yn=np.inner(self.datalist,self.hn)
 
        return yn    
  
 
#TASK 1.2
#creates python class which inherits FIR class from task 1-1 and cut off frequencies can be specified by user
class firflt(FIR_filter):
    def __init__(self,name,taps,fs,leftthreshold,stopband1,stopband2,rightthreshold,pass0,window):
        FIR_filter.__init__(self,name,taps,fs)
 
        self.nq=self.fs/2.0 #Nyquist Frequency
        if leftthreshold!=[]:
            self.lthr=float(leftthreshold)
            self.dolthr=1
        else:
            self.lthr=0
            self.dolthr=0
        if rightthreshold!=[]:
            self.rthr=float(rightthreshold)
            self.dorthr = 1
        else:
            self.rthr=self.nq
            self.dorthr = 0
        if stopband1!=[] and stopband2!=[]:
            self.stopband1=stopband1
            self.stopband2=stopband2
            self.dostopband = 1
        else:
            self.stopband2=self.rthr
            self.stopband1=self.lthr
            self.dostopband = 0
        self.pass0=pass0             
        self.window=window
 
        if not (self.lthr<=self.stopband1<=self.stopband2<=self.rthr<=self.nq):
            raise ValueError("Such filter does not exist")
 
        #create impulse response
        x=np.ones(taps)
        if self.dolthr:
            k1=int(float(self.lthr)/float(fs)*taps)
            x[0:k1+1]=0
            x[taps-k1+1:taps-1]=0
             
        if self.dorthr:
            k2=int(self.rthr/float(fs)*taps)
            x[(k2-1):(taps-k2+1)]=0
             
        if self.dostopband:
            k1=int(stopband1/float(fs)*taps)
            k2=int(stopband2/float(fs)*taps)
            x[k1:(k2+1)]=0
            x[(taps-k2):(taps-k1+1)]=0
             
        if not (pass0):
            x[0]=0
            x[taps-1]=0
         
        self.impres=x
        self.fm=np.linspace(0,fs,taps)
 
         
        #inverse fourier transform
        h1=np.fft.ifft(self.impres)
        self.h1=h1
         
        #shuffle
        h2=np.zeros(taps)
        h2[0:int(taps/2)]=h1[int(taps/2):taps]
        h2[int(taps/2):taps]=h1[0:int(taps/2)]
 
        self.h0=h2
         
        #do window
        if self.window==[]:
            self.wd=np.ones(taps)
            self.hn=h2            
        elif self.window=='hamming':
            self.wd=np.hamming(taps)
            self.hn=h2*self.wd
        elif self.window=='hanning':
            self.wd=np.hanning(taps)
            self.hn=self.wd*h2
        else:
            raise ValueError("such window is not in this class")
             
        self.hn=np.real(self.hn)
        fhn=np.abs(np.fft.fft(self.hn))
        self.datalist=np.zeros(taps)
        self.cleandata=[]
        return
     
 
       
     
#this function finds the channel with highest signal for task 1 & 2
def highestchannel(data): #finding highest channel
    highestchannel=0
    average=[]
    for i in range(1,4):
        average.append(np.average(np.abs(data[:,i])))
    highestchannel=np.argmax(average)+1
        #pick up the highest channel
    data=data[:,highestchannel]
    return data
 
 
#create filter for task 1-1
taps=1000
fs=1000
EEGfilter0=FIR_filter("EEG_FILTER0",taps,fs)
 
 
#plot the task 1-1 filter in dB in frequency domain
plt.figure()
plt.plot(np.linspace(0,fs,taps),20*np.log10(np.abs(np.fft.fft(EEGfilter0.hn))))
plt.title("Task 1-1: FIR filter in frequency domain")
plt.xlabel("f/HZ")
plt.ylabel("A/dB")
plt.show()
 
 
#load data
data=np.loadtxt('ecg_1.dat')
data2=highestchannel(data)
data2=data2
 
#filter the value in real time way ie. one by one
cleandata2=[]
size=len(data2)
# do filter of the data in real time, ie. one by one
for i in range(0,size):
    nextcleandata=EEGfilter0.dofilter(data2[i])
    #print(nextcleandata)
    cleandata2.append(nextcleandata)
#after using the filter, please clear the datalist
EEGfilter0.datalist=np.zeros(EEGfilter0.taps)
 
 
#plot the first 1000 points (first second) of data before and after filter
plt.figure()    
plt.plot(np.linspace(0,1000/fs,1000),data2[0:1000])
plt.plot(np.linspace(0,1000/fs,1000),cleandata2[0:1000],'r-')
plt.title("Task 1-1: ECG Before And After Filtering, 0&45-55HZ removed, first second")
plt.xlabel("t/s")
plt.ylabel("A/mV")
plt.show()
 
 
#task 1-2 with multiple group of specified frequency
 
#fir filter, low pass filter, Fc=35 HZ
lthr=0
stopband1=[]
stopband2=[]
rthr=35
pass0=0
#the coefficient can be specified in a very user-friendly way
EEGfilter2=firflt("EEG_FILTER2",taps,fs,lthr,stopband1,stopband2,rthr,pass0,'hamming')
 
 
#filter the value in real time way ie. one by one
cleandata2=[]
size=len(data2)
for i in range(0,size):
    nextcleandata=EEGfilter2.dofilter(data2[i])
    #print(nextcleandata)
    cleandata2.append(nextcleandata)
 
#after using the filter, please clear the datalist
EEGfilter2.datalist=np.zeros(EEGfilter2.taps)
 
 
plt.figure()    
plt.plot(np.linspace(0,1000/fs,1000),data2[0:1000])
plt.plot(np.linspace(0,1000/fs,1000),cleandata2[0:1000],'r-')
plt.title("Task 1-2: ECG Before And After Filtering, 0&35-HZ removed, first second")
plt.xlabel("t/s")
plt.ylabel("A/mV")
plt.show()
 
#next filter
#removing DC and 45-55HZ
lthr=[]
stopband1=45
stopband2=55
rthr=[]
pass0=0
EEGfilter1=firflt("EEG_FILTER1",taps,fs,lthr,stopband1,stopband2,rthr,pass0,'hamming')
 
#load data
#remove DC and 50Hz interference by targeting 45-55HZ
data=np.loadtxt('ecg_1.dat')
data2=highestchannel(data)
 
 
#filter the value in real time way ie. one by one
cleandata2=[]
size=len(data2)
for i in range(0,size):
    nextcleandata=EEGfilter1.dofilter(data2[i])
    #print(nextcleandata)
    cleandata2.append(nextcleandata)
 
#after using the filter, please clear the datalist
EEGfilter1.datalist=np.zeros(EEGfilter1.taps)
 
 
plt.figure()    
plt.plot(np.linspace(0,1000/fs,1000),data2[0:1000])
plt.plot(np.linspace(0,1000/fs,1000),cleandata2[0:1000],'r-')
plt.title("Task 1-2: ECG Before And After Filtering, 0&45-55HZ removed, first 1 second")
plt.xlabel("t/s")
plt.ylabel("A/mV")
plt.show()
