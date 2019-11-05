import numpy as np
import matplotlib.pyplot as plt
 
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
     
 
def highestchannel(data):
        #finding highest channel
    highestchannel=0
    average=[]
    for i in range(1,4):
        average.append(np.average(np.abs(data[:,i])))
    highestchannel=np.argmax(average)+1
        #pick up the highest channel
    data=data[:,highestchannel]
    return data
 
class matchfit:
    def __init__(self,matflt,fs):
        self.matflt=matflt
        self.taps=len(matflt)
        self.heartbeat=firflt("heart_beat_detector",self.taps,fs,[],[],[],[],0,[])
        #install match filter into fir
        self.heartbeat.hn=matflt
        self.matchdata=[]
        self.fs=fs
        return
    def domatch(self,v):
        v=v
        nextmatchdata=self.heartbeat.dofilter(v)
        self.matchdata.append(nextmatchdata)
        return nextmatchdata
     
class hbtcounter:
    def __init__(self,matchdata,fs):
        self.fs=fs
        self.matchdata=matchdata
        hbtindex=[]
        hbttime=[]
        hbttimeinterval=[]
        temhbr=[]
        hbttimes=0
        index=[]
        mininterval=60/200*fs
        for i in range(2,size):
            if matchdata[i]>8 and ((matchdata[i]-matchdata[i-1])*(matchdata[i-1]-matchdata[i-2]))<0:
                #the condition consists of a threshold and a peak condition
                index.append(i)
                lastindex=index[len(index)-2]
                if (i-lastindex)<mininterval:
                    #this is a fault beat/counting the same beat multiple times, otherwise heart rate over 200, unlikely to occur in this experiment
                    continue
                else:
                    hbtindex.append(i)
                    hbttimes=hbttimes+1
                    time_of_this_hbt=float(i)/float(fs)
                    hbttime.append(time_of_this_hbt)
                    if hbttimes>=2:
                        this_hbt_timeinterval=hbttime[hbttimes-1]-hbttime[hbttimes-2]
                        hbttimeinterval.append(this_hbt_timeinterval)
                        this_hbr=60/this_hbt_timeinterval
                        temhbr.append(this_hbr)
        self.hbttime=hbttime
        self.temhbr=temhbr
 
 
#task 2.1
#this plot is commented, if uncommented, can help find matchfilter.
#plt.figure()
#plt.plot(np.linspace(900,2200,1300),cleandata2[900:2200],'r-')
#plt.show()
taps=1000
fs=1000
lthr=[]
stopband1=45
stopband2=55
rthr=[]
pass0=0
EEGfilter1=firflt("EEG_FILTER1",taps,fs,lthr,stopband1,stopband2,rthr,pass0,'hamming')
#load data
 
 
#filter the data off DC and 45-55HZ
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
#template
EEGfilter1.datalist=np.zeros(EEGfilter1.taps)
matflt=cleandata2[1200:1800]
 
#2.1.1 create match filter
plt.figure()
plt.title('Match Filter of Heart rate ECG')
plt.plot(matflt)
matflt.reverse()
plt.plot(matflt,'r-')
plt.show()
 
#2.1.2 filter the data, which will be then analysed for heart rate
data=np.loadtxt('ecg_2.dat')
data3=highestchannel(data)
 
#dofilter
cleandata3=[]
size=len(data3)
for i in range(0,size):
    nextcleandata=EEGfilter1.dofilter(data3[i])
    #print(nextcleandata)
    cleandata3.append(nextcleandata)
EEGfilter1.datalist=np.zeros(EEGfilter1.taps)
 
#remove the noisy start 0-7000
data=cleandata3[7000:]
size=len(data)
 
#2.1.3 in this class, the data is match-filtered with the match filter in real time way
#feed in data     
hbt=matchfit(matflt,fs)
matchdata=[]
for i in range(0,size):
    nextmatchdata=hbt.domatch(data[i])
    matchdata.append(nextmatchdata)
    #record heartbeat
hbt.heartbeat.datalist=np.zeros(hbt.heartbeat.taps)
#task 2.2
matchdata=np.square(matchdata)
plt.plot(matchdata[0:5000])
counthbt=hbtcounter(matchdata,fs)
#remove the first time of heart beat, which don't have heart rate
hbttime2=counthbt.hbttime[1:len(counthbt.hbttime)]
plt.figure()
plt.title('Heart rate- time relation')
plt.ylabel('Heart rate/(times/min)')
plt.xlabel('time/s')
 
plt.plot(hbttime2,counthbt.temhbr)
plt.show()
