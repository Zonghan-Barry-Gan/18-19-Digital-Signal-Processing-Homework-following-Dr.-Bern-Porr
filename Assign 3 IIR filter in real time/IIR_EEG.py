

import numpy as np 

import matplotlib.pyplot as plt 

import matplotlib.animation as animation 

from scipy import signal as signal 

import pyusbdux as c 

from scipy import signal as signal 

  

# Realtime oscilloscope at a sampling rate of 50Hz 

# It displays analog channel 0. 

# You can plot multiple chnannels just by instantiating 

# more RealtimePlotWindow instances and registering 

# callbacks from the other channels. 

class IIRfilter: 

    def __init__(self,sos): 

        sos1=np.split(sos,2) 

  

        self.b10,self.b11,self.b12=sos1[0] 

        rubbin,self.a11,self.a12=sos1[1] 

  

  

        self.buffer01=0 

        self.buffer02=0 

          

  

  

          

    def filter(self, x):  

      

        acc_input1 = x - self.buffer01*self.a11- self.buffer02*self.a12 

        #acc_output is input multipied by IIR parts 

        acc_output1 = acc_input1*self.b10 + self.buffer01*self.b11+ self.buffer02*self.b12 

        #buffer of first filter input 

        self.buffer02 = self.buffer01 

        self.buffer01 = acc_input1 

        #buffer of first filter output 

  

      

  

  

  

        return acc_output1 

  

  

#!/usr/bin/python3 

  

class IIR2filter: 

    def __init__(self,_b0,_b1,_a1,_a2): 

        self.b0=_b0 

        self.b1=_b1 

        self.a1=_a1 

        self.a2=_a2 

          

  

        self.buffer1=0 

        self.buffer2=0 

          

    def filter(self, x): #x is input signal 

        #take signal and subtract buffers and a coefficients because -ve 

        acc_input = x - self.buffer1*self.a1 - self.buffer2*self.a2 

        #acc_output is input multipied by IIR parts 

        acc_output = acc_input*self.b0 + self.buffer1*self.b1 

        self.buffer2 = self.buffer1 

        self.buffer1 = acc_input 

        return acc_output #gives output signal 

  

# 

# Creates a scrolling data display 

class RealtimePlotWindow: 

  

    def __init__(self): 

        # create a plot window 

        self.fig, self.ax = plt.subplots() 

        # that's our plotbuffer 

        self.plotbuffer = np.zeros(500) 

        self.plotbuffer2 = np.zeros(500) 

        self.plotbufferfilter=np.zeros(500) 

        # create an empty line 

        self.line, = self.ax.plot(self.plotbuffer) 

        self.line2, = self.ax.plot(self.plotbuffer) 

        self.linefilter, = self.ax.plot(self.plotbufferfilter) 

        # axis 

        self.ax.set_ylim(-10, 10) 

        # That's our ringbuffer which accumluates the samples 

        # It's emptied every time when the plot window below 

        # does a repaint 

#        self.ringbuffer = [] 

#        # start the animation 

#        self.iffilter=0 

        self.ani1 = animation.FuncAnimation(self.fig, self.update,self.data_gen, interval=500) 

#        self.iffilter=1 

#        self.ani2 = animation.FuncAnimation(self.fig, self.updatefilter,self.data_gen, interval=500) 

          

        f=0.1 #normalised freq 

        q=10 #q factor 

          

        si = np.complex(-np.pi/f*q, np.pi/f/np.sqrt(4-(1/q**2))) 

          

        b0 = 1 

        b1 = -1 

        a1 = np.real(-(np.exp(si)+np.exp(np.conjugate(si)))) 

        a2 = np.exp(2*np.real(si)) 

          

          

        self.f2 =IIR2filter(b0,b1,a1,a2) # created our filter 

          

        #y=IIR2filter.filter(x) 

        self.order=4 

        sosz=signal.butter(self.order,[0.1/250*2,48/250*2],'bandpass', output = 'sos') 

        self.fltset=[]     

        for i in range(0,self.order): 

            fltseti=IIRfilter(sosz[i]) 

            self.fltset.append(fltseti) 

  

    def data_gen(self): 

        #endless loop which gets data 

        while True: 

            data = np.zeros(0) 

            filterdata1 = np.zeros(0) 

            filterdata2 = np.zeros(0) 

            while (not c.hasSampleAvilabale() == 0): 

                sample = c.getSampleFromBuffer() 

                filtersample2=self.f2.filter(sample[4]) 

                  

  

                filtersample1=sample[4] 

                for i in range(0,self.order): 

                    filtersample1=self.fltset[i].filter(filtersample1) 

                filterdata1 = np.append(filterdata1,filtersample1) 

                filterdata2 = np.append(filterdata2,filtersample2)             

                data = np.append(data,sample[4]) 

#                self.data=data 

                self.data2=filterdata2 

                self.f1data= filterdata1 

                yield data 

  

    # updates the plot 

    # receives the data from the generator below 

    def update(self,data): 

        self.plotbuffer=np.append(self.plotbuffer,data) 

        # only keep the 500 newest ones and discard the old ones 

        self.plotbuffer=self.plotbuffer[-500:] 

        # set the new 500 points of channel 9 

        self.line.set_ydata(self.plotbuffer) 

  

# 

        self.plotbufferfilter=np.append(self.plotbuffer,self.f1data) 

        # only keep the 500 newest ones and discard the old ones 

        self.plotbufferfilter=self.plotbufferfilter[-500:] 

        # set the new 500 points of channel 9 

        self.linefilter.set_ydata(self.plotbufferfilter) 

# 

# 

        self.plotbuffer2=np.append(self.plotbuffer2,self.data2) 

#        # only keep the 500 newest ones and discard the old ones 

        self.plotbuffer2=self.plotbuffer2[-500:] 

# 

#        # set the new 500 points of channel 9 

        self.line2.set_ydata(self.plotbuffer2) 

        return 

      

    # this checks in an endless loop if there is data in the ringbuffer 

    # of there is data then emit it to the update funnction above 

  

  

  

# Create an instance of an animated scrolling window 

# To plot more channels just create more instances and add callback handlers below 

  

  

# sampling rate: 50Hz 

#samplingRate = 50 

# 

## create a 2nd order Butterworth lowpass with 

## cutoff frequency at = 10 Hz and a gain of 5 

#cutoff = 10 

#gain = 5 

#lp_b, lp_a = signal.butter(2, cutoff/samplingRate*2.0) 

#lp_b = lp_b * gain 

#lp_z = signal.lfiltic(lp_b, lp_a, [0]) 

  

# our callback where we filter the data 

  

  

#realtimePlotWindow1.addData(data) 

  

  

  

# show the plot and start the animation 

  

# 

## receives the data from the generator below 

#def update(data): 

#    global plotbuffer 

#    plotbuffer=np.append(plotbuffer,data) 

#    # only keep the 500 newest ones and discard the old ones 

#    plotbuffer=plotbuffer[-500:] 

#    # set the new 500 points of channel 9 

#    line.set_ydata(plotbuffer) 

#    return line, 

# 

## this checks in an endless loop if there is data in the ringbuffer 

## of there is data then emit it to the update funnction above 

# 

  

  

c.open() 

c.start(8,250) 

  

realtimePlotWindowRaw = RealtimePlotWindow() 

  

  

  

  

# show it 

plt.show() 

  

c.stop() 

c.close() 

  

print("finished") 
