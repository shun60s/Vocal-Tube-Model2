#coding:utf-8

# Tube model check:
# Comparison computed frequency response to white noise input frequency response using FFT analysis.
#

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal # version > 1.2.0
from scipy.io.wavfile import write as wavwrite

class Class_Tube_Check(object):
    def __init__(self, tube):
        #
        self.tube= tube
        self.sr= tube.sr
        #
        self.freq_low=100
        self.freq_high=6000
        self.resolution=1
        self.f_list=np.linspace(self.freq_low, self.freq_high,int((self.freq_high - self.freq_low)/self.resolution + 1 ) )
        #
        self.H1_linear()
        self.get_peaks()
        #
        self.make_white_noise(sinpuku=0.01, length=1)
        self.yout=self.tube.process(self.xin)
        self.fft_ana(self.yout, N_sample=4096*4)
        self.draw()
        #
        if 1:
            self.save_wav(self.xin,'test_' + str(self.tube.num_of_tube) +'tube_xin.wav')
            self.save_wav(self.yout,'test_'+ str(self.tube.num_of_tube) +'tubu_yout.wav')
        
    def H1_linear(self,):
        self.f_amp=self.tube.fone(self.f_list * 2.0 * np.pi)
        self.amp1=np.log10(self.f_amp) * 20
        self.freq=self.f_list
        
    def get_peaks(self,):
        MIN_HIGH= 0.4   # ピークの最小高さ
        #MIN_DIS= 0.000001    # 最小の周辺距離
        MIN_WIDTH= 1 # 最小の周辺幅
        self.peaks, _ = signal.find_peaks(self.amp1, height= MIN_HIGH * max(self.amp1),width= MIN_WIDTH) #, distance= MIN_DIS * n , width= MIN_WIDTH * n)
        print('peaks', self.freq[self.peaks])
 
    def draw(self,):
        # draw
        fig = plt.figure()
        
        plt.plot(self.freq, self.amp1,'r')
        plt.plot(self.freq[self.peaks], self.amp1[self.peaks], "x")
        
        if 1:
            id0=np.where( self.y_freq > self.f_list[0])[0][0]
            id1=np.where( self.y_freq > self.f_list[-1])[0][0]
            plt.plot(self.y_freq[id0:id1], self.yf[id0:id1],'b')
            plt.title('frequency response ' + str(self.tube.num_of_tube) + ' : red computed vs blue white noise input')
        else:
            plt.title('frequency response')
        
        plt.xlabel('Hz')
        plt.ylabel('dB')
        plt.grid(which='both', axis='both')
        
        #
        fig.tight_layout()
        plt.show()

    def fft_ana(self, yin, N_sample=4096*4):
        N_sample_half= int(N_sample/2)
        x=yin[ int(len(yin)/2 - N_sample_half) : int(len(yin)/2 + N_sample_half) ]
        self.yf= np.log10(np.abs(np.fft.fft( x *  signal.hann(N_sample)))) * 20  # unit dB
        self.y_freq= np.fft.fftfreq(N_sample,1 / self.sr)
        
    def make_white_noise(self, sinpuku=0.05, length=1):
        self.xin=np.random.randn(self.sr * length ) * sinpuku
        
    def save_wav(self, yout, wav_path, wav_dir='wav_white_noise_in_out'):
        if not os.path.isdir(wav_dir):
            os.makedirs(wav_dir)
        out_path=os.path.join(wav_dir,wav_path)
        wavwrite( out_path, self.sr, ( yout * 2 ** 15).astype(np.int16))
        print ('save ', out_path) 

if __name__ == '__main__':
    
    # Length & Area value, from problems 3.8 in "Digital Processing of Speech Signals" by L.R.Rabiner and R.W.Schafer
    #
    # /a/
    L1_a=9.0    # set list of 1st tube's length by unit is [cm]
    A1_a=1.0    # set list of 1st tube's area by unit is [cm^2]
    L2_a=8.0    # set list of 2nd tube's length by unit is [cm]
    A2_a=7.0    # set list of 2nd tube's area by unit is [cm^2]
    
    # /u/
    L1_u=10.0   # set list of 1st tube's length by unit is [cm]
    A1_u=7.0    # set list of 1st tube's area by unit is [cm^2]
    L2_u=7.0    # set list of 2nd tube's length by unit is [cm]
    A2_u=3.0    # set list of 2nd tube's area by unit is [cm^2]
    
    # /o/: L3,A3 is  extend factor to /a/ connecting as /u/
    L3_o= L2_a * (L2_u / L1_u)     # set list of 3rd tube's length by unit is [cm]
    A3_o= A2_a * (A2_u / A1_u)     # set list of 3rd tube's area by unit is [cm^2]
    
    #
    from twotube import *
    from threetube import *
    from fourtube import *
    from fivetube import *
    
    
    # insatnce Two tube model example
    tube  =  Class_TwoTube(L1_a,L2_a,A1_a,A2_a, sampling_rate=48000*4)
    Class_Tube_Check(tube)
    
    
    # insatnce Three tube model example
    tube  =  Class_ThreeTube(L1_a, L2_a, L3_o, A1_a, A2_a, A3_o, sampling_rate=48000*4)
    Class_Tube_Check(tube)
    
    
    # insatnce Four tube model example
    tube  =  Class_FourTube(3., 8.8, 1.5, 1.4, 1.0, 8.4, 237., 11.6 , sampling_rate=48000*4)
    Class_Tube_Check(tube)
    
    
    # insatnce Five tube model example
    tube  =  Class_FiveTube(2.9, 8.7, 1.4, 2.9, 1.5, 1.0, 21.9, 438., 166., 62., sampling_rate=48000*4)
    Class_Tube_Check(tube)
    