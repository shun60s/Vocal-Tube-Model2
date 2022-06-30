#coding:utf-8

# Tube model:
# draw frequency response, cross-sectional view (area), and waveform, considering glottal voice source and mouth radiation
# save generated waveform as a wav file
#

import os
import argparse
import copy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from scipy.io.wavfile import write as wavwrite


# Check version
#  Python 3.10.4, 64bit on Win32 (Windows 10)
#  numpy 1.21.6
#  matplotlib  3.5.2
#  scipy 1.8.0


def show_figure1(tube, glo, hpf, yout=None):
    #
    fig = plt.figure()
    
    # draw frequency response
    ax1 = fig.add_subplot(311)
    amp0, freq=glo.H0(freq_low=100, freq_high=6000, Band_num=1024)
    amp1, freq=tube.H0(freq_low=100, freq_high=6000, Band_num=1024)
    amp2, freq=hpf.H0(freq_low=100, freq_high=6000, Band_num=1024)
    
    ax1.set_title('frequency response: ' )
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Amplitude [dB]')
    ax1.plot(freq, (amp0+amp1+amp2))
    
    plt.grid()
    
    # draw cross-sectional view
    ax2 = fig.add_subplot(312)
    L1=tube.L1
    L2=tube.L2
    L3=0.
    L4=0.
    L5=0.
    A1=tube.A1
    A2=tube.A2
    A3=0.
    A4=0.
    A5=0.
    if tube.num_of_tube >= 3:
        L3=tube.L3
        A3=tube.A3
    if tube.num_of_tube >= 4:
        L4=tube.L4
        A4=tube.A4
    if tube.num_of_tube >= 5:
        L5=tube.L5
        A5=tube.A5
    
    
    ax2.add_patch( patches.Rectangle((0, -0.5* A1), L1, A1, hatch='/', fill=False))
    ax2.add_patch( patches.Rectangle((L1, -0.5* A2), L2, A2, hatch='/', fill=False))
    ax2.add_patch( patches.Rectangle((L1+L2, -0.5* A3), L3, A3, hatch='/', fill=False))
    ax2.add_patch( patches.Rectangle((L1+L2+L3, -0.5* A4), L4, A4, hatch='/', fill=False))
    ax2.add_patch( patches.Rectangle((L1+L2+L3+L4, -0.5* A5), L5, A5, hatch='/', fill=False))
    ax2.set_xlim([0, L1+L2+L3+L4+L5+5])
    ax2.set_ylim([(max(A1,A2,A3,A4,A5)*0.5+5)*-1, max(A1,A2,A3,A4,A5)*0.5+5 ])
    
    ax2.set_title('cross-section area')
    plt.xlabel('Length [cm]')
    plt.ylabel('Cross-section area [ratio]')
    plt.grid()
    
    
    # draw generated waveform
    ax3 = fig.add_subplot(313)
    if yout is not None:
        ax3.set_title('generated waveform')
        plt.xlabel('msec')
        plt.ylabel('level')
        plt.plot( (np.arange(len(yout)) * 1000.0 / glo.sr) , yout)
        plt.grid()
    
    plt.tight_layout()
    plt.show()


def generate_waveform1(tube, glo, hpf, repeat_num=50):
    yg_repeat=glo.make_N_repeat(repeat_num=repeat_num)
    y2tm=tube.process(yg_repeat)
    yout=hpf.iir1(y2tm)
    return yout
    
def down_sample(xin, sampling_rate, over_sampling_ratio, cutoff=15000):
    if over_sampling_ratio == 1:
        return xin   # return xin itself, it's not over sample.
    
    # simple down sampler by FFT
    y= np.fft.fft(xin)
    freq= np.fft.fftfreq(len(xin),1 / sampling_rate)
    id0=np.where( freq > cutoff)[0][0]
    id1=len(xin) - id0
    y[id0:int(id1+1)]=0.
    z= np.real(np.fft.ifft(y))
    return  z.reshape(int(len(xin)/over_sampling_ratio),over_sampling_ratio,)[:,0].copy()

def save_wav(yout, wav_path, sampling_rate=48000, wav_dir='generated_waveform'):
    if not os.path.isdir(wav_dir):
        os.makedirs(wav_dir)
    out_path=os.path.join(wav_dir,wav_path)
    wavwrite( out_path, sampling_rate, ( yout * 2 ** 15).astype(np.int16))
    print ('save ', out_path) 


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='generate waveform using tube model')
    parser.add_argument('--osratio', '-r', type=int, default=4, help='specify over sampling ratio. if 1, not over sampling.')
    args = parser.parse_args()
    
    
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
    from glottal import *
    from HPF import *
    
    over_sampling_ratio= args.osratio
    
    
    # insatnce Two tube model example
    tube_2=  Class_TwoTube(L1_a,L2_a,A1_a,A2_a, sampling_rate=48000*over_sampling_ratio)
    
    # insatnce Three tube model example
    tube_3=  Class_ThreeTube(L1_a, L2_a, L3_o, A1_a, A2_a, A3_o, sampling_rate=48000*over_sampling_ratio)
    
    # insatnce Four tube model example
    tube_4p1=  Class_FourTube(3., 8.8, 1.5, 1.4, 1.0, 8.4, 237., 11.6 , sampling_rate=48000*over_sampling_ratio)
    
    # insatnce Five tube model examples
    tube_5p1= Class_FiveTube(2.9, 8.7, 1.4, 2.9, 1.5, 1.0, 21.9, 438., 166., 62., sampling_rate=48000*over_sampling_ratio)
    tube_5p2=Class_FiveTube(5.0, 2.6, 2.7, 2.6, 2.5, 19.1, 1.0, 1.3, 10.8, 69.0, sampling_rate=48000*over_sampling_ratio)
    tube_5p3= Class_FiveTube(1.3, 4.9, 1.1, 2.5, 4.9, 19.8, 2.8, 15.8, 1.0, 1.0, sampling_rate=48000*over_sampling_ratio)
    tube_5p4= Class_FiveTube(2.7, 5.4, 2.7, 1.4, 4.0, 1.0, 5.6, 2.1, 1.1, 3.2, sampling_rate=48000*over_sampling_ratio)
    tube_5p5= Class_FiveTube(3.3, 1.6, 5.9, 6.3, 1.3, 2.1, 3.0, 1.1, 8.8, 1.0, sampling_rate=48000*over_sampling_ratio)
    
    # choice one tube model to generate waveform
    tube=tube_5p1
    
    # specify output wav file name
    outFile="tube_5p1.wav"
    
    # instance as glottal voice source
    glo=Class_Glottal(sampling_rate=48000*over_sampling_ratio)
    
    # instance for mouth radiation effect
    hpf=Class_HPF(sampling_rate=48000*over_sampling_ratio)
    
    # generate waveform
    yout=generate_waveform1(tube, glo, hpf, repeat_num=50)
    
    # draw
    yout_early_portion= yout[0: int(len(glo.yg) * 3)]  # draw only 3 plus length early portion
    show_figure1(tube,glo,hpf,yout_early_portion)
    
    # save generated waveform as a wav file
    yout2=down_sample(yout, 48000*over_sampling_ratio, over_sampling_ratio)
    save_wav(yout2, outFile, sampling_rate=48000)