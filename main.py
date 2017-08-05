#! python3
from SNR13 import *

def step1_2(name, perctl, perpc, single, picn = None, cloudn = None):
    '''plotting the extinction per kpc on all distance'''
    xr, yr, image_data =  radioar.readfits(name, ext = 0)
    xgrid, ygrid, ar, realdis = radioar.readsnrsav(name)
    radioar.radio_ar(name, ar, realdis, image_data, xr, yr, xgrid, ygrid, perctl, perpc, single = single, picn = picn, cloudn = cloudn)

def step3(name, ncloud):
    for i in range(ncloud):
        dis_ext.datainborder(name, cloudn = i+1)
        pfit, perror = dis_ext.disarMCMCall(name, strategy='MC', cloudn = i+1)
        print('Fitted paras for cloud {0}: '.format(i+1), pfit)
        print('Errors for cloud {0}: '.format(i+1), perror)
if __name__=='__main__':
    name = 'snr169'
    perctl = 98
    perpc=True
    single = False
    picn = 4
    cloudn = 1

    step1(name, perctl, perpc)
    step2(name, perctl, perpc, picn, cloudn)
    step3(name, cloudn)