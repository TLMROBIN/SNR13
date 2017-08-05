#! python3
from SNR13 import *

def step1_2(name, ext, perctl, perpc, single, picn = None, cloudn = None, dirname = None):
    '''plotting the extinction per kpc on all distance
    Args:
        dirname: a string representing the dir storing observation data, passed to  readfits
    '''
    xr, yr, image_data =  radioar.readfits(name, dirname = dirname)
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

    step1_2(name, perctl, perpc)
    step1_2(name, perctl, perpc, single, picn, cloudn)
    step3(name, cloudn)