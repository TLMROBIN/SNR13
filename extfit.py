#! python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
from scipy.special import erf
import matplotlib.ticker as ticker
import sys
from astropy.io.fits import getdata
import _pickle
from SNR13.dis_ext import binMedMean
from scipy.optimize import curve_fit
import matplotlib


def fitmodel(x, a, b, d0, d_ar):
    model = a * x + b * x ** 2 + d_ar / 2.0 * (1 + erf((x - d0) / np.sqrt(2) / (45/60./360*2*np.pi)/d0))
    return model


def exfit(ax, snrNum, cloudNum):
    
    pointsfile_path = '../../data/map2018/borders/snr{0}_border_{1}.pkl'.format(snrNum, cloudNum)
    pointsfile = open(pointsfile_path, 'rb')  # if there is a file, load it
    points = _pickle.load(pointsfile)
    pointsfile.close()
    p = path.Path(points)  # form a path object , the border of the region selected
    stardata = getdata('../../data/map2018/snrs/E3/snr{0}new.fits'.format(snrNum))
    coordinates = np.vstack([stardata.L, stardata.B]).T
    EGK = stardata.EGK
    dis = stardata.D / 1000
    eEGK = 0.5 * (stardata.EGKHI - stardata.EGKLO)
    edis = 0.5 * (stardata.DHI - stardata.DLO) / 1000
    index = p.contains_points(coordinates) & (EGK >0)
    # extract data 
    EGK = EGK[index]
    dis = dis[index]
    eEGK = eEGK[index]
    edis = edis[index]
    N = 300
    paras = np.zeros((N,4))
    for i in range(N):
        disr = np.random.normal(dis, edis)
        EGKr = np.random.normal(EGK, eEGK)
        mean, med = binMedMean([disr, EGKr], 0.2, 3.2, 0.1, 3, 3)
        popt, pcov = curve_fit(f=fitmodel, xdata=med[0], ydata=med[1],
                            p0=[0.81, -0.07, 1.0, 0.2],
                            bounds=([0, -10, 0.2, 0], [10, 10, 1.5, 10]))
        paras[i] = popt
    results = np.percentile(paras, [16, 50, 84], axis=0)
    pfit = results[1]
    perror = (results[2], results[0])
    # plotting the result
    matplotlib.rcdefaults()
    p = matplotlib.rcParams
    p["figure.subplot.left"] = 0.125
    p["figure.subplot.right"] = 0.95   
    p["figure.subplot.bottom"] = 0.12  
    p["figure.subplot.top"] = 0.95   
    p["figure.subplot.wspace"] = 0.05
    p["figure.subplot.hspace"] = 0.05
    plt.rc('text', usetex=True)
    stars = ax.plot(dis, EGK, '.k', ms=3, label='Sources in selected area')
    x = np.linspace(0, 3.5)
    for para in paras[np.random.randint(len(paras), size=100)]:
        y = fitmodel(x, *para)
        _ = ax.plot(x, y, color="b", alpha=0.1)
    mean, med = binMedMean([dis, EGK], 0.2, 3.2, 0.1, 3, 3)
    medValue = ax.plot(med[0], med[1], 'r.', ms=8, label='Medians')
    a, b = pfit[0], pfit[1]
    fakeline = ax.plot(x, a * x + b * x ** 2, 'r--', label='EGK(d) without MC')
    result = r'$a$, $b$, $d_0$, $\delta E$ = {0:4.2f}, {1:4.3f}, {2:4.2f}, {3:4.2f}'.format(
        pfit[0], pfit[1], pfit[2], pfit[3])
    ax.text(0.02, 0.9, result, va='top', ha='left', transform=ax.transAxes, color='red', fontsize=12, zorder=10)
    ax.set_xlim([0, 3.5])
    ax.set_ylim([0, 2.5])
    ax.set_ylabel(r'$\rm E(G-K_s)$', fontsize=12)
    ax.set_xlabel('Dis/kpc', fontsize=12)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.tick_params(labelsize=10)
    uplimit = perror[0] - pfit
    lowlimit = pfit - perror[1]
    print('a, b, d0, dE:{0:4.2f}, {1:4.3f}, {2:4.3f}, {3:4.2f}'.format(
        pfit[0], pfit[1], pfit[2], pfit[3]))
    print('upperlimit:{0:4.2f}, {1:4.3f}, {2:4.3f}, {3:4.2f}'.format(
        uplimit[0], uplimit[1], uplimit[2], uplimit[3]))
    print('lowerlimit:{0:4.2f}, {1:4.3f}, {2:4.3f}, {3:4.2f}'.format(
        lowlimit[0], lowlimit[1], lowlimit[2], lowlimit[3]))

def fitmodel2(x, a, b, d01, d_ar1, d02, d_ar2):
    model = a * x + b * x ** 2 + \
    d_ar1 / 2.0 * (1 + erf((x - d01) / np.sqrt(2) / (45/60./360*2*np.pi)/d01)) + \
    d_ar2 / 2.0 * (1 + erf((x - d02) / np.sqrt(2) / (45/60./360*2*np.pi)/d02))
    return model


def exfit2(ax, snrNum, cloudNum):
    
    pointsfile_path = '../../data/map2018/borders/snr{0}_border_{1}.pkl'.format(snrNum, cloudNum)
    pointsfile = open(pointsfile_path, 'rb')  # if there is a file, load it
    points = _pickle.load(pointsfile)
    pointsfile.close()
    p = path.Path(points)  # form a path object , the border of the region selected
    stardata = getdata('../../data/map2018/snrs/E3/snr{0}new.fits'.format(snrNum))
    coordinates = np.vstack([stardata.L, stardata.B]).T
    EGK = stardata.EGK
    dis = stardata.D / 1000
    eEGK = 0.5 * (stardata.EGKHI - stardata.EGKLO)
    edis = 0.5 * (stardata.DHI - stardata.DLO) / 1000
    index = p.contains_points(coordinates) & (EGK >0)
    # extract data 
    EGK = EGK[index]
    dis = dis[index]
    eEGK = eEGK[index]
    edis = edis[index]
    N = 300
    paras = np.zeros((N,6))
    for i in range(N):
        disr = np.random.normal(dis, edis)
        EGKr = np.random.normal(EGK, eEGK)
        mean, med = binMedMean([disr, EGKr], 0.2, 3.2, 0.1, 3, 3)
        popt, pcov = curve_fit(f=fitmodel2, xdata=med[0], ydata=med[1],
                            p0=[0.81, -0.07, 1.0, 0.2, 1.5, 0.2],
                            bounds=([0, -10, 0.5, 0, 1.0, 0], [10, 10, 1.5, 10, 2.5, 10]))
        paras[i] = popt
    results = np.percentile(paras, [16, 50, 84], axis=0)
    pfit = results[1]
    perror = (results[2], results[0])
    # plotting the result
    matplotlib.rcdefaults()
    p = matplotlib.rcParams
    p["figure.subplot.left"] = 0.125
    p["figure.subplot.right"] = 0.95   
    p["figure.subplot.bottom"] = 0.12  
    p["figure.subplot.top"] = 0.95   
    p["figure.subplot.wspace"] = 0.05
    p["figure.subplot.hspace"] = 0.05
    plt.rc('text', usetex=True)
    stars = ax.plot(dis, EGK, '.k', ms=3, label='Sources in selected area')
    x = np.linspace(0, 3.5)
    for para in paras[np.random.randint(len(paras), size=100)]:
        y = fitmodel2(x, *para)
        _ = ax.plot(x, y, color="b", alpha=0.1)
    mean, med = binMedMean([dis, EGK], 0.2, 3.2, 0.1, 3, 3)
    medValue = ax.plot(med[0], med[1], 'r.', ms=8, label='Medians')
    a, b = pfit[0], pfit[1]
    fakeline = ax.plot(x, a * x + b * x ** 2, 'r--', label='EGK(d) without MC')
    result = r'$a$, $b$, $d_0$, $\delta E$ = {0:4.2f}, {1:4.3f}, {2:4.2f}, {3:4.2f}, {4:4.2f}, {5:4.2f}'.format(
        pfit[0], pfit[1], pfit[2], pfit[3], pfit[4], pfit[5])
    ax.text(0.02, 0.9, result, va='top', ha='left', transform=ax.transAxes, color='red', fontsize=12, zorder=10)
    ax.set_xlim([0, 3.5])
    ax.set_ylim([0, 2.5])
    ax.set_ylabel(r'$\rm E(G-K_s)$', fontsize=12)
    ax.set_xlabel('Dis/kpc', fontsize=12)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.tick_params(labelsize=10)
    uplimit = perror[0] - pfit
    lowlimit = pfit - perror[1]
    print('a, b, d0, dE:{0:4.2f}, {1:4.3f}, {2:4.3f}, {3:4.2f}, {4:4.3f}, {5:4.2f}'.format(
        pfit[0], pfit[1], pfit[2], pfit[3], pfit[4], pfit[5]))
    print('upperlimit:{0:4.2f}, {1:4.3f}, {2:4.3f}, {3:4.2f}, {4:4.3f}, {5:4.2f}'.format(
        uplimit[0], uplimit[1], uplimit[2], uplimit[3], uplimit[4], uplimit[5]))
    print('lowerlimit:{0:4.2f}, {1:4.3f}, {2:4.3f}, {3:4.2f}, {4:4.3f}, {5:4.2f}'.format(
        lowlimit[0], lowlimit[1], lowlimit[2], lowlimit[3], lowlimit[4], lowlimit[5]))


if __name__ == '__main__':
    snrNum = sys.argv[1]
    cloudNum = sys.argv[2]
    switch = int(sys.argv[3])
    fig = plt.figure(figsize=(4.5, 4))
    ax = fig.add_subplot(111)
    if switch == 1:
        exfit(ax, snrNum, cloudNum)
    else:
        exfit2(ax, snrNum, cloudNum)

plt.show()