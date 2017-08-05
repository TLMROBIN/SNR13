#! python3
# this program is designed to plot the relation between extiction and distance in the selected region
import _pickle
#import corner
#import re
import emcee
from os.path import join as pjoin
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
from matplotlib import path
from scipy.special import erf
from scipy import stats
import numpy.ma as ma
from pandas import  Series, DataFrame
import pandas as pd
import matplotlib.ticker as ticker


#Todo4 choose the stars in that region
def inborder(whichborder, points = None):
#to judge if points are located in the cloud border we choose
    #read the border coordinates
    #bordrname = re.compile(r'snr\d\d\d')
    #mat = bordrname.search(whichborder)
    name = whichborder
    if points == None:
        pointsfile_path = pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_points.pkl'.format(name))
        pointsfile = open(pointsfile_path, 'rb')#if there is a file, load it
        points = _pickle.load(pointsfile)
        pointsfile.close()
    p = path.Path(points)#form a path object, which represent the border
    #judge if  stars are inside   the border
    #get coordinates imformation
    stardata = readsav(pjoin('..', '..', 'Data', 'extin3d', 'discalib', 'result', '{0}.sav'.format(name)), python_dict = True)
    stardata = stardata['res']
    coordinates = np.vstack(stardata.rd)[:,0:2]
    index = p.contains_points(coordinates)
    return index
def exdis_inborder(index, name, giant ='True'):
#get ar and dis of those stars in the border and may delete giants from them
    stardata = readsav(pjoin('..', '..', 'Data', 'extin3d', 'discalib', 'result', '{0}.sav'.format(name)), python_dict = True)
    stardata = stardata['res']
    data = stardata[index]
    data = data[data.ar > 0] #delete the data where ar = 0
    index2 = np.logical_or(data.dis > 1, data.ar < 0.5 + 1.5 * data.dis)
    if giant == 'False':
        data = data[index2]
    return data
def disarinborder(name, imgpath, img2 = False):
    borderpoints = readpointsx(imgpath, img2 = img2)#read the coordinates of clouds's border
    index = inborder(name, borderpoints)#select stars in choosen border,return their index
    stardata = exdis_inborder(index, name, giant = 'True')#ar,dis in the border
    if img2== False:
        datapath = pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_diar.pkl'.format(name))
    else:
        datapath = pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_diar_2.pkl'.format(name))
    datafile = open(datapath, 'wb')
    _pickle.dump(stardata, datafile, 2)
    datafile.close()
def binArray(data, axis, start, end, length, func=np.nanmean):
# bin the data in given axis with start, end and binsize, apply the func on bined data
    #transform the data into a uniform shape for further manipulate
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis]= argdims[axis], argdims[0]
    data = data.transpose(argdims)
    bin_number = (end - start) // length#decide how many bins will be choosen
    # use for loop in a list, get the index of selected data, take them with np.take method in given axis, then apply the given function on them
    data = [func(np.take(data,np.arange(len(data))[np.logical_and(data[:,0]>start + i * length, data[:,0]<start + (i + 1) * length)],0),0) for i in np.arange(bin_number)] 
    data = np.array(data).transpose(argdims)#tranform back
    return data
def medifit(dis, ar, start, end, step, delstd = False):
    data = np.array([dis,ar])
    #data_nogiant=np.array([dis_nogiant, ar_nogiant])
    if delstd == True:
        data = bin3std(data, axis, start, end, length)
    medi = binArray(data, 1, start, end, step, np.median)
    medi_sem = binArray(data, 1, start, end, step, stats.sem)
    #medi_nogiant=binArray(data_nogiant,1,0.5,4.0,step,np.median)
    return medi, medi_sem#, medi_nogiant
def sigmaclips(data, start, end, length, sigmanum):
    #sigmaclip the samples in bins and return their mean and error
    data=np.array(data)
    bin_number = int((end - start) // length)#decide how many bins will be choosen
    binmed = np.zeros(bin_number)
    binmed_err = np.zeros(bin_number)
    bindis = np.zeros(bin_number)
    arlist = []
    dislist = []
    for i in np.arange(bin_number):
        index = np.logical_and(data[0] > start + i * length, data[0] < start + (i + 1) * length)
        if np.count_nonzero(index) >= 1:
            binar = data[1][index]
            bindis = data[0][index]
            #arb = sigma_clip(binar, sigma = sigmanum, cenfunc = np.mean)
            #disb = np.ma.masked_array(bindis, arb.mask)
            c, l, u= stats.sigmaclip(binar, sigmanum, sigmanum)
            ix = np.where((l <= binar) & (binar <= u))[0]
            arb = binar[ix]
            disb = bindis[ix]
        else:
            arb = []
            disb = []
        arlist.append(arb)
        dislist.append(disb)
    ar = np.concatenate(arlist)
    dis = np.concatenate(dislist)
    ar = np.ma.compressed(ar)
    dis = np.ma.compressed(dis)
    return dis, ar

def tryall(name,Number = (6,13), perctl = 98, perpc=True, fill = False, level = 10, lowlev=2):
    xr, yr, image_data = readfits(name, ext=0)
    xgrid, ygrid, ar, realdis = readsnrsav(name)
    radioar(name, ar, realdis, image_data, xr, yr, xgrid, ygrid, perctl= perctl, perpc = perpc, fill = fill, Number = Number, level = level, lowlev=lowlev)

def trysingle(name,Number=(10,11), img2 = False, perctl = 98, perpc=True, fill = False, level = 10, lowlev=2):
    xr, yr, image_data = readfits(name, ext=0)
    xgrid, ygrid, ar, realdis = readsnrsav(name)
    imgpath=radioar(name, ar, realdis, image_data, xr, yr, xgrid, ygrid, perctl= perctl, perpc=perpc, fill = fill, Number = Number, level = level, lowlev=lowlev)
    if img2 == False:
        disarinborder(name, imgpath, img2 = False)
    elif img2 ==True:
        disarinborder(name, imgpath, img2 = False)
        disarinborder(name, imgpath, img2 = True)

def lnlike(theta, x, y, yerr, D):
    a, b, d0, d_ar= theta
    model = a * x + b * x ** 2 + d_ar / 2.0 * (1 + erf((x - d0) / np.sqrt(2) / (D/60./360*2*np.pi)/d0))
    return np.sum(np.log(1 / (np.sqrt(2 * np.pi) * yerr)) + (-(y - model) ** 2 / (2 * yerr ** 2)))

def lnprior(theta):
    a, b, d0, d_ar= theta
    if -10 < a < 10 and -1 < b < 1 and 0 < d0 < 2 and 0 < d_ar < 3. :
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr, D):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr, D)

def MCMCbasic(dis, ar, ndim, nwalkers, start, end, step, p0, diameter, sigmanum):
    ##
    x, y = sigmaclips(np.array([dis,ar]), start, end, step, sigmanum)
    #x, y, yerr = np.insert(x, 0, 0), np.insert(y, 0, 0), np.insert(yerr, 0, 0.01)
    yerr = np.zeros_like(x)
    D = np.zeros_like(x)
    yerr[:] = 0.15
    D[:] = diameter
    pos = [np.array(p0) + 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr, D))
    sampler.run_mcmc(pos, 2000)
    samples = sampler.chain[:, 200:, :].reshape((-1, ndim))
    #samfig = corner.corner(samples, labels = ["$a$", "$b$", r'$\d_0\$', r'$\deltaA_r\$' ])
    results = np.percentile(samples, [16, 50, 84], axis = 0)
    pfit = results[1]
    lower = results[0] - results[1]
    uper = results[2] - results[1]
    return pfit, samples, lower, uper

def disarMCMCall(name, ndim = 4, nwalkers = 100, start = 0.3, end = 3.0 ,step = 0.1, p0 = [0.81, -0.07, 1.2, 0.2], 
                 diameter = 45, img2 = False, sigmanum = 2, strategy = None, N = 100):
    if img2 == False:
        datapath = pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_diar.pkl'.format(name))
    else:
        datapath = pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_diar_2.pkl'.format(name))
    datafile = open(datapath, 'rb')
    data = _pickle.load(datafile)
    dis, ar = data.dis, data.ar
    datafile.close()
    if strategy == 'Bootstrap':
        paras = np.zeros((N,4))
        for i in range(N):
            ix = np.random.choice(len(dis), int(len(dis)*0.9))
            disi = dis[ix]
            ari = ar[ix]
            paras[i], _, _, _= MCMCbasic(disi, ari, ndim, nwalkers, start, end, step, p0, diameter, sigmanum)
        pfit = np.mean(paras, axis = 0)
        perror = np.std(paras, axis = 0)
    elif strategy == 'MC':
        paras = np.zeros((N,4))
        for i in range(N):
            disi = dis * (1+np.random.normal(loc=0.,scale=0.2, size=len(dis)))
            ari = ar
            paras[i], _, _, _= MCMCbasic(disi, ari, ndim, nwalkers, start, end, step, p0, diameter, sigmanum)
        pfit = np.mean(paras, axis = 0)
        perror = np.std(paras, axis = 0)
    else:
        pfit, paras, l, u = MCMCbasic(dis, ar, ndim, nwalkers, start, end, step, p0, diameter, sigmanum)
        perror = (u - l) / 2
    
    x1 = np.arange(0,5,0.01)
    fig = plt.figure()
    plt.rc('text', usetex=True)
    ax = fig.add_subplot(111)
    stars = ax.plot(dis,ar,'.k', ms = 3, label = r'$Sources~in~selected~area$')
    for a, b, d0, d_ar  in paras[np.random.randint(len(paras), size=100)]:
        _=ax.plot(x1, a * x1 + b * x1 ** 2 + d_ar / 2.0 * (1 + erf((x1 - d0) / np.sqrt(2) / (diameter/60./360*2*np.pi)/d0)) ,
                  color="b", alpha=0.1)
    #star_sem = ax.errorbar(x, y, yerr = yerr, fmt = '.', color = 'r', capsize = 2, zorder=10)
    #star_sem.set_label('Medians and errors')
    ax.set_xlim([0,3.5])
    ax.set_ylim([0,4.0])
    xla = ax.set_xlabel(r"$\rm d~(kpc)$", fontsize=14)
    yla = ax.set_ylabel(r"$\rm A_r~(mag)$", fontsize=16)
    a, b = pfit[0], pfit[1]
    fakeline = ax.plot(x1, a * x1 + b * x1 ** 2, 'r--', label = r'$A_r(d)~without~MC$')
    ax.legend()
    #tit = ax.set_title('extinction vs distance of {0}'.format(name))
    if img2 == False:
        eximg = pjoin('..', '..', 'Data', 'extin3d', 'results', 'snrext_{0}.png'.format(name))
        samimg = pjoin('..', '..', 'Data', 'extin3d', 'results', 'corner_{0}.png'.format(name))
    else:
        eximg = pjoin('..', '..', 'Data', 'extin3d', 'results', 'snrext_{0}_2.png'.format(name))
        samimg = pjoin('..', '..', 'Data', 'extin3d', 'results', 'corner_{0}_2.png'.format(name))
    fig.savefig(eximg, dpi = 1000, transparent = True)
    #samfig.savefig(samimg)
    plt.show()
    return pfit, perror

if __name__=='__main__':
    name = 'snr169'
    tryall(name)
    trysingle(name)
    disarMCMC(name)