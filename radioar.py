#! python3
# this program is designed to plot the relation between extiction and distance in the selected region
import _pickle
from astropy.io import fits
from os.path import join as pjoin
from os.path import isfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
from scipy.ndimage.filters import gaussian_filter
import matplotlib
import numpy.ma as ma
import pandas as pd
import matplotlib.ticker as ticker
# Todo1 use the radio diagram as background


def dismod(dis):
    """convert distance to distance module"""
    mod=-5+5*np.log10(dis*1000)
    return mod


def moddis(mod):
    """convert distance module to distance"""
    dis=10**(0.2*mod+1)/1000
    return dis


def readfits(name, dirname = None, picdir = None):
    """this function can read the fits picture and return its xrange, yrange and imagedata
    Args:
        name: name of target SNR, eg.'snr169'
        dirname: a string representing the dir storing observation data
        picdir: a path of the dir storing observation data
    Returns:
        xranges: 1-d array, range of coordinates in x axis
        yranges: 1-d array, range of coordinates in y axis
        image_data: 2-d array, may contain nan in it; observation data
    """
    # radio pic path
    if dirname is None:
        dirname = 'NEW6CM'
    if picdir == None:
        image_file = pjoin('..', '..', 'Data', 'snrpic', '{0}'.format(dirname), '{0}.fits'.format(name))
    else:
        image_file = pjoin(picdir, '{0}.fits'.format(name))
    #read the fits pic and get the imagedata
    image_data=fits.getdata(image_file)
    #image_data=np.nan_to_num(image_data)
    if dirname == 'NEW6CM' or dirname == 'Effel11CM':
        ext = 0
        image_data=image_data[ext]
    elif dirname == '6cmSNR':
        ext = 1
    hdulist=fits.open(image_file)
    nx=hdulist[ext].header['NAXIS1']
    ny=hdulist[ext].header['NAXIS2']
    dx=hdulist[ext].header['CDELT1']
    dy=hdulist[ext].header['CDELT2']
    x0=hdulist[ext].header['CRVAL1']
    y0=hdulist[ext].header['CRVAL2']
    zx=hdulist[ext].header['CRPIX1']
    zy=hdulist[ext].header['CRPIX2']
    xranges=(np.arange(nx)-zx)*dx+x0
    yranges=(np.arange(ny)-zy)*dy+y0 
    image_data = ma.masked_invalid(image_data)
    return xranges, yranges, image_data

def readsnrsav(name, extsavdir = None):
    """read the extinction data from .sav files of Chen. return its xygrid ar data and disbins
    Args:
        name: name of target SNR, eg.'snr169'
        extsavdir：a path of the dir storing .sav files
    Returns:
        xgrid, ygrid: 1-d arrray, coordinates of centers of every grids
        ar: 1-d array with 2-d array as elements, extinction data loaded from .sav files
        realdis: 1-d array, distance of centers of each bin
    """
    #set the datapath
    if extsavdir == None:
        data_file = pjoin('..', '..', 'Data', 'extin3d', '{0}extin3d015.sav'.format(name))
    else:
        data_file = pjoin(extsavdir, '{0}extin3d015.sav'.format(name))
    #read cordinate, distance and extinction for each star from sav files
    res=readsav(data_file, python_dict=True)
    res=res['exta']
    xgrid=res.gl[0]
    ygrid=res.gb[0]
    ar=res.dar[0]
    dismo=res.dis[0]  #distant modulu
    realdis=10.**(dismo/5.+1)/1000.
    return xgrid, ygrid, ar, realdis

def readpoints(name, fig, num):
    """read coordinates of the clouds's border
    Args: 
        name: name of target SNR, eg.'snr169'
        fig: a figure object, on which we need to select border
        num: int, number of clouds
    Return:
        points: a list of coordinates 
    """
    pointsfile_path = pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_points_{1}.pkl'.format(name, num))
    if isfile(pointsfile_path):  #to see if points have been selected
        pointsfile = open(pointsfile_path, 'rb')  #if there is a file, load it
        points = _pickle.load(pointsfile)
        pointsfile.close()
    else: #get coordinates from function ginput
        points = fig.ginput(n=0, timeout = 0)
        #save the points
        pointsfile = open(pointsfile_path, 'wb')
        _pickle.dump(points, pointsfile, 2)
        pointsfile.close()
    return points




def radio_ar(name, ar, realdis, image_data, xr, yr, xgrid, ygrid, perpc, single, nradius, vmax, 
             vmin= None, fill = False, Number = None, level = None, lowlev = None, picn = None, cloudn = None):
    """draw contours of 3-d extinction map on diffirent distance ranges with radio picture as background, save as .pdf file.
    Args:
        name: name of target SNR, eg.'snr169'
        ar: 1-d array with 2-d array as elements, extinction data loaded from .sav files
        realdis: 1-d array, distance of centers of each bin
        image_data: 2-d array, may contain nan in it; observation data
        xr: 1-d array, range of coordinates in x axis
        yr: 1-d array, range of coordinates in y axis
        xgrid, ygrid: 1-d arrray, coordinates of centers of every grids
        perpc: Boolean value, if True , extinction data will be divided by distance of bin
        single: Boolean value, if False, contours will be ploted on all distance, else, only the distance selected
        nradius: float, the plotting size of the region
        vmin: float, lower limit of the intensity of radio picture, default is 0
        vmax: float, upper limit of the intensity of radio picture, unit is K
        fill : optional, Boolean value,  True, contourf; False, contour, Default is False
        Number: optional, (a, b), numbers of bins that need to be ploted, Default is (6,13)
        levels: optional, number of levels for contour, Default is 10
        lowlev: optional, levels could be hidden at head and tail, Default is 2
        picn: optional, number of the picture wanted, start from 0, Default is 4, the last ax will be deleted
        cloudn: optional, int, how many clouds we want to select, Default is 1
    Returns:
    """
    #默认参数
    if Number == None:
        Number = (6, 13)
    if level == None:
        level = 10
    if lowlev == None:
        lowlev = 2
    if picn == None:
        picn = 4
    if cloudn ==None:
        cloudn = 1
    if vmin == None:
        vmin = 0
    #将背景数据转化为DataFrame，并进行整理
    #ImageFrame = DataFrame(image_data, index=yr, columns=xr)
    #ImageFrame = ImageFrame.fillna(0)
    #ImageFrame.index.name = 'b'
    #ImageFrame.columns.name = 'l'
    #ImageFrame[ImageFrame<0] = 0.
    #ImageFrame = ImageFrame/1000
    imgdata = image_data/1000
    #读取SNR坐标和范围
    SNRs = pd.read_csv(pjoin('..', '..', 'Data', 'snrlist.csv'))
    SNRs.set_index('codename',inplace=True)
    xleft = SNRs.loc[name]['l']+SNRs.loc[name]['size']/60*nradius
    xright = SNRs.loc[name]['l']-SNRs.loc[name]['size']/60*nradius
    yup = SNRs.loc[name]['b']+SNRs.loc[name]['size']/60*nradius
    ydown = SNRs.loc[name]['b']-SNRs.loc[name]['size']/60*nradius

    #整理射电图数据
    xi = (xr<xleft) & (xr>xright)
    yi = (yr<yup) & (yr>ydown)
    imgdata = imgdata[yi,:]
    imgdata = imgdata[:,xi]
    extent = [xr[xi].max(), xr[xi].min(), yr[yi].min(), yr[yi].max()]
    
    #整理消光数据
    disbins = moddis((dismod(realdis) + 0.25)) - moddis((dismod(realdis) - 0.25))
    if perpc == True:
        for i in range(len(disbins)):
            ar.dar[i] /= disbins[i]
    adata = np.vstack(ar.dar)
    madata = ma.masked_where(adata<0, adata)
    allar = DataFrame(madata, index=[np.repeat(np.arange(len(ar.dar)),len(xgrid)), np.tile(ygrid, len(ar.dar))], columns=xgrid)
    allar.index.names = ['cell','gb']
    allar.columns.name = 'gl'
    xii = (xgrid<xleft) & (xgrid>xright) 
    yii = (ygrid<yup) & (ygrid>ydown)
    zmax = np.percentile(allar.loc[np.tile(yii,16), xii], 98)
    levs = zmax / level * np.arange(level+1)
    
    #绘图
    matplotlib.rcdefaults()
    if single == False:
        p = matplotlib.rcParams
        # 配置绘图区域的大小和位置，下面的值是基于图标的宽和高的比例
        p["figure.subplot.left"] = 0.05  # 左边距
        p["figure.subplot.right"] = 0.98   # 右边距
        p["figure.subplot.bottom"] = 0.1  # 下边距
        p["figure.subplot.top"] = 0.95   # 上边距
        # 配置subplots之间的间距（水平间距和垂直间距），也是基于图标的宽和高的比例
        p["figure.subplot.wspace"] = 0.05
        p["figure.subplot.hspace"] = 0.05
        fig, axes = plt.subplots(nrows=2, ncols=4, sharex='col', sharey='row', figsize=(7.75, 4.25))
        fig.delaxes(axes.flat[-1])
        for ii, ax in zip(range(Number[0], Number[1]), axes.flat[:-1]):
            #画射电灰度图
            ax.imshow(imgdata, cmap = 'gray_r', vmax = vmax, vmin = vmin, origin = 'lower', extent=extent)

            #画消光contour图
            iar = allar.loc[ii]
            iar = iar.loc[yii, xii]
            araa = gaussian_filter (iar ,0.68)
            if fill == False:
                cont = ax.contour(xgrid[xii], ygrid[yii], araa, levs[lowlev:-lowlev])
            else:
                cont = ax.contourf(xgrid[xii], ygrid[yii], araa, levs[lowlev:-lowlev])
            tit = r'{0:4.2f}~{1:4.2f}kpc'.format(realdis[ii] - disbins[ii]/2, realdis[ii] + disbins[ii]/2)
            ax.text(0.02, 0.99, tit, va='top', ha='left',transform=ax.transAxes, color='red', fontsize=10, zorder=10)

            #调整刻度
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            ax.tick_params(labelsize=9)
            
        #调整xylabel
        fig.text(0.5, 0., 'Galactic Latitude', ha='center', va='bottom', fontsize=14)
        fig.text(0., 0.5, 'Galactic Longitude', ha='left', va='center', rotation='vertical', fontsize=14)

        print('contour levels:',levs[lowlev:-lowlev],'mag/kpc')
        fig.savefig(pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_all.pdf'.format(name)), transparent = True, dpi=360)
    else:
        ii = Number[0]+picn
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111)
        ax.imshow(imgdata, cmap = 'gray_r', vmax = vmax, vmin = vmin, origin = 'lower', extent = extent)
        iar = allar.loc[ii]
        iar = iar.loc[yii, xii]
        araa = gaussian_filter (iar ,0.68)
        if fill == False:
            cont = ax.contour(xgrid[xii], ygrid[yii], araa, levs[lowlev:-lowlev])
        else:
            cont = ax.contourf(xgrid[xii], ygrid[yii], araa, levs[lowlev:-lowlev])
        tit = r'{0:4.2f}~{1:4.2f}kpc'.format(realdis[ii] - disbins[ii]/2, realdis[ii] + disbins[ii]/2)
        ax.text(0.02, 0.99, tit, va='top', ha='left',transform=ax.transAxes, color='red', fontsize=12, zorder=10)
        for k in range(cloudn):
            points = readpoints(name, fig, k+1)
            points = np.array(points)
            x = points[:, 0]
            y = points[:, 1]
            x = np.append(x, x[0])  # to join the start point and end point on the picture
            y = np.append(y, y[0])
            ax.plot(x, y, 'r')

        #调整刻度
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.tick_params(labelsize=12)
        
        #调整xylabel
        fig.text(0.5, 0., 'Galactic Latitude', ha='center', va='bottom', fontsize=12)
        fig.text(0., 0.5, 'Galactic Longitude', ha='left', va='center', rotation='vertical', fontsize=12)
        fig.savefig(pjoin('..', '..', 'Data', 'extin3d', 'results', '{0}_single.pdf'.format(name)), transparent = True, dpi=360)
    plt.show()
if __name__=='__main__':
    name = 'snr169'
    perctl = 98
    perpc=True
    single = False
    xr, yr, image_data = readfits(name, ext = 0)
    xgrid, ygrid, ar, realdis = readsnrsav(name)
    radio_ar(name, ar, realdis, image_data, xr, yr, xgrid, ygrid, perctl, perpc, single)