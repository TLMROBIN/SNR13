from astropy.io import fits
from scipy.io import readsav
from scipy.ndimage.filters import gaussian_filter
import matplotlib.ticker as ticker
import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.io.fits import getdata
from os.path import isfile
import _pickle
from matplotlib import path


def readpoints(name, fig, num):
    """read coordinates of the clouds's border
    Args: 
        name: name of target SNR, eg.'snr169'
        fig: a figure object, on which we need to select border
        num: int, number of clouds
    Return:
        points: a list of coordinates 
    """
    pointsfile_path = '../../data/map2018/borders/snr{0}_border_{1}.pkl'.format(name, num)
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


def datainborder(name, cloudn, points):
    """to judge which points are located in the cloud border we choose, then get and save their data
    Args:
        name: name of target SNR, eg.'snr169'
        points: a list of number pairs, if None, points will be loaded from file
        cloudn : number of clouds, if None, only one cloud's resources will be selected
    Return:
    """
    #judge which stars are located in the region
    p = path.Path(points)  #form a path object, which represent the border of the region selected 
    stardata = fits.getdata('../../data/map2018/snrs/snr{0}new.fits'.format(name))
    coordinates = np.vstack([stardata.l, stardata.b]).T
    AG = stardata.REGK + 1.987 * stardata.REHK
    dis = stardata.DIS
    index = p.contains_points(coordinates) & (AG >0)
    # extract data 
    AG = AG[index]
    dis = dis[index]
    #save the data, so time can be saved for the next time
    datapath = '../../data/map2018/snrsext/{0}_disar_{1}.pkl'.format(name, cloudn)
    datafile = open(datapath, 'wb')
    _pickle.dump([dis, AG], datafile, 2)
    datafile.close()


def dgplot1(ax1, snrNum, level, lC, bC, hl, num):
    # 只画消光数据
    # 整理消光数据
    res = readsav('../../data/map2018/IDLs/extin3d015{0}.sav'.format(snrNum))
    midMum = res.exta.dis[0] # 距离间隔的中间值，距离模数
    gl = res.exta.gl[0]  # 网点的银经
    gb = res.exta.gb[0]  # 网点的银维
    # debr = np.stack(res.exta.debr[0].dar)  # 区间内的E(B-V)h值，三维数组，23层，每层是一个二维矩阵 
    degk = np.stack(res.exta.degk[0].dar)
    AG = degk    # + 1.987 * dehk  # G波段消光
    disBins = [[i - 0.125, i + 0.125] for i in midMum]
    AG_new = np.zeros([int(AG.shape[0] / 2), AG.shape[1], AG.shape[2]])
    for i in range(int(AG.shape[0] / 2)):
        AG_new[i, :, :] = np.sum(AG[2 * i:2 * i +2, :, :], axis=0)
    AG = AG_new
    # 绘制消光/距离等值线图
    zmax = np.percentile(AG, 98)
    levs = zmax / level * (np.arange(level) + 1)[1:]
    extin = gaussian_filter(AG[num, :, :], 0.68)
    cont = ax1.contour(gl, gb, extin, levs)
    ax1.set_xlim(lC + hl, lC - hl)
    ax1.set_ylim(bC - hl, bC + hl)
    tit = r'{0:4.2f}-{1:4.2f}'.format(disBins[num][0] * 2, disBins[num][1] * 2)
    ax1.text(0.02, 1.05, tit, va='top', ha='left', transform=ax1.transAxes, color='red', fontsize=12, zorder=10)
    # 调整刻度
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax1.tick_params(labelsize=10)
    # 调整xylabel
    ax1.set_xlabel('Galactic Latitude', fontsize=12)
    ax1.set_ylabel('Galactic Longitude', fontsize=12)

    

def dgplot2(fig, ax2, snrNum, level, lC, bC, hl, num, cloudNum, colorMethod, vmax):
    # 在射电图上画消光数据，并圈出云的范围
    # 整理射电图数据
    image_file = '../../data/snrpic/6cmSNR/snr{0}.fits'.format(snrNum)
    image_data, hdr = getdata(image_file, header=True)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    dx = hdr['CDELT1']
    dy = hdr['CDELT2']
    x0 = hdr['CRVAL1']
    y0 = hdr['CRVAL2']
    zx = hdr['CRPIX1']
    zy = hdr['CRPIX2']
    xranges = (np.arange(nx)-zx)*dx+x0
    yranges = (np.arange(ny)-zy)*dy+y0 
    # image_data = ma.masked_invalid(image_data)
    imadata = np.zeros_like(image_data)
    imadata[image_data > 0] = np.log(image_data[image_data > 0])
    extent = [xranges.max() - 0.5 * dx, xranges.min() + 0.5 * dx,
            yranges.min() - 0.5 * dy, yranges.max() + 0.5 * dy]
    # image_data = ma.masked_less_equal(image_data,0)
    # 整理消光数据
    res = readsav('../../data/map2018/IDLs/extin3d015{0}.sav'.format(snrNum))
    midMum = res.exta.dis[0] # 距离间隔的中间值，距离模数
    gl = res.exta.gl[0]  # 网点的银经
    gb = res.exta.gb[0]  # 网点的银维
    # debr = np.stack(res.exta.debr[0].dar)  # 区间内的E(B-V)h值，三维数组，23层，每层是一个二维矩阵 
    degk = np.stack(res.exta.degk[0].dar)
    dehk = np.stack(res.exta.dehk[0].dar)
    AG = degk    # + 1.987 * dehk  # G波段消光
    disBins = [[i - 0.125, i + 0.125] for i in midMum]
    AG_new = np.zeros([int(AG.shape[0] / 2), AG.shape[1], AG.shape[2]])
    for i in range(int(AG.shape[0] / 2)):
        AG_new[i, :, :] = np.sum(AG[2 * i:2 * i +2, :, :], axis=0)
    AG = AG_new

    # 绘制消光/距离等值线图
    zmax = np.percentile(AG, 98)
    levs = zmax / level * (np.arange(level) + 1)[1:]
    extin = gaussian_filter(AG[num, :, :], 0.68)

    # 画叠加图
    if colorMethod == 0:
        im = ax2.imshow(imadata, origin='lower', interpolation='nearest', vmin=0, vmax=np.sqrt(vmax), cmap='jet', extent=extent)
    elif colorMethod == 1:
        im = ax2.imshow(imadata, origin='lower', interpolation='nearest', vmin=0, vmax=vmax, cmap='jet', extent=extent)
    elif colorMethod == 2:
        im = ax2.imshow(imadata, origin='lower', interpolation='nearest', vmin=0, vmax=np.log(vmax), cmap='jet', extent=extent)
    ax2.set_xlim(lC + hl, lC - hl)
    ax2.set_ylim(bC - hl, bC + hl)
    #plt.colorbar(im, ax=ax2, extend='max', fraction=0.05, pad=0.01, panchor=(2.0, 0.5))
    extin = gaussian_filter(AG[num, :, :], 0.68)
    cont = ax2.contour(gl, gb, extin, levs, cmap='gray')
    tit = r'{0:4.2f}-{1:4.2f}'.format(disBins[num][0] * 2, disBins[num][1] * 2)
    ax2.text(0.02, 1.05, tit, va='top', ha='left', transform=ax2.transAxes, color='red', fontsize=12, zorder=10)
    # 调整刻度
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax2.tick_params(labelsize=10)
    # 调整xylabel
    ax2.set_xlabel('Galactic Latitude', fontsize=12)
    ax2.set_ylabel('Galactic Longitude', fontsize=12)
    # 读取边界点
    points = readpoints(snrNum, fig, cloudNum)
    datainborder(snrNum, cloudNum, points)
    points = np.array(points)
    x = points[:, 0]
    y = points[:, 1]
    x = np.append(x, x[0])  # to join the start point and end point on the picture
    y = np.append(y, y[0])
    ax2.plot(x, y, 'r', lw=2)

def dgplot3(ax2, snrNum, lC, bC, hl, colorMethod, vmax):
    # 只画射电图像
    # 整理射电图数据
    image_file = '../../data/snrpic/6cmSNR/snr{0}.fits'.format(snrNum)
    image_data, hdr = getdata(image_file, header=True)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    dx = hdr['CDELT1']
    dy = hdr['CDELT2']
    x0 = hdr['CRVAL1']
    y0 = hdr['CRVAL2']
    zx = hdr['CRPIX1']
    zy = hdr['CRPIX2']
    xranges = (np.arange(nx)-zx)*dx+x0
    yranges = (np.arange(ny)-zy)*dy+y0 
    # image_data = ma.masked_invalid(image_data)
    imadata = np.zeros_like(image_data)
    imadata[image_data > 0] = np.log(image_data[image_data > 0])
    extent = [xranges.max() - 0.5 * dx, xranges.min() + 0.5 * dx,
            yranges.min() - 0.5 * dy, yranges.max() + 0.5 * dy]
    # image_data = ma.masked_less_equal(image_data,0)
    # 画叠加图
    if colorMethod == 0:
        im = ax2.imshow(imadata, origin='lower', interpolation='nearest', vmin=0, vmax=np.sqrt(vmax), cmap='jet', extent=extent)
    elif colorMethod == 1:
        im = ax2.imshow(imadata, origin='lower', interpolation='nearest', vmin=0, vmax=vmax, cmap='jet', extent=extent)
    elif colorMethod == 2:
        im = ax2.imshow(imadata, origin='lower', interpolation='nearest', vmin=0, vmax=np.log(vmax), cmap='jet', extent=extent)
    ax2.set_xlim(lC + hl, lC - hl)
    ax2.set_ylim(bC - hl, bC + hl)
    #plt.colorbar(im, ax=ax2, extend='max', fraction=0.05, pad=0.01, panchor=(2.0, 0.5))
    # 调整刻度
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax2.tick_params(labelsize=10)
    # 调整xylabel
    ax2.set_xlabel('Galactic Latitude', fontsize=12)
    ax2.set_ylabel('Galactic Longitude', fontsize=12)


if __name__ == '__main__':
    snrNum = sys.argv[1]
    level = int(sys.argv[2])
    lC = float(sys.argv[3])
    bC = float(sys.argv[4])
    hl = float(sys.argv[5])
    num = int(sys.argv[6])
    cloudNum = int(sys.argv[7])
    colorMethod = int(sys.argv[8])
    vmax = int(sys.argv[9])
    fig = plt.figure(figsize=(4.5, 9))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    dgplot1(ax1, snrNum, level, lC, bC, hl, num)
    dgplot2(fig, ax2, snrNum, level, lC, bC, hl, num, cloudNum, colorMethod, vmax)
    plt.show()
#fig.savefig('../../data/map2018/snrOverlaps/snr{0}_overlap.pdf'.format(snrNum),
#            transparent=True, dpi=360)