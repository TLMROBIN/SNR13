from astropy.io import fits
from scipy.io import readsav
from SNR13.radioar import moddis
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


snrNum = sys.argv[1]
level = int(sys.argv[2])
lC = float(sys.argv[3])
bC = float(sys.argv[4])
hl = float(sys.argv[5])
num = int(sys.argv[6])
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
image_data[image_data < 0] = 0
image_data = np.sqrt(image_data)
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
# 初步整理消光的距离间隔
AG = degk  + 1.987 * dehk  # G波段消光
disBins = [[moddis(i - 0.25), moddis(i+0.25)] for i in midMum]
# 讲前4个bin，第五到第八个bin合并
AG_new = np.zeros([AG.shape[0]-6, AG.shape[1], AG.shape[2]])
AG_new[0, :, :] = np.sum(AG[0:4, :, :], axis=0)
AG_new[1, :, :] = np.sum(AG[4:8, :, :], axis=0)
AG_new[2:, :, :] = AG[8:, :, :]
AG = AG_new
disBin1 = [disBins[0][0], disBins[3][1]]
disBin2 = [disBins[4][0], disBins[7][1]]
del disBins[0:8]
disBins.insert(0, disBin1)
disBins.insert(1, disBin2)
# 绘制消光/距离等值线图
binValues = [j - i for i, j in disBins] # 计算间隔大小
AG_kpc = np.zeros_like(AG)
for i in range(AG.shape[0]):
    AG_kpc[i, :, :] = AG[i, :, :] / binValues[i]
zmax = np.percentile(AG_kpc, 98)
levs = zmax / level * (np.arange(level) + 1)[1:]

# 画叠加图
fig = plt.figure()
ax = fig.add_subplot(111)
im = ax.imshow(image_data, origin='lower', interpolation='nearest', vmin=0, vmax=np.sqrt(1187), cmap='jet', extent=extent)
ax.set_xlim(lC + hl, lC - hl)
ax.set_ylim(bC - hl, bC + hl)
fig.colorbar(im, ax=ax, extend='max')
extin = gaussian_filter(AG_kpc[num, :, :], 0.68)
cont = ax.contour(gl, gb, extin, levs, cmap='gray')
tit = r'{0:4.2f}~{1:4.2f}'.format(disBins[num][0], disBins[num][1])
ax.text(0.02, 1.05, tit, va='top', ha='left', transform=ax.transAxes, color='red', fontsize=13, zorder=10)
# 调整刻度
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax.tick_params(labelsize=10)
# 调整xylabel
fig.text(0.5, 0.02, 'Galactic Latitude', ha='center', va='bottom', fontsize=14)
fig.text(0.1, 0.5, 'Galactic Longitude', ha='left', va='center',
         rotation='vertical', fontsize=14)
# 读取边界点
#points = readpoints(snrNum, fig, 1)
#datainborder(snrNum, 1, points)
plt.show()
#fig.savefig('../../data/map2018/snrOverlaps/snr{0}_overlap.pdf'.format(snrNum),
#            transparent=True, dpi=360)
