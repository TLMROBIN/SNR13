import sys 
from astropy.io import fits
from astropy.io.fits import getdata
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from SNR13.radioar import moddis
import matplotlib.ticker as ticker
snrNum = sys.argv[1]
snrShape = [int(sys.argv[2]), int(sys.argv[2])]
level = int(sys.argv[3])
lC = float(sys.argv[4])
bC = float(sys.argv[5])
hl = float(sys.argv[6])
data, hdr = getdata('../data/map2018/cemap{0}.fits'.format(snrNum), header=True)
l, b = data.L.reshape(snrShape), data.B.reshape(snrShape)
mu = np.arange(4.25, 13.75, 0.5)
dis = moddis(mu)
disbin = np.zeros(19)
extinbin = np.zeros_like(data.AG)
for i in range(19):
    if i == 0:
        disbin[i] = dis[i]
        extinbin[:, i] = data.AG[:, i] / disbin[i]
    else:
        disbin[i] = dis[i] - dis[i - 1]
        extinbin[:, i] = (data.AG[:, i] - data.AG[:, i - 1]) / disbin[i]
zmax = np.percentile(extinbin, 98)
levs = zmax / level * np.arange(level)
fig, axes = plt.subplots(
    nrows=6, ncols=3, sharex='col', sharey='row', figsize=(5,9))
for i, ax in zip(range(19), axes.flat):
    if i == 18:
        continue
    extin = extinbin[:, i+1]
    extin = gaussian_filter(extin.reshape(snrShape), 0.68)
    cont = ax.contour(l, b, extin, levs[1:])
    tit = r'{0:4.2f}~{1:4.2f}'.format(dis[i], dis[i+1])
    ax.text(0.02, 1.125, tit, va='top', ha='left',transform=ax.transAxes, 
    color='red', fontsize=8, zorder=10)
    ax.set
    # 调整刻度
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax.tick_params(labelsize=7)
    ax.set_xlim(lC - hl, lC + hl)
    ax.set_ylim(bC - hl, bC + hl)
# 调整xylabel
fig.text(0.5, 0.05, 'Galactic Latitude', ha='center', va='bottom', fontsize=14)
fig.text(0.03, 0.5, 'Galactic Longitude', ha='left', va='center', rotation='vertical', fontsize=14)
plt.show()
# 画上距离