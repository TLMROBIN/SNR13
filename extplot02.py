from scipy.io import readsav
import sys
from SNR13.radioar import moddis
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import matplotlib.ticker as ticker
snrNum = sys.argv[1]
level = int(sys.argv[2])
lC = float(sys.argv[3])
bC = float(sys.argv[4])
hl = float(sys.argv[5])
res = readsav('../../data/map2018/IDLs/extin3d015{0}.sav'.format(snrNum))
# 初步读取数据
midMum = res.exta.dis[0] # 距离间隔的中间值，距离模数
gl = res.exta.gl[0]  # 网点的银经
gb = res.exta.gb[0]  # 网点的银维
# debr = np.stack(res.exta.debr[0].dar)  # 区间内的E(B-V)h值，三维数组，23层，每层是一个二维矩阵 
degk = np.stack(res.exta.degk[0].dar)
dehk = np.stack(res.exta.dehk[0].dar)
# 初步整理消光的距离间隔
AG = degk + 1.987 * dehk  # G波段消光
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
print(levs)
fig, axes = plt.subplots(
    nrows=4, ncols=3, sharex='col', sharey='row', figsize=(5, 9))
for i, ax in zip(range(12), axes.flat[:12]):
    extin = gaussian_filter(AG_kpc[i, :, :], 0.68)
    cont = ax.contour(gl, gb, extin, levs)
    tit = r'{0:4.2f}~{1:4.2f}'.format(disBins[i][0], disBins[i][1])
    ax.text(0.02, 1.125, tit, va='top', ha='left', transform=ax.transAxes, color='red', fontsize=8, zorder=10)
    # 调整刻度
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax.tick_params(labelsize=7)
    ax.set_xlim(lC - hl, lC + hl)
    ax.set_ylim(bC - hl, bC + hl)
    ax.invert_xaxis()
# 调整xylabel
fig.text(0.5, 0.05, 'Galactic Latitude', ha='center', va='bottom', fontsize=14)
fig.text(0.03, 0.5, 'Galactic Longitude', ha='left', va='center', rotation='vertical', fontsize=14)
plt.show()