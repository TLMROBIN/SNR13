from scipy.io import readsav
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import matplotlib.ticker as ticker


snrNum = sys.argv[1]
level = int(sys.argv[2])
lC = float(sys.argv[3])
bC = float(sys.argv[4])
hl = float(sys.argv[5])
switch = int(sys.argv[6])
res = readsav('../../data/map2018/IDLs/extin3d015{0}.sav'.format(snrNum))
# 初步读取数据
midMum = res.exta.dis[0]  # 距离间隔的中间值，距离
gl = res.exta.gl[0]  # 网点的银经
gb = res.exta.gb[0]  # 网点的银维
# debr = np.stack(res.exta.debr[0].dar)  # 区间内的E(B-V)h值，三维数组，23层，每层是一个二维矩阵 
degk = np.stack(res.exta.degk[0].dar)
dehk = np.stack(res.exta.dehk[0].dar)
# 初步整理消光的距离间隔
AG = degk  # + 1.987 * dehk  # G波段消光
disBins = [[i - 0.125, i + 0.125] for i in midMum]
# 讲前4个bin，第五到第八个bin合并
AG_new = np.zeros([int(AG.shape[0] / 2), AG.shape[1], AG.shape[2]])
for i in range(int(AG.shape[0] / 2)):
    AG_new[i, :, :] = np.sum(AG[2 * i:2 * i +2, :, :], axis=0)
AG = AG_new

# 绘制消光/距离等值线图
zmax = np.percentile(AG, 98)
levs = zmax / level * (np.arange(level) + 1)[1:]
print(levs)
fig, axes = plt.subplots(
    nrows=4, ncols=2, sharex='col', sharey='row', figsize=(4.2, 6.2))
for i, ax in zip(range(8), axes.flat[:8]):
    extin = gaussian_filter(AG[i, :, :], 0.68)
    cont = ax.contour(gl, gb, extin, levs)
    tit = r'{0:4.2f}~{1:4.2f}'.format(disBins[i][0] * 2, disBins[i][1] * 2)
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
fig.text(0.5, 0.01, 'Galactic Latitude', ha='center', fontsize=10)
fig.text(0.01, 0.5, 'Galactic Longitude', va='center', rotation='vertical', fontsize=10)

plt.tight_layout()
if switch == 1:
    plt.show()
else:
    fig.savefig('../../data/map2018/snrSlices/snr{0}_slices.pdf'.format(snrNum),
           transparent=True, dpi=360, overwrite=True)
