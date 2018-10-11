from SNR13.extfit import exfit
from SNR13.dgplot03 import dgplot1
from SNR13.dgplot03 import dgplot2
from SNR13.dgplot03 import dgplot3
import sys
import matplotlib.pyplot as plt
snrNum = sys.argv[1]
level = int(sys.argv[2])
lC = float(sys.argv[3])
bC = float(sys.argv[4])
hl = float(sys.argv[5])
num = int(sys.argv[6])
cloudNum = int(sys.argv[7])
colorMethod = int(sys.argv[8])
vmax = int(sys.argv[9])

fig = plt.figure(figsize=(22, 4.5))
ax1 = fig.add_subplot(151)
ax2 = fig.add_subplot(152)
ax3 = fig.add_subplot(153)

ax4 = fig.add_subplot(154)
ax5 = fig.add_subplot(155)
#ax6 = fig.add_subplot(336)

#ax8 = fig.add_subplot(338)
#ax9 = fig.add_subplot(339)

dgplot3(ax1, snrNum, lC, bC, hl, colorMethod, vmax)
dgplot1(ax2, snrNum, level, lC, bC, hl, num)
dgplot1(ax3, snrNum, level, lC, bC, hl, num + 1)
dgplot1(ax4, snrNum, level, lC, bC, hl, num + 2)
dgplot1(ax5, snrNum, level, lC, bC, hl, num + 3)

#exfit(ax3, snrNum, cloudNum)

#dgplot1(ax1, snrNum, level, lC, bC, hl, num)
#dgplot2(fig, ax2, snrNum, level, lC, bC, hl, num, cloudNum, colorMethod, vmax)
#exfit(ax3, snrNum, cloudNum)


#dgplot2(fig, ax8, snrNum, level, lC, bC, hl, num, cloudNum + 2, colorMethod, vmax)
#exfit(ax9, snrNum, cloudNum + 2)

plt.tight_layout()
fig.savefig('../../data/map2018/snrOverlaps/snr{0}_overlap.pdf'.format(snrNum),
            transparent=True, dpi=360, overwrite=True)
