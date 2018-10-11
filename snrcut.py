from astropy.io import fits
data, hdr = fits.getdata('../../data/map2018/catwithdis_galactic.fits', header=True)
# 1.选取ic443附近数据，保存
lc = 189
bc = 3
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr169new.fits', field, hdr, overwrite=True)
# 2.选取SNR 170 G190.9-2.2附近数据，保存
lc = 190.9
bc = -2.2
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr170new.fits', field, hdr, overwrite=True)
# 3.SNR172 Monecerous附近数据，保存
lc = 205.5
bc = 0.5
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr172new.fits', field, hdr, overwrite=True)
# 4.SNR 167 G182.4+4.3附近数据，保存
lc = 182.4
bc = 4.3
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr167new.fits', field, hdr, overwrite=True)
# 5.SNR 174 G213-0.6附近数据，保存
lc = 213
bc = -0.6
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr174new.fits', field, hdr, overwrite=True)
# 6.SNR160 G156.2+5.7附近数据，保存
lc = 156.2
bc = 5.7
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr160new.fits', field, hdr, overwrite=True)
# 7.SNR 159 G152.4-2.1附近数据，保存
lc = 152.4
bc = -2.1
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr159new.fits', field, hdr, overwrite=True)
# 8.SNR 162 G160.9+2.6附近数据，保存
lc = 160.9
bc = 2.6
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr162new.fits', field, hdr, overwrite=True)
# 9.SNR 163 G166+4.3附近数据，保存
lc = 166
bc = 4.3
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr163new.fits', field, hdr, overwrite=True)
# 10.SNR 164 G178.2-4.2附近数据，保存
lc = 178.2
bc = -4.2
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr164new.fits', field, hdr, overwrite=True)
# 11.SNR 165 G179+2.6附近数据，保存
lc = 179
bc = 2.6
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr165new.fits', field, hdr, overwrite=True)
# 12.SNR173 G206.9+2.3附近数据，保存
lc = 206.9
bc = 2.3
r = 3
index = (data.l >= lc - r) & (data.l <= lc + r) & (data.b >= bc - r) & (data.b <= bc + r)
field = data[index]
fits.writeto('../../data/map2018/snrs/snr173new.fits', field, hdr, overwrite=True)
