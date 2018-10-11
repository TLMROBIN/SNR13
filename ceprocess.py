from astropy.io import fits
from astropy.io.fits import getdata
# 读取cemap计算AG并添加后保存
cedata = getdata('../data/map2018/cemap.fits')
AG = cedata.EGK + cedata.EHK * 1.987
col_AG = fits.Column(name='AG', format='19E', array=AG)
orig_cols = cedata.columns
hdu = fits.BinTableHDU.from_columns(orig_cols + col_AG)
hdu.writeto('../data/map2018/cemap_AG.fits', overwrite=True)
# 读取添加了消光的星表，裁剪ic443附近的数据并保存
data, hdr = getdata('../data/map2018/cemap_AG.fits', 1, header=True)
index = (data.L >= 186) & (data.L <= 192) & (data.B >= 0) & (data.B <= 6)
fits.writeto('../data/map2018/cemap169.fits', data[index], hdr, overwrite=True)
