from astropy.io import fits
from scipy.io import readsav
import numpy as np
res = readsav('../data/map2018/cemap_res15.sav')['cat']
l = fits.Column(name='L', format='D', array=res.l)
b = fits.Column(name='B', format='D', array=res.b)
EGK = fits.Column(name='EGK', format='19E', array=np.vstack(res.EGK))
DEGK = fits.Column(name='DEGK', format='19E', array=np.vstack(res.DEGK))
EBR = fits.Column(name='EBR', format='19E', array=np.vstack(res.EBR))
DEBR = fits.Column(name='DEBR', format='19E', array=np.vstack(res.DEBR))
EHK = fits.Column(name='EHK', format='19E', array=np.vstack(res.EHK))
DEHK = fits.Column(name='DEHK', format='19E', array=np.vstack(res.DEHK))
cols = fits.ColDefs([l, b, EGK, DEGK, EBR, DEBR, EHK, DEHK])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('../data/map2018/cemap_res15.fits', overwrite=True)
