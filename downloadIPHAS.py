from astroquery.vo_conesearch.conesearch import conesearch
from astropy import units as u
from astropy.coordinates import SkyCoord
center = SkyCoord(189.1 * u.deg, 3 * u.deg, frame='galactic')
search = conesearch(center=center,
                    radius=1,
                    verb=1,
                    catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
search.to_table().write('../data/snrpic/iphas-data.fits', format='fits', overwrite=True)