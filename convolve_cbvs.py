import numpy as np
from tqdm import trange
import os, sys
from astropy.io import fits
import eleanor
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut
from glob import glob

sectors = [int(sys.argv[1])]

for sector in sectors:
    """
    # Gets FFI time from MAST TESSCut
    coord = SkyCoord('04:35:50.330 -64:01:37.33', unit=(u.hourangle, u.deg))
    sector_table = Tesscut.get_sectors(coord)
    manifest = Tesscut.download_cutouts(coord, 31, sector=sector)
    cutout = fits.open(manifest['Local Path'][0])
    time = cutout[1].data['TIME'] - cutout[1].data['TIMECORR']
    """
    # get times from the postcards
    postdir = f'/data/tessraid/data/local_eleanor/postcards/s{sector:04d}'

    cbv_dir = './cbvs/'
    files = glob(os.path.join(cbv_dir, f'*s{sector:04d}*cbv.fits'))

    for c in trange(len(files)):
        cbv = fits.open(files[c])
        camera = cbv[1].header['CAMERA']
        ccd = cbv[1].header['CCD']
        cbv_time = cbv[1].data['Time']

        postfile = os.path.join(postdir, f'{camera}-{ccd}')
        if not os.path.exists(postfile):
            print(f'dir {postfile} does not exist; skipping')
            continue

        postfile = glob(os.path.join(postfile, '*[0-9].fits'))[0]
        hdu = fits.open(postfile)
        tstart = hdu[1].data['tstart']
        tend = hdu[1].data['tstop']

        new_fn = './metadata/s{0:04d}/cbv_components_s{0:04d}_{1:04d}_{2:04d}.txt'.format(
            sector, camera, ccd)

        convolved = np.zeros((len(tstart), 16))
        inds = np.array([], dtype=int)
        for i in range(len(tstart)):
            g = np.where((cbv_time >= tstart[i]) & (cbv_time <= tend[i]))[0]
            assert len(g) == 15
            for j in range(16):
                index = 'VECTOR_{0}'.format(j + 1)
                convolved[i, j] = np.mean(cbv[1].data[index][g])
        np.savetxt(new_fn, convolved)
        cbv.close()
        hdu.close()
