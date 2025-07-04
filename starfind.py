import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.datasets import load_star_image
from photutils.detection import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from photutils.datasets import load_star_image
from XRBID.WriteScript import WriteReg
from astropy.wcs import WCS

import warnings
warnings.filterwarnings('ignore')

M66_image_path = '/Users/undergradstudent/Research/XRB-Analysis/Galaxies/M66/HST/M66_mosaic_uvis_f555w_drc_sci.fits'
hdu = fits.open(M66_image_path)
data = fits.getdata(M66_image_path)
subset = data[3500:4000, 3500:4000]

def show_galaxy():
    '''
    Plot an image of the M66 galaxy.
    '''
    plt.figure(figsize=(10, 10))
    plt.imshow(
        data, 
        cmap='viridis',
        vmin=0, 
        vmax=0.3,
        origin='lower',
        interpolation='nearest'
    )

def show_zoomed_image(data=subset, vmin=0, vmax=0.5):
    '''
    Show a zoomed in image of the M66 galaxy.
    '''
    plt.figure(figsize=(10, 10))
    plt.imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    plt.show()

def find_stars(
    size, 
    sensitivity,
    data=data,
    std_multiple=5,
    sigma=3,
    plot=True,
    radius=4,
    cmap='Greys',
    vmin=0,
    vmax=0.3,
    aperture_color='blue',
    create_regions=True
):
    '''Find all the stars in the M66 galaxy.'''
    fwhm = {
        'very low': 0.07,
        'low': 0.09,
        'moderately low': 0.02, 
        'moderate': 0.2,
        'moderately high': 0.9,
        'high': 1.5,
        'very high': 3
    }

    threshold = {
        'very low': 0.015,
        'low': 0.02,
        'moderately low': 0.08, 
        'moderate': 0.2,
        'moderately high': 0.2,
        'high': 0.5,
        'very high': 3
    }
    mean, median, std = sigma_clipped_stats(data, sigma=sigma)
    
    if sensitivity:
        daofind = DAOStarFinder(fwhm=fwhm[size], threshold=threshold[sensitivity])
    else:
        daofind = DAOStarFinder(fwhm=fwhm[size], threshold=5*std)
    objects = daofind(data)

    print(f"Found {len(objects)} stars")

    positions = np.transpose((objects['xcentroid'], objects['ycentroid']))

    # Create apertures around sources
    if plot:
        plt.figure(figsize=(10, 10))
        apertures = CircularAperture(positions, r=radius)
        plt.imshow(
            data,
            cmap=cmap,
            vmin=vmin, 
            vmax=vmax,
            origin='lower',
            interpolation='nearest'
        )
        apertures.plot(color=aperture_color)
        plt.show()
    
    if create_regions == True:
        try: wcs = WCS(hdu['SCI'].header)
        except: wcs = WCS(hdu['PRIMARY'].header)
        xcoord_img = objects['xcentroid'].tolist()
        ycoord_img = objects['ycentroid'].tolist()
        xcoords_fk5, ycoords_fk5 = wcs.wcs_pix2world(xcoord_img, ycoord_img, 1)
        WriteReg(sources=[xcoords_fk5, ycoords_fk5], coordsys="fk5", \
                outfile="M66_stars.reg", \
                radius=0.15, radunit="arcsec", label=objects["id"].tolist(), color='blue')

def find_stars_in_zoomed_image(
    size, 
    sensitivity,
    data=subset,
    std_multiple=5,
    sigma=3,
    plot=True,
    radius=4,
    cmap='Greys',
    vmin=0,
    vmax=0.3,
    aperture_color='blue'
):
    '''Find all the stars in the zoomed in image of the M66 galaxy.'''
    fwhm = {
        'very low': 0.07,
        'low': 0.09,
        'moderately low': 0.02, 
        'moderate': 0.2,
        'moderately high': 0.9,
        'high': 1.5,
        'very high': 3
    }

    threshold = {
        'very low': 0.015,
        'low': 0.02,
        'moderately low': 0.08, 
        'moderate': 0.2,
        'moderately high': 0.2,
        'high': 0.5,
        'very high': 3
    }
    mean, median, std = sigma_clipped_stats(data, sigma=sigma)
    
    if sensitivity:
        daofind = DAOStarFinder(fwhm=fwhm[size], threshold=threshold[sensitivity])
    else:
        daofind = DAOStarFinder(fwhm=fwhm[size], threshold=5*std)
    objects = daofind(data)

    print(f"Found {len(objects)} stars")

    positions = np.transpose((objects['xcentroid'], objects['ycentroid']))

    # Create apertures around sources
    if plot:
        plt.figure(figsize=(10, 10))
        apertures = CircularAperture(positions, r=radius)
        plt.imshow(
            data,
            cmap=cmap,
            vmin=vmin, 
            vmax=vmax,
            origin='lower',
            interpolation='nearest'
        )
        apertures.plot(color=aperture_color)
        plt.show()
