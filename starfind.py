import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.datasets import load_star_image
from photutils.detection import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture

def show(image, vmin=0, vmax=10, colourmap='viridis', figsize=(10, 8)):
    '''A function to plot an image in the fits file'''
    data = fits.getdata(image)
    plt.figure(figsize=figsize)
    plt.imshow(data, vmin=vmin, vmax=vmax, origin='lower', cmap=colourmap)
    plt.show()

def show_example_image():
    '''A function to plot an example image of stars.'''
    hdu = fits.getdata('M66_hubble.fits')[6000:6500, 2000:2500]
    plt.figure(figsize=(8, 8))
    plt.imshow(hdu, origin='lower', vmin=0, vmax=0.3)
    plt.show()



def find_stars(data, size=3, sensitivity=None, cmap='Greys'):
    '''A utility function to find stars in an image.'''
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    if not sensitivity:
        daofind = DAOStarFinder(fwhm=size, threshold=5.*std)
    else: daofind = DAOStarFinder(fwhm=size, threshold=sensitivity)
    sources = daofind(data - median)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    print(f"Your code found {len(positions)} stars...\nLet's plot them!")
    apertures = CircularAperture(positions, r=4.0)
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.figure(figsize=(8, 8))
    plt.imshow(data, cmap=cmap, origin='lower', norm=norm,
               interpolation='nearest')
    apertures.plot(color='blue', lw=1.5, alpha=0.5)


if __name__ == "__main__":
    show_example_image()
    plt.show()
    data = load_star_image().data[0:401, 0:401]
    find_stars(data, size=3)
    plt.show()
