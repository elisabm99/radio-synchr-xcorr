# This is a sample Python script.
import sys
import logging
import numpy as np
import scipy as sp
import pandas as pd

from pytictoc import TicToc
from astropy.io import fits
from astropy.wcs import WCS

from astropy.coordinates import SkyCoord
import treecorr
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, least_squares
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import stats

from parameters import *


plt.rc('text', usetex = True)
plt.rc('font', **{'family' : "sans-serif"})
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

pi = np.pi
################################################################################################
################################################################################################
################################################################################################

###########################################################
#                 Handle catalogs                         #
###########################################################
# Put SDSSDR12 into Treecorr format
def format_sdssdr12_catalog():
    
    print('Loading big catalog')
    data = np.loadtxt(filename_full_catalog, delimiter=',')
    ras = data[:, 1]
    decs = data[:, 2]
    zs = data[:, 8]
    print(ras.shape)
    
    # idxs = np.where(decs >= 70)
    # ras = ras[idxs]
    # decs = decs[idxs]
    # zs = zs[idxs]
    # print(ras.shape)

    print('\nSaving big catalog')
    np.savetxt(filename_full_catalog_formatted, np.transpose([ras, decs, zs]), fmt=['%.12f', '%.12f', '%.6f'])

# Put SDSSDR16Q into Treecorr format
def format_sdssdr16q_catalog():

    hdul = fits.open(filename_full_catalog1)
    data = hdul[1].data
    ras1  = data['RA']
    decs1  = data['DEC']
    zs1 = data['Z']
    
    hdul = fits.open(filename_full_catalog2)
    data = hdul[1].data
    ras2  = data['RA']
    decs2  = data['DEC']
    zs2 = data['Z']
    
    ras = np.concatenate((ras1, ras2))
    decs = np.concatenate((decs1, decs2))
    zs = np.concatenate((zs1, zs2))

    idx = np.where(zs >= 1)[0]
    ras = ras[idx]
    decs = decs[idx]
    zs = zs[idx]

    arr = np.column_stack((ras, decs, zs))
    # print(arr.shape)
    np.savetxt(filename_full_catalog_formatted, arr, fmt=['%.12f', '%.12f', '%.12f'])

# Put UNWISE into Treecorr format
def strip_line(line):
    line = line.strip()
    line = " ".join(line.split())
    return line.split(' ')

def format_unwise_catalog():
    print('Loading big catalog')

    ras = []
    decs = []

    i = 0
    with open(filename_full_catalog1, "r") as fin:
        for line in fin:
            if i > 21:
                line = strip_line(line)
                ras.append(float(line[0]))
                decs.append(float(line[1]))
            i += 1
    i = 0
    with open(filename_full_catalog2, "r") as fin:
        for line in fin:
            if i > 21:
                line = strip_line(line)
                ras.append(float(line[0]))
                decs.append(float(line[1]))
            i += 1
    i = 0
    with open(filename_full_catalog3, "r") as fin:
        for line in fin:
            if i > 21:
                line = strip_line(line)
                ras.append(float(line[0]))
                decs.append(float(line[1]))
            i += 1
    
    ras = np.array(ras)
    decs = np.array(decs)
    zs = np.zeros(len(ras))
    
    arr = np.column_stack((ras, decs, zs))
    print(arr.shape)
    np.savetxt(filename_full_catalog_formatted, arr, fmt=['%.12f', '%.12f', '%.12f'])

# Put 2MPZ into Treecorr format
def format_2mpz_catalog():

     data = np.genfromtxt(filename_full_catalog, skip_header=1, delimiter=',')
     ras = data[:,9]
     decs = data[:,10]
     zs = data[:, 34]

     # idx = np.where(zs >= 0.1)[0]
     # ras = ras[idx]
     # decs = decs[idx]
     # zs = zs[idx]
     
     arr = np.column_stack((ras, decs, zs))
     np.savetxt(filename_full_catalog_formatted, arr, fmt=['%.12f', '%.12f', '%.12f'])

# Put PSZ2 into Treecorr format
def format_psz2_catalog():
    print('Loading big catalog')

    # idxs = np.where(decs >= 70)
    # ras = ras[idxs]
    # decs = decs[idxs]
    # zs = zs[idxs]
    # print(ras.shape)
    
    hdul = fits.open(filename_full_catalog)
    data = hdul[1].data
    ras  = data['RA']
    decs  = data['DEC']
    zs = data['REDSHIFT']

            
    ras = np.array(ras)
    decs = np.array(decs)
    zs = np.array(zs)
    
    print(ras[:5])
    print(decs[:5])
    print(zs[:5])
    print(ras.shape)

    idxs = np.where(zs > -1)
    ras = ras[idxs]
    decs = decs[idxs]
    zs = zs[idxs]
    print(ras.shape)
    
    
    
    print('\nSaving big catalog')
    np.savetxt(filename_full_catalog_formatted, np.transpose([ras, decs, zs]), fmt=['%.12f', '%.12f', '%.6f'])
   
# Create a smaller catalog to xcorrelate with small image in Treecorr format.
# This method is only ok for small images. Neglects projection effects
def create_small_catalog(isColdest):
    
    print('isColdest =', isColdest)
    if isColdest:
        filename_fits_small = filename_fits_small_coldest
        filename_small_catalog_formatted = filename_small_catalog_formatted_coldest
    else:
        filename_fits_small = filename_fits_small_reference
        filename_small_catalog_formatted = filename_small_catalog_formatted_reference
    
    c_dec, delta_dec, c_ra, delta_ra, wcs, f = read_header(filename_fits_small)

    # print('c_dec, delta_dec, c_ra, delta_ra =', c_dec, delta_dec, c_ra, delta_ra)
    
    print('Loading big catalog')
    data = np.loadtxt(filename_full_catalog_formatted, delimiter=' ')
    ras = data[:, 0]
    decs = data[:, 1]
    zs = data[:, 2]

    
    print('\nFinding incides')
    print('c_ra + delta_ra,  c_ra - delta_ra =', c_ra + delta_ra, c_ra - delta_ra)
    print(' c_dec + delta_dec,   c_dec - delta_dec =', c_dec + delta_dec, c_dec - delta_dec)
    idx = np.where(
        (ras < c_ra + delta_ra) & (ras > c_ra - delta_ra) & (decs < c_dec + delta_dec) & (decs > c_dec - delta_dec))
    print('length indices =', len(idx[0]))
    
    print('\nSaving small catalog')

    np.savetxt(filename_small_catalog_formatted, np.transpose([ras[idx], decs[idx], zs[idx]]), fmt=['%.12f', '%.12f', '%.6f'])

# Make an histogram of redshift of the sources in the full catalog, or the subset corresponding to the
# codest spot or reference field
def make_redshift_histo(isColdest, isFull=False):
    if isFull:
        data = np.loadtxt(filename_full_catalog_formatted, delimiter=' ')
        filename_zhist = filename_zhist_full
    elif isColdest:
        data = np.loadtxt(filename_small_catalog_formatted_coldest, delimiter=' ')
        filename_zhist = filename_zhist_coldest
    else:
        data = np.loadtxt(filename_small_catalog_formatted_reference, delimiter=' ')
        filename_zhist = filename_zhist_reference
    
    zs = data[:, 2]
    
    idx = np.where(zs >= -0.1)[0]
    # print(len(idx), len(zs), len(idx)/len(zs))
    zs = zs[idx]
    print(np.amin(zs), np.amax(zs))
    
    if catalog == 0:
        zmin = np.amin(zs) - dz / 2
    elif catalog == 1:
        zmin = - dz / 2
    elif catalog == 3:
        zmin = - dz / 2
    elif catalog == 4:
        zmin = - dz / 2
    
    binedges = np.arange(zmin, np.amax(zs) + dz / 2, dz)
    bins = (binedges[1:] + binedges[:-1]) / 2
    
    counts, binedges, patches = plt.hist(zs, binedges)
    
    np.savetxt(filename_zhist, np.transpose([bins, counts]), fmt=['%.2f', '%d'], delimiter='\t', header='z\t\t counts')

###########################################################
#                   Handle image                          #
###########################################################
def get_Omega_psf(isMosaic, isColdest):
    if isMosaic:
        BMAJ = 0.0571624149291658  # degrees
        BMIN = 0.055289331960945  # degree
    elif isColdest:
        BMAJ = 0.0626009017077947  # degrees
        BMIN = 0.054289551715023  # degrees
    else:
        BMAJ = 0.055774335667134  # degrees
        BMIN = 0.050173110977948  # degrees
    
    # return 9 / 4 * pi * BMIN * BMAJ / (8 * np.log(2) ** 2) * (pi / 180) ** 2  # sr^2
    return pi * BMIN * BMAJ / (4 * np.log(2)) * (pi / 180) ** 2  # sr^2


# Convert units of correlator to Kelvin^2
def get_T_factor(isMosaic, isColdest):
    if isMosaic:
        freq = 120603942.871094  # Hz
    else:
        freq = 142283630.371094  # Hz
    
    return 1e-26 * 2.99792458e8 ** 2 / (2 * 1.380649e-23 * freq ** 2 * get_Omega_psf(isMosaic, isColdest))


# Get CRVAL and range of ra, dec of image. Only works for small image (flat approx)
def read_header(file):
    f = fits.open(file)
    header = f[0].header
    # CRPIX1 = header['CRPIX1']
    # CRPIX2 = header['CRPIX2']
    NAXIS1 = header['NAXIS1']
    NAXIS2 = header['NAXIS2']
    CRVAL1 = header['CRVAL1']
    CRVAL2 = header['CRVAL2']
    CDELT1 = header['CDELT1']
    CDELT2 = header['CDELT2']
    
    delta_dec = np.abs(CDELT2 * NAXIS2)
    delta_ra = np.abs(CDELT1 * NAXIS1)
    
    wcs = WCS(header, naxis=2)
    
    return CRVAL2, delta_dec, CRVAL1, delta_ra, wcs, f


def load_image(file, isMosaic, field):
    
    f = fits.open(file)
    header = f[0].header
    w = WCS(header)  # first pixel is [0,0]
    NAXIS1 = header['NAXIS1']  # ra
    NAXIS2 = header['NAXIS2']  # dec
    print('NAXIS1, NAXIS2 =', NAXIS1, NAXIS2)
    if isMosaic and (field == -1 or field == 0):
        data = f[0].data.transpose()  # transpose so first dimension is ra, second dimension is dec
    else:
        data = f[0].data[0, 0].transpose()  # transpose so first dimension is ra, second dimension is dec
    
    print('data.shape =', data.shape)
    
    ras = np.empty(NAXIS1 * NAXIS2)
    decs = np.empty_like(ras)
    pixels_ra = np.empty_like(ras)
    pixels_dec = np.empty_like(decs)
    
    dec_pixels = np.arange(NAXIS2)
    
    print('Reading data')
    for i in np.arange(NAXIS1):
        # print(i*NAXIS2 , (i+1) * NAXIS2)
        # try:
        sky = w.pixel_to_world(i, dec_pixels, 0, 0)
        # print(i, sky[0].ra.deg, sky[0].dec.deg)
        
        # except AstropyUserWarning:
        #     print('##### ', sky)
        
        ras[i * NAXIS2: (i + 1) * NAXIS2] = sky[0].ra.deg
        decs[i * NAXIS2: (i + 1) * NAXIS2] = sky[0].dec.deg
        
        pixels_ra[i * NAXIS2: (i + 1) * NAXIS2] = i * np.ones(NAXIS2)
        pixels_dec[i * NAXIS2: (i + 1) * NAXIS2] = dec_pixels
        
    return data, ras, decs, pixels_ra, pixels_dec
    
# Create a text file from image for Treecorr. For each pixel there is one row: ra, dec, flux
def create_kappa_catalog(isMosaic, isColdest, field=-1):
    
    print('isMosaic, field =', isMosaic, field)
    if isMosaic:
        if field == -1:
            file = filename_fits_big
        elif field == 0:
            file = filename_fits_big_f0
        elif field == 1:
            file = filename_fits_big_f1
    else:
        if isColdest:
            file = filename_fits_small_coldest
        else:
            file = filename_fits_small_reference

        
    data, ras, decs, pixels_ra, pixels_dec = load_image(file, isMosaic, field)
    fluxes = data.flatten()  # flatten in row major order. The first NAXIS2 entries have the same RA (up to distorsion due to the projection)

    # add weights
    if isMosaic:
        # print('\nGetting weigths')
        # weights = np.empty_like(ras)
        # foot = fits.open(file_footprint_big)
        # foot_data = foot[0].data.transpose()
        # # print(foot_data.shape)
        # for i in np.arange(NAXIS1):
        #     weights[i*NAXIS2 : (i+1) * NAXIS2] = foot_data[i]**2
        # weights = weights / np.nanmax(weights)
        # # idxs = np.where(decs >  85)
        # # weights[idxs] = 0
        weights = np.ones_like(ras)
        # idxs = np.where(decs >  85)
        # weights[idxs] = 0

    else:
        weights = np.ones_like(ras)
        
    idx = np.where(~np.isnan(fluxes))
    print(idx[0].shape)
    ras = ras[idx]
    decs = decs[idx]
    fluxes = fluxes[idx]
    weights = weights[idx]
    pixels_ra = pixels_ra[idx]
    pixels_dec = pixels_dec[idx]

    
    print('Saving data to file ',file[:-4] +'dat')
    np.savetxt(file[:-4] +'dat', np.transpose([ras, decs, fluxes, weights, pixels_ra, pixels_dec]), fmt=['%.12f', '%.12f', '%.12f', '%.12f', '%d', '%d'])

# Get a rotation invariant PSF profile from an image. Also normalize integral of PSF to 1
def get_psf_profile(isMosaic, isColdest):
    if isMosaic:
        filename_psf_fits = filename_psf_big
        filename_psf_profile = filename_psf_profile_big
    else:
        if isColdest:
            filename_psf_fits = filename_psf_small_coldest
            filename_psf_profile = filename_psf_profile_small_coldest
        else:
            filename_psf_fits = filename_psf_small_reference
            filename_psf_profile = filename_psf_profile_small_reference
        
    print('PSF file:', filename_psf_fits)
    f = fits.open(filename_psf_fits)
    header = f[0].header
    w = WCS(header) # first pixel is [0,0]
    NAXIS1 = header['NAXIS1'] # ra in deg
    NAXIS2 = header['NAXIS2'] # dec in deg
    CRPIX1 = int(header['CRPIX1']) - 1 # astropy says pixel coordinates should be zero-based
    CRPIX2 = int(header['CRPIX2']) - 1
    
    psf = f[0].data.transpose()[:, :, 0, 0]
    # psf = f[0].data[0, 0]
    print('psf.shape =', psf.shape)
    center = w.pixel_to_world(CRPIX1, CRPIX2, 0, 0)[0]
    print('Center of the image\n', center)
    print('Value of PSF at central pixel =', psf[CRPIX1, CRPIX2])

    print('\nCalculating distances from the center of the image')
    dists = np.empty_like(psf)
    
    dec_pixels = np.arange(NAXIS2)
    for i in np.arange(NAXIS1): # loop over ras
        sky = w.pixel_to_world(i, dec_pixels, 0, 0)[0]
        dists[i] = sky.separation(center) # deg
    center_value =  psf[CRPIX1, CRPIX2]
    print('center value =', center_value)

    # Convert distances to radians
    dists = dists * pi / 180 # rad

    # Calculate average flux in concentric rings around center
    # Choose  bins
    max_dist = max_sep_psf / 60 * pi / 180 # rad
    min_dist = min_sep_psf / 60 * pi / 180 # rad
    print(min_dist, max_dist)
    binedges = np.array(np.concatenate(([0], np.logspace(np.log10(min_dist), np.log10(max_dist), num=nbins_psf+1, endpoint=True) )) ) # rad
    # binedges = np.logspace(np.log10(min_dist), np.log10(max_dist), num=nbins_psf+1, endpoint=True) # rad
    bins = 10**((np.log10(binedges)[1:] + np.log10(binedges)[:-1]) / 2) # rad


    print('binedges.shape =', binedges.shape )
    print('bins.shape =', bins.shape )

    print('\nGetting profile')
    print('Digitizing')
    dists_flat = dists.flatten()
    psf_flat = psf.flatten()

    mean_flux, binedges, binnumber = stats.binned_statistic(dists_flat, values=psf_flat, statistic ='mean', bins=binedges)
    
    # Normalize integral of PSF to one
    s = pd.Series(mean_flux)
    s = s.interpolate(limit_area=None,  method='quadratic') # using this because it interpolates over NaNs automatically
    mean_flux = s.to_numpy()
    
    # Throw away zeroth bin after using it for interpolation over nans
    binedges = binedges[1:]
    bins = bins[1:]
    mean_flux = mean_flux[1:]
    
    norm = 2*pi * np.trapz(np.sin(bins) * mean_flux, x=bins)
    mean_flux_norm = mean_flux / norm
    center_value = center_value / norm
    print('integral of PSF =', 2*pi * np.trapz(np.sin(bins) * mean_flux_norm, x=bins))


    spline = InterpolatedUnivariateSpline(bins,
                                          2 * pi * np.sin(bins) * mean_flux_norm,
                                          bbox=[0, pi], k=1)
    print('integral of PSF with spline=', spline.integral(bins[0], bins[-1]) )


    plt.plot(bins * 180/pi * 60, mean_flux_norm, marker='.')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel(r'$\theta~[\mathrm{arcmin}]$')
    plt.ylabel(r'$\mathrm{PSF~normalized}$')
    plt.show()
    print('Saving profile')
    np.savez(filename_psf_profile % (min_sep_psf, max_sep_psf, nbins_psf), profile=mean_flux_norm, binedges=binedges, center_flux=center_value, profile_nonorm=mean_flux)
  
# Xcorrelate image (K) and source catalog (N)
def xcorrelate(isMosaic, isColdest, field=-1):
    
    if isMosaic:
        metric = 'Arc' #'Euclidean' #
        my_min_sep = min_sep_big

        filename_N_catalog = filename_full_catalog_formatted
        filename_cov_cross = filename_cov_cross_big

        if field == -1:
            filename_fits_catalog = filename_fits_catalog_big
            filename_cross = filename_cross_big
        elif field == 0:
            filename_fits_catalog = filename_fits_catalog_big_f0
            filename_cross = filename_cross_big_f0
        elif field == 1:
            filename_fits_catalog = filename_fits_catalog_big_f1
            filename_cross = filename_cross_big_f1
    else:
        metric = 'Euclidean'
        my_min_sep = min_sep
        if isColdest:
            filename_fits_catalog = filename_fits_catalog_small_coldest
            filename_N_catalog = filename_small_catalog_formatted_coldest_no_extra
            filename_cross = filename_cross_small_coldest
            filename_cov_cross = filename_cov_cross_small_coldest
        else:
            filename_fits_catalog = filename_fits_catalog_small_reference
            filename_N_catalog = filename_small_catalog_formatted_reference
            filename_cross = filename_cross_small_reference
            filename_cov_cross = filename_cov_cross_small_reference
    print('Reading catalog from FITS')
    cat_K = treecorr.Catalog(filename_fits_catalog, config_K, npatch=npatch, logger=logger)
    
    print('Reading SDSS catalog')
    cat_N = treecorr.Catalog(filename_N_catalog, config_N, patch_centers=cat_K.patch_centers, logger=logger)

    nk = treecorr.NKCorrelation(min_sep=my_min_sep, max_sep=max_sep, nbins=nbins, sep_units='arcmin',  metric=metric, var_method=var_method, verbose=2, logger=logger, brute=False)
    nk.process(cat_N, cat_K)
    np.savez(filename_cov_cross % (my_min_sep, max_sep, nbins, var_method, npatch), cov=nk.cov)
    nk.write(filename_cross % (my_min_sep, max_sep, nbins, var_method, npatch))
    np.savetxt(filename_cov_cross[:-4] % (my_min_sep, max_sep, nbins, var_method, npatch) + '.dat', nk.cov, delimiter='\t')

def auto_correlate(isMosaic, isColdest, field=-1):
    if isMosaic:
        filename_cov_auto = filename_cov_auto_big
        my_min_sep = min_sep_big
        if field == -1:
            filename_fits_catalog = filename_fits_catalog_big
            filename_auto = filename_auto_big
        elif field == 0:
            filename_fits_catalog = filename_fits_catalog_big_f0
            filename_auto = filename_auto_big_f0
        elif field == 1:
            filename_fits_catalog = filename_fits_catalog_big_f1
            filename_auto = filename_auto_big_f1
    else:
        my_min_sep = min_sep

        if isColdest:
            filename_fits_catalog = filename_fits_catalog_small_coldest
            filename_cov_auto = filename_cov_auto_coldest
            filename_auto = filename_auto_coldest
        else:
            filename_fits_catalog = filename_fits_catalog_small_reference
            filename_cov_auto = filename_cov_auto_reference
            filename_auto = filename_auto_reference
    

    cat_K = treecorr.Catalog(filename_fits_catalog, config_K, npatch=npatch, logger=logger)

    kk = treecorr.KKCorrelation(min_sep=my_min_sep, max_sep=max_sep, nbins=nbins, sep_units='arcmin', var_method=var_method)#, bin_slop=0.5)
    kk.process(cat_K)
    print(kk.xi)

    np.savez(filename_cov_auto % (my_min_sep, max_sep, nbins, var_method, npatch), cov=kk.cov)
    kk.write(filename_auto % (my_min_sep, max_sep, nbins, var_method, npatch))

def gauss(x, *p):
        A, sigma = p
        return A * np.exp(-x ** 2 / (2. * sigma ** 2))

def normalizeToZero(xi, rs):
    print('Normalizing')
    norm = np.trapz(xi, x=rs)
    # print('Integral before subtraction =', np.trapz(xi, x=rs))
    # print(xi)
    xi = xi - norm / (rs[-1] - rs[0])
    # print(xi)
    print('Integral after subtraction =', np.trapz(xi, x=rs))
    return xi

# Method to get Table 1
def fit_psf(isColdest):
    print('\n\n**************************')
    print('Catalog =', catname)
    print('isColdest =', isColdest)
    if isColdest:
        filename_psf_profile = filename_psf_profile_small_coldest_marco
        filename_cross_small = filename_cross_small_coldest
        filename_cov_cross = filename_cov_cross_small_coldest
    else:
        filename_psf_profile = filename_psf_profile_small_reference_marco
        filename_cross_small = filename_cross_small_reference
        filename_cov_cross = filename_cov_cross_small_reference

    print(filename_psf_profile)
    print(filename_cross_small% (min_sep, max_sep, nbins, var_method, npatch))
    print(filename_cov_cross% (min_sep, max_sep, nbins, var_method, npatch))
    # print('Loading measurement')
    data = np.loadtxt(filename_cross_small % (min_sep, max_sep, nbins, var_method, npatch), skiprows=2)
    rs = data[:, 0] / 60 * pi / 180  # rad
    xi = data[:, 3]
    cov = np.load(filename_cov_cross % (min_sep, max_sep, nbins, var_method, npatch))['cov']#[:2, :2]
    print(cov.shape)
    covavg = np.mean(cov)
    # print(cov/np.mean(cov))
    # cov = np.diag(np.diag(cov))
    print('Inverting covariance matrix')
    covinv = np.linalg.inv(cov)
    # covinv = sp.linalg.inv(cov)

    data = np.loadtxt(filename_psf_profile)
    bins = data[:,0] / 60 * pi / 180  # rad
    psf = data[:,1]
    # psf = np.abs(psf)
    psf = psf * xi[0] / psf[0]
    psf_func = interp1d(bins, psf, bounds_error=False, fill_value=0)
    npar_psf = 1


    print('Loading Template')
    data = np.loadtxt(filename_template)
    angles = data[:, 0] / 60 * pi / 180  # rad
    curve = data[:, 1]

    curve = curve * xi[0] / curve[0]
    template_func = interp1d(angles, curve,  bounds_error=False, fill_value=0)
    A = np.array([psf_func(rs).T,  template_func(rs).T]).T
    npar_templ = 2

    
    # Fitting PSF only
    print('\n* Fitting PSF only')
    def func_chi(theta):
        x = xi - psf_func(rs) * theta
        chi2 = np.dot(x.T, np.dot(covinv, x))
        return np.sqrt(chi2)
    
    result = least_squares(func_chi, [1], ftol=1e-8)
    print(result)
    
    norm = np.dot(psf_func(rs).T, np.dot(covinv, psf_func(rs)))
    thetabest_psf = np.dot(psf_func(rs).T, np.dot(covinv, xi)) / norm
    print('Best fit PSF', thetabest_psf)
    x = xi - psf_func(rs) * thetabest_psf
    chi2_best_psf = np.dot(x.T, np.dot(covinv, x))
    print('chi2-best =', chi2_best_psf, 2 * result.cost)
    print()
    # print( (xi[0] - xi[1]) / (psf_func(rs[0]) - psf_func(rs[1])) )
    
    # print(np.dot(psf_func(rs).T, np.dot(covinv, xi)))
    # print(norm)
    # for theta in np.linspace(2e-4, 2.5e-4, num=3):
    #     print()
    #     print(theta, func_chi(theta)**2)
    #     print(np.sum(xi - theta * psf_func(rs)))
    #     print(xi - theta * psf_func(rs))
    #     print((np.dot(covinv, xi - theta * psf_func(rs))))

    # for i in range(10):
    #     plt.plot(covinv[i], label=str(i))
    # for i in range(10,20):
    #     plt.plot(covinv[i], linestyle=':', label=str(i))
    # plt.legend(loc=4)
    # plt.show()
    #
    # plt.plot(rs, xi, label='data')
    # plt.plot(rs, thetabest_psf * psf_func(rs), label='b.f. PFS')
    # # plt.plot(rs, psf_func(rs), label='PFS ini')
    # # plt.plot(rs, Abest[0]*psf_func(rs) + Abest[1] * template_func(rs), label='b.f. PFS + template')
    # plt.xscale('log')
    # plt.legend()
    # # plt.yscale('log')
    # # plt.show()
    # plt.close()

    # Fitting PSF + template
    
    
    print('\n* Fitting PSF + template 2-halo')
    def func_chi_templ(theta):
        x = xi - theta[0] * psf_func(rs) - theta[1] * template_func(rs)
        chi2 = np.dot(x.T, np.dot(covinv, x))
        # chi2 = np.dot(x.T, x  /  np.diag(cov))
        return np.sqrt(chi2)

    result = least_squares(func_chi_templ, [1, 1], ftol=1e-8)
    print(result)
    
    
    norm = np.dot(A.T, np.dot(covinv, A))
    Abest = np.dot(np.linalg.inv(norm), np.dot(A.T, np.dot(covinv, xi)))
    print('Best fit PSF + template 2-halo', Abest)
    x = xi - Abest[0] * psf_func(rs) - Abest[1] * template_func(rs)
    chi2_best_templ = np.dot(x.T, np.dot(covinv, x))
    print('chi2-best =', chi2_best_templ, 2 * result.cost)


    print('\nUsing Wilks\' theorem')
    deltachi2 = chi2_best_psf - chi2_best_templ
    c = stats.chi2.cdf(deltachi2, 1)
    print('deltachi2 =', deltachi2)
    print('p-value =', 1-c)
    
    
    if catalog == 3:
        print('\n* Fitting PSF + template 2h + template 1h')
        data = np.loadtxt(filename_template_1h)
        angles = data[:, 0] / 60 * pi / 180  # rad
        curve = data[:, 1]
        # curve = curve * xi[0] / curve[0]
        template_func_1h = interp1d(angles, curve, bounds_error=False, fill_value=0)
        A = np.array([psf_func(rs).T, template_func(rs).T, template_func_1h(rs).T]).T

        def func_chi_templ_1h(A):
            x = xi - A[0] * psf_func(rs) - A[1] * template_func(rs) -  A[2] * template_func_1h(rs)
            chi2 = np.dot(x.T, np.dot(covinv, x))
            return np.sqrt(chi2)
    
        result = least_squares(func_chi_templ_1h, [0.5, 0.1, 0.5], ftol=1e-8)
        print(result)
        
        
        norm = np.dot(A.T, np.dot(covinv, A))
        # print(norm)
        Abest = np.dot(np.linalg.inv(norm), np.dot(A.T, np.dot(covinv, xi)))
        print('Best fitPSF + template 2h + template 1h', Abest)
        lambd = xi - Abest[0] * psf_func(rs) - Abest[1] * template_func(rs) - Abest[2] * template_func_1h(rs)
        chi2_best_templ = np.dot(lambd.T, np.dot(covinv, lambd))
        print('chi2-best =', chi2_best_templ, 2 * result.cost)

    

        print('\nUsing Wilks\' theorem')
        deltachi2 = chi2_best_psf - chi2_best_templ
        c = stats.chi2.cdf(deltachi2, 2)
        print('deltachi2 =', deltachi2)
        print('p-value =', 1 - c)


    # plt.plot(rs, xi, label='data')
    # plt.plot(rs, thetabest_psf * psf_func(rs), label='b.f. PFS')
    # # plt.plot(rs, Abest[0]*psf_func(rs) + Abest[1] * template_func(rs), label='b.f. PFS + template')
    # plt.xscale('log')
    # plt.legend()
    # # plt.yscale('log')
    # # plt.show()
    # plt.close()
    

###########################################################
#                    Plotting                            #
###########################################################
# Plots xcorrelator, psf, and error band
def plot_simple(isAuto, isMosaic, isColdest, field=-1):
    print('isMosaic =', isMosaic)
    print('field =', field)
    if not isAuto:
        if isMosaic:
            my_min_sep = min_sep_big
            filename_psf_profile = filename_psf_profile_big
            if field == -1:
                filename_in = filename_cross_big
                filename_plot = filename_plot_cross_big
            elif field == 0:
                filename_in = filename_cross_big_f0
                filename_plot = filename_plot_cross_big_f0
            elif field == 1:
                filename_in = filename_cross_big_f1
                filename_plot = filename_plot_cross_big_f1
        else:
            my_min_sep = min_sep
            if isColdest:
                filename_in = filename_cross_small_coldest
                filename_plot = filename_plot_cross_small_coldest
                filename_psf_profile = filename_psf_profile_small_coldest
            else:
                filename_in = filename_cross_small_reference
                filename_plot = filename_plot_cross_small_reference
                filename_psf_profile = filename_psf_profile_small_reference
    else:
        if isMosaic:
            my_min_sep = min_sep_big
            filename_psf_profile = filename_psf_profile_big
            if field == -1:
                filename_in = filename_auto_big
                filename_psf_profile = filename_psf_profile_big
                filename_plot = filename_plot_auto_big
            elif field == 0:
                filename_in = filename_auto_big_f0
                filename_psf_profile = filename_psf_profile_big
                filename_plot = filename_plot_auto_big_f0
            elif field == 1:
                filename_in = filename_cross_big_f1
                filename_psf_profile = filename_psf_profile_big
                filename_plot = filename_plot_auto_big_f1
        else:
            my_min_sep = min_sep
            if isColdest:
                filename_in = filename_auto_coldest
                filename_plot = filename_plot_auto_coldest
                filename_psf_profile = filename_psf_profile_small_coldest
            else:
                filename_in = filename_auto_reference
                filename_plot = filename_plot_auto_reference
                filename_psf_profile = filename_psf_profile_small_reference
    
    print('filename_in =', filename_in% (my_min_sep, max_sep, nbins, var_method, npatch))
    print('filename_psf_profile =', filename_psf_profile % (min_sep_psf, max_sep_psf, nbins_psf))
    # Load correlator
    data = np.loadtxt(filename_in % (my_min_sep, max_sep, nbins, var_method, npatch), skiprows=2)
    rs = data[:, 0]  # arcmin
    xi = data[:, 3]
    varxi = data[:, 4]
    # print('rs =', rs)
    # print('corr = ', xi)
    
    # Load psf
    data_psf = np.load(filename_psf_profile % (min_sep_psf, max_sep_psf, nbins_psf))
    binedges = data_psf['binedges']  # rad
    bins = (np.log10(binedges[:-1]) + np.log10(binedges[1:])) / 2
    bins = 10 ** bins * 180 / pi * 60  # arcmin
    
    profile = data_psf['profile']
    factor = np.amax(xi) / profile[0]
    profile = profile * factor
    
    label = 'auto corr' if isAuto else 'xcorr'
    plt.plot(rs, xi, label=label)
    # plt.plot(rs, corr_fit, linestyle=':')
    plt.fill_between(rs, xi - varxi, xi + varxi, alpha=0.5, facecolor='cyan')
    
    plt.plot(bins, profile, label='psf')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel(r'$r~\mathrm{[arcmin]}$')
    ylabel = r'$\xi(r)~\mathrm{[(Jy/beam)^2]}$' if isAuto else r'$\xi(r)~\mathrm{[Jy/beam]}$'
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig(filename_plot % (my_min_sep, max_sep, nbins, var_method, npatch), dpi=300, bbox_inches="tight")
    plt.close()

# Compare PSF for the various images
def plot_compare_psf():
    
    data = np.load(filename_psf_profile_small_coldest % (min_sep_psf, max_sep_psf, nbins_psf))
    binedges_coldest = data['binedges'] * 180 / pi * 60  # arcmin
    bins_coldest = (np.log10(binedges_coldest[:-1]) +np.log10( binedges_coldest[1:])) / 2
    bins_coldest = 10**bins_coldest
    profile_coldest = data['profile_nonorm']
    
    data = np.load(filename_psf_profile_small_reference % (min_sep_psf, max_sep_psf, nbins_psf))
    binedges_reference = data['binedges'] * 180 / pi * 60  # arcmin
    bins_reference = (np.log10(binedges_reference[:-1]) +np.log10( binedges_reference[1:])) / 2
    bins_reference = 10**bins_reference
    profile_reference = data['profile_nonorm']
    
    data = np.load(filename_psf_profile_big % (min_sep_psf, max_sep_psf, nbins_psf))
    binedges_mosaic = data['binedges'] * 180 / pi * 60  # arcmin
    bins_mosaic = (np.log10(binedges_mosaic[:-1]) +np.log10( binedges_mosaic[1:])) / 2
    bins_mosaic = 10**bins_mosaic
    profile_mosaic = data['profile_nonorm']

    plt.plot(bins_coldest, profile_coldest, label='Coldest spot')
    plt.plot(bins_reference, profile_reference, label='Reference')
    plt.plot(bins_mosaic, profile_mosaic, label='Mosaic')

    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel(r'$\theta~\mathrm{[arcmin]}$')
    plt.ylabel(r'PSF$~\mathrm{[Jy/beam]}$')
    plt.legend()
    plt.show()
    

################################################################################################
################################################################################################
################################################################################################
if __name__ == '__main__':
    
    ## Catalog
    # format_sdssdr12_catalog() # Keeps full catalog, just format
    # format_sdssdr16q_catalog()
    # format_unwise_catalog()
    # format_2mpz_catalog()
    # format_psz2_catalog()

    # create_small_catalog(isColdest=True) # Reduces dimension of the catalog
    # create_small_catalog(isColdest=False) # Reduces dimension of the catalog

    # make_redshift_histo(isColdest=True, isFull=False)


    ## Image. (Mosaic was not used in the paper)
    # get_psf_profile(isMosaic=False, isColdest=True)
    # plot_compare_psf()

    myfield = -1 # choose a specific field for the mosaic. If -1, use whole mosaic
    # create_kappa_catalog(isMosaic=False, isColdest=True, field=myfield)


    ## Auto correlation image
    # auto_correlate(isMosaic=False, isColdest=True, field=myfield)
    # plot_simple(isAuto=True, isMosaic=False, isColdest=True, field=myfield)

    ## Cross correlation
    # xcorrelate(isMosaic=False, isColdest=True, field=myfield)
    # plot_simple(isAuto=False, isMosaic=False, isColdest=True, field=myfield)

    ## Method to get Table 1
    fit_psf(isColdest=True)

