config_N = {'ra_col' : 1,
           'dec_col' : 2,
           'ra_units' : 'deg',
           'dec_units' : 'deg'
            ,'verbose' : 1
            }

config_K = {'ra_col' : 1,
           'dec_col' : 2,
           'ra_units' : 'deg',
           'dec_units' : 'deg',
            'k_col' : 3
    # ,        'w_col ' : 4
            ,'wpos_col' :4
            ,'verbose' : 1
            }

# Binning for correlator
min_sep = 0.12 # arcmin
max_sep = 7 * 60.  # arcmin
nbins = 20
min_sep_big = 0.3   # arcmin

# Binning for PSF
min_sep_psf = 0.10   # arcmin
max_sep_psf =  19 #7 * 60. # arcmin
nbins_psf = 100

# Correlator options
npatch=300 #1 #300
var_method='jackknife'#'sample' #'marked_bootstrap' #'bootstrap' # #'jackknife' #'shot'


###################################### FILENAMES ##########################################################
# Input
data_dir = 'data/'

### Images ###
# Mosaic
filename_fits_big = data_dir + 'mosaic-image.fits'
filename_psf_big = data_dir + 'mosaic-psf.fits'
file_footprint_big = data_dir + 'mosaic-footprint.fits'
filename_fits_big_f0 = data_dir + 'field0_final.fits'
filename_fits_big_f1 = data_dir + 'field1_final.fits'

# Small images
filename_fits_small_coldest = data_dir + 'coldest-spot-residual-natural.fits'
filename_psf_small_coldest = data_dir + 'coldest-spot-psf-natural.fits'
filename_psf_profile_small_coldest_marco = data_dir + 'LoFAR_psf_coldest_marco.dat'

filename_fits_small_reference = data_dir + 'reference-residual-natural.fits'
filename_psf_small_reference = data_dir + 'reference-psf-natural.fits'
filename_psf_profile_small_reference_marco = data_dir + 'LoFAR_psf_reference_marco.dat'

# Template curves
template_dir = 'rexcorrconaltricataloghi/'


### Catalogs ###
catalog = 0 # 0: sdss12, 1: sdss16q, 2: unwise, 3: 2MPZ, 4: PSZ2


# Small images
if catalog == 0:
    catname = 'SDSS_DR12'
    filename_full_catalog = '../sdssdr12/dr12.dat'
    filename_full_catalog_formatted = data_dir + 'sdssdr12_big.dat'
    filename_small_catalog_formatted_coldest = data_dir + 'sdssdr12_small_coldest.dat'
    filename_small_catalog_formatted_reference = data_dir + 'sdssdr12_small_reference.dat'
    filename_small_catalog_formatted_coldest_no_extra = data_dir + 'sdssdr12_small_coldest_no_extra.dat'
    filename_small_catalog_formatted_reference_no_extra = data_dir + 'sdssdr12_small_reference_no_extra.dat'
    filename_zhist_full = data_dir + 'sdssdr12_zhist_full.dat'
    filename_zhist_coldest = data_dir + 'sdssdr12_zhist_coldest.dat'
    filename_zhist_reference = data_dir + 'sdssdr12_zhist_reference.dat'
    filename_template = template_dir + 'CCF_template_2h_SDSSmain.dat'
    dz = 0.01
elif catalog == 1:
    catname = 'SDSS_DR16Q'
    filename_full_catalog1 = '../sdssdr16q/DR16Q_v4.fits'
    filename_full_catalog2 = '../sdssdr16q/DR16Q_Superset_v3.fits'
    filename_full_catalog_formatted = data_dir + 'sdssdr16q_big.dat'
    filename_small_catalog_formatted_coldest = data_dir + 'sdssdr16q_small_coldest.dat'
    filename_small_catalog_formatted_reference = data_dir + 'sdssdr16q_small_reference.dat'
    filename_small_catalog_formatted_coldest_no_extra = data_dir + 'sdssdr16q_small_coldest_no_extra.dat'
    filename_small_catalog_formatted_reference_no_extra = data_dir + 'sdssdr16q_small_reference_no_extra.dat'
    filename_zhist_full = data_dir + 'sdssdr16q_zhist_full.dat'
    filename_zhist_coldest = data_dir + 'sdssdr16q_zhist_coldest.dat'
    filename_zhist_reference = data_dir + 'sdssdr16q_zhist_reference.dat'
    filename_template = template_dir + 'CCF_template_2h_SDSSqso.dat'
    dz = 0.1
elif catalog == 2:
    catname = 'UNWISE'
    filename_full_catalog1 = '../unwise/unwise.unwise_2019_29558.tbl.txt'
    filename_full_catalog2 = '../unwise/unwise.unwise_2019_24615.tbl.txt'
    filename_full_catalog3 = '../unwise/unwise.unwise_2019_5765.tbl.txt'
    filename_full_catalog_formatted = data_dir + 'unwise_big.dat'
    filename_small_catalog_formatted_coldest = data_dir + 'unwise_small_coldest.dat'
    filename_small_catalog_formatted_reference = data_dir + 'unwise_small_reference.dat'
    filename_small_catalog_formatted_coldest_no_extra = data_dir + 'unwise_small_coldest_no_extra.dat'
    filename_small_catalog_formatted_reference_no_extra = data_dir + 'unwise_small_reference_no_extra.dat'
    filename_template = template_dir + 'CCF_template_2h_unWISE.dat'
elif catalog == 3:
    catname = '2MPZ'
    filename_full_catalog = '../2MPZ/results1_11_45_29_2.csv'
    filename_full_catalog_formatted = data_dir + '2mpz_big.dat'
    filename_small_catalog_formatted_coldest = data_dir + '2mpz_small_coldest.dat'
    filename_small_catalog_formatted_reference = data_dir + '2mpz_small_reference.dat'
    filename_small_catalog_formatted_coldest_no_extra = data_dir + '2mpz_small_coldest_no_extra.dat'
    filename_small_catalog_formatted_reference_no_extra = data_dir + '2mpz_small_reference_no_extra.dat'
    filename_zhist_full = data_dir + '2mpz_zhist_full.dat'
    filename_zhist_coldest = data_dir + '2mpz_zhist_coldest.dat'
    filename_zhist_reference = data_dir + '2mpz_zhist_reference.dat'
    filename_template_1h = template_dir + 'CCF_template_1hext_2MPZ.dat'
    filename_template = template_dir + 'CCF_template_2h_2MPZ.dat'
    dz = 0.01
elif catalog == 4:
    catname = 'PSZ2'
    filename_full_catalog = '../PSZ2/HFI_PCCS_SZ-union_R2.08.fits'
    filename_full_catalog_formatted = data_dir + 'psz2_big.dat'
    filename_small_catalog_formatted_coldest = data_dir + 'psz2_small_coldest.dat'
    filename_small_catalog_formatted_reference = data_dir + 'psz2_small_reference.dat'
    filename_zhist_full = data_dir + 'psz2_zhist_full.dat'
    dz = 0.01

### Images converted to catalogs ##
# Small images
filename_fits_catalog_small_coldest = filename_fits_small_coldest[:-4] + 'dat'
filename_fits_catalog_small_reference = filename_fits_small_reference[:-4] + 'dat'

# Mosaic
filename_fits_catalog_big = filename_fits_big[:-4] + 'dat'
filename_fits_catalog_big_f0 = filename_fits_big_f0[:-4] + 'dat'
filename_fits_catalog_big_f1 = filename_fits_big_f1[:-4] + 'dat'




### Output
if catalog == 0:
    corr_dir = 'corrs_sdssdr12/'
    plot_dir = 'plots_sdssdr12/'
elif catalog == 1:
    corr_dir = 'corrs_sdssdr16q/'
    plot_dir = 'plots_sdssdr16q/'
elif catalog == 2:
    corr_dir = 'corrs_unwise/'
    plot_dir = 'plots_unwise/'
elif catalog == 3:
    corr_dir = 'corrs_2mpz/'
    plot_dir = 'plots_2mpz/'
elif catalog == 4:
    corr_dir = 'corrs_psz2/'
    plot_dir = 'plots_psz2/'
    
# Cross
filename_cross_small_coldest = corr_dir + catname + '_cross_coldest_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_cov_cross_small_coldest = corr_dir + catname + '_cov_cross_coldest_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_plot_cross_small_coldest = plot_dir + catname + '_xcorr_coldest_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_psf_profile_small_coldest = data_dir + 'psf_profile_small_coldest_min_sep=%.2f_max_sep=%.2f_nbins=%g.npz'

filename_cross_small_reference = corr_dir + catname + '_cross_reference_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_cov_cross_small_reference = corr_dir + catname + '_cov_cross_reference_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_plot_cross_small_reference = plot_dir + catname + '_xcorr_reference_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_psf_profile_small_reference = data_dir + 'psf_profile_small_reference_min_sep=%.2f_max_sep=%.2f_nbins=%g.npz'


filename_cross_big = corr_dir + catname + '_cross_mosaic_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_cov_cross_big = corr_dir + catname + '_cov_mosaic_cross_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_plot_cross_big = plot_dir + catname + '_xcorr_mosaic_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_psf_profile_big = data_dir + 'psf_profile_mosaic_min_sep=%.2f_max_sep=%.2f_nbins=%g.npz'

filename_cross_big_f0 = corr_dir + catname + '_cross_mosaic_f0_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_plot_cross_big_f0 = plot_dir + catname + '_xcorr_mosaic_f0_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_cross_big_f1 = corr_dir + catname + '_cross_mosaic_f1_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_plot_cross_big_f1 = plot_dir + catname + '_xcorr_mosaic_f1_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'

filename_plot_coldest_reference = plot_dir + catname + '_xcorr_compare_with_errors_sdssdr16q.png'


# Auto
filename_auto_coldest = corr_dir + 'auto_coldest_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_cov_auto_coldest = corr_dir + 'cov_auto_coldest_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_plot_auto_coldest = plot_dir + 'corr_auto_coldest_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_Cls_auto_coldest = corr_dir + 'auto_coldest_Cls_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_Cls_psf_coldest = corr_dir + 'auto_coldest_Cls_psf_min_sep=%.2f_max_sep=%.2f_nbins=%g.npz'

filename_auto_reference = corr_dir + 'auto_reference_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_cov_auto_reference = corr_dir + 'cov_auto_reference_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_plot_auto_reference = plot_dir + 'corr_auto_reference_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_Cls_auto_reference = corr_dir + 'auto_reference_Cls_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_Cls_psf_reference = corr_dir + 'auto_reference_Cls_psf_min_sep=%.2f_max_sep=%.2f_nbins=%g.npz'

filename_plot_Cls_psf_compare = plot_dir + 'PSF_small_Cls_compare_min_sep=%.2f_max_sep=%.2f_nbins=%g.png'

filename_auto_big = corr_dir + 'auto_mosaic_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_cov_auto_big = corr_dir + 'cov_auto_mosaic_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'
filename_auto_big_f0 = corr_dir + 'auto_mosaic_f0_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'
filename_auto_big_f1 = corr_dir + 'auto_mosaic_f1_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.dat'

filename_plot_auto_big = plot_dir + 'corr_auto_mosaic_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_plot_auto_big_f0 = plot_dir + 'corr_auto_mosaic_f0_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_plot_auto_big_f1 = plot_dir + 'corr_auto_mosaic_f1_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.png'
filename_Cls_psf_big = corr_dir + 'auto_mosaic_Cls_psf_min_sep=%.2f_max_sep=%.2f_nbins=%g.npz'
filename_Cls_auto_big = corr_dir + 'auto_mosaic_Cls_min_sep=%.2f_max_sep=%.2f_nbins=%g_var_method=%s_npatch=%d.npz'