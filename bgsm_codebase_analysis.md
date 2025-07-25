Project Path: Bayesian-Global-Sky-Model-B-GSM-Paper-1

Source Tree:

```txt
Bayesian-Global-Sky-Model-B-GSM-Paper-1
├── gen_EDGES_beams.py
├── mk24_mixed_models_v4_mock_data_gen_comp_amps_final.py
├── mock_data_v4_mk24_mixed_models_nested_sampling_perfect_noise.py
├── pix_by_pix_mk24_final2.py
├── plot_mean_comps_and_skys_v2.py
├── plot_posterior_results_for_paper_mixed_model_with_crosshair_v2.py
└── test_approximation2.py

```

`Bayesian-Global-Sky-Model-B-GSM-Paper-1/gen_EDGES_beams.py`:

```py
   1 | import matplotlib.pyplot as plt
   2 | import numpy as np
   3 | from scipy.optimize import minimize
   4 | from scipy.interpolate import CubicSpline
   5 | import os
   6 | import healpy as hp
   7 | 
   8 | 
   9 | def gen_EDGES_beam(freq,Max_Nside,azimuth=0,lat_FWHM=None,long_FWHM=None,beam_gain=None):
  10 | 
  11 |     Npix = hp.nside2npix(Max_Nside)
  12 |     pixel_thetas,pixel_phis = hp.pix2ang(Max_Nside, np.arange(Npix))
  13 | 
  14 |     #mask out below the EDGES obs horizon
  15 |     obs_lat = -26.7
  16 |     rot3=hp.Rotator(rot=(0,90-obs_lat))
  17 |     mask=np.zeros(Npix)
  18 |     mask[(pixel_thetas<=np.pi/2)] = 1
  19 |     mask=rot3.rotate_map_alms(mask)
  20 |     
  21 |     pixel_thetas,phis = hp.pix2ang(Max_Nside, np.arange(Npix),lonlat=True)
  22 | 
  23 |     #lat_FWHM = 112
  24 |     #long_FWHM = 72
  25 |     #beam_gain = 7
  26 | 
  27 |     if beam_gain==None:
  28 |         beam_gain = np.interp(freq,np.array([45,150]),np.array([7,5.89]))
  29 |     #lat_FWHM = np.interp(freq,np.array([45,150]),np.array([98,110]))
  30 |     #long_FWHM = np.interp(freq,np.array([45,150]),np.array([68,72]))
  31 |     if long_FWHM==None:
  32 |         long_FWHM = np.interp(freq,np.array([45,150]),np.array([98,112]))
  33 |     if lat_FWHM==None:
  34 |         lat_FWHM = np.interp(freq,np.array([45,150]),np.array([68,72]))
  35 | 
  36 |     #print ("for freq",freq,"gain",beam_gain,"lat FWHM",lat_FWHM,"long_FWHM",long_FWHM)
  37 | 
  38 |     beam = np.exp(-0.5*(phis/(lat_FWHM/2.355))**2)*np.exp(-0.5*((pixel_thetas-180)/(long_FWHM/2.355))**2)
  39 |     
  40 |      #center the peak of the beam to a longditude of 0 deg
  41 |     rot=hp.Rotator(rot=(180,0))
  42 |     beam = rot.rotate_map_alms(beam)
  43 |     #hp.mollview(beam,title="beam straight up")
  44 |     #plt.show()
  45 | 
  46 |     #rotate so the beams peak is at the north pole 
  47 |     rot=hp.Rotator(rot=(0,-90))
  48 |     beam = rot.rotate_map_alms(beam)
  49 |     #hp.mollview(beam,title="beam center at north pole")
  50 |     #plt.show()
  51 |     
  52 |     #rotate to the correct azimuth angle (rotate around the north pole)
  53 |     rot1=hp.Rotator(rot=(azimuth,0))
  54 |     beam = rot1.rotate_map_alms(beam)
  55 |     #hp.mollview(beam,title="rotated to azimuth="+str(azimuth))
  56 |     #plt.show()
  57 |     
  58 |     #rotate back to having the beam point straight up but at the correct azimuth angle
  59 |     rot=hp.Rotator(rot=(0,90))
  60 |     beam = rot.rotate_map_alms(beam)
  61 |     #hp.mollview(beam,title="beam straight up at azimuth angle"+str(azimuth))
  62 |     #plt.show()
  63 | 
  64 |     
  65 |     #rotate the beam to correct for our observation latitude 
  66 |     rot2=hp.Rotator(rot=(0,0-obs_lat))
  67 |     beam = rot2.rotate_map_alms(beam)
  68 | 
  69 |     beam *=beam_gain
  70 | 
  71 |     
  72 |     return beam, mask
  73 | 
  74 | 
  75 | def gen_EDGES_beams_at_LSTs(freq,LSTs,Nside,azimuth=0,test=False,lat_FWHM=None,long_FWHM=None,beam_gain=None):
  76 |     beam_at_LST0, horizon0 = gen_EDGES_beam(freq,Nside,azimuth=azimuth,lat_FWHM=lat_FWHM,long_FWHM=long_FWHM,beam_gain=beam_gain) #the beam for LST 0 in equatorial coords
  77 |     rot4 = hp.Rotator(coord=["C","G"])
  78 |     beams = []
  79 |     
  80 |     if test==True:
  81 |         hp.mollview(beam_at_LST0,title="beam at LST=0")
  82 |         plt.show()
  83 |         ts =[]
  84 |         m = hp.read_map(os.getcwd()+"/dataset_Haslam_LWA1_uncal_Guz_LW_Monsalve_cal_smoothed_FWHM=5_at_Nside=32/map_"+str(freq)+"MHz_FWHM=5_nside=32.fits")
  85 |         m[(m==-32768)]=float("NaN")
  86 |         for LST in LSTs:
  87 |             rotation = hp.Rotator(rot=(-15*LST,0))
  88 |             beam = rotation.rotate_map_alms(beam_at_LST0)
  89 |             horizon = rotation.rotate_map_alms(horizon0)
  90 |             #hp.mollview(beam,title="beam for LST="+str(LST))
  91 |             #plt.show()
  92 | 
  93 |             #rotate the beam into galactic coords to match the map, do the same to the horizon
  94 |             #rotate to galactic coords
  95 |             beam = rot4.rotate_map_alms(beam)
  96 |             horizon = rot4.rotate_map_alms(horizon)
  97 |             
  98 |             mask2use=np.zeros(len(beam))
  99 |             mask2use[(horizon<0.5)] = 0
 100 |             mask2use[(horizon>=0.5)] = 1
 101 |             #mask below the horizon
 102 |             beam *= mask2use
 103 |             beams.append(beam)
 104 |             #hp.mollview(beam,title="beam for LST="+str(LST)+" Galactic coords")
 105 |             #plt.show()
 106 | 
 107 |             m1=np.zeros(len(beam))
 108 |             m2=np.zeros(len(beam))
 109 | 
 110 |             #m1[(beam>0.495*np.max(beam))]=1
 111 |             #m2[(beam<0.505*np.max(beam))]=1
 112 |             #hp.mollview(m1*m2)
 113 |             #plt.show()
 114 |             #print (np.nansum((1/(4*np.pi))*beam*m*hp.nside2pixarea(Nside)))
 115 |             ts.append(np.nansum((1/(4*np.pi))*beam*m*hp.nside2pixarea(Nside)))
 116 |         return ts,beams
 117 |     
 118 |     else:
 119 |         beams = np.zeros((len(beam_at_LST0),len(LSTs)))
 120 |         i=0
 121 |         for LST in LSTs:
 122 |             rotation = hp.Rotator(rot=(-15*LST,0))
 123 |             beam = rotation.rotate_map_alms(beam_at_LST0)
 124 |             horizon = rotation.rotate_map_alms(horizon0)
 125 |             #hp.mollview(beam,title="beam for LST="+str(LST))
 126 |             #plt.show()
 127 | 
 128 |             #rotate the beam into galactic coords to match the map, do the same to the horizon
 129 |             #rotate to galactic coords
 130 |             beam = rot4.rotate_map_alms(beam)
 131 |             horizon = rot4.rotate_map_alms(horizon)
 132 |             
 133 |             mask2use=np.zeros(len(beam))
 134 |             mask2use[(horizon<0.5)] = 0
 135 |             mask2use[(horizon>=0.5)] = 1
 136 |             #mask below the horizon
 137 |             beam *= mask2use
 138 | 
 139 |             beams[:,i] = beam
 140 |             i+=1
 141 |         return beams
 142 | #lst = np.linspace(0,24,97)
 143 | #t,b=gen_EDGES_beams_at_LSTs(45,np.array(lst),32,azimuth=-5,test=True)#,lat_FWHM=71.6,long_FWHM=110)
 144 | #print ("simulated ant temps for 45MHz")
 145 | #print (t)
 146 | #zs,z2s=gen_EDGES_high_T_LST_trace(45,lst)
 147 | 
 148 | 
 149 | 
 150 | #fig = plt.figure(figsize=(10,10))
 151 | #plt.errorbar(lst,zs,z2s,label="gen from spec index")#
 152 | #plt.plot(lst,t,label="simulated obs using EDGES beam model")
 153 | #plt.legend(loc="upper left")
 154 | #plt.title("EDGES high 45 MHz")
 155 | #plt.show()
 156 | 
 157 | #t=np.array(t)
 158 | #zs=np.array(zs)
 159 | #func = lambda x: np.sum(((zs-x[0]*t-x[1])/np.array(z2s))**2)
 160 | #res=minimize(func,[1,0])
 161 | #print (res)
 162 | 
 163 | #print ("reduced chi sqr after corrections:",res.fun/(len(lst)-2))
 164 | 
 165 | #fig = plt.figure(figsize=(10,10))
 166 | #plt.errorbar(lst,zs,z2s,label="gen from spec index")#
 167 | #plt.plot(lst,res.x[0]*t+res.x[1],label="simulated obs using EDGES beam model")
 168 | #plt.legend(loc="upper left")
 169 | #plt.title("EDGES high 45 MHz")
 170 | #plt.show()
 171 | 
 172 | 

```

`Bayesian-Global-Sky-Model-B-GSM-Paper-1/mk24_mixed_models_v4_mock_data_gen_comp_amps_final.py`:

```py
   1 | from random import sample
   2 | import matplotlib as mpl
   3 | #mpl.use('Agg')
   4 | from itertools import combinations_with_replacement
   5 | from itertools import product
   6 | import numpy as np
   7 | import healpy as hp
   8 | import matplotlib.pyplot as plt
   9 | from mpl_toolkits.axes_grid1 import make_axes_locatable
  10 | from matplotlib.colors import SymLogNorm
  11 | from matplotlib.colors import LogNorm
  12 | import os
  13 | import pandas
  14 | from numpy import pi, log, sqrt
  15 | import pix_by_pix_mk24_final2 as likelihood
  16 | 
  17 | import scipy.optimize as so
  18 | import math
  19 | import scipy
  20 | 
  21 | import matplotlib.cm as cm
  22 | from scipy.optimize import minimize
  23 | 
  24 | try:
  25 |     import pypolychord
  26 |     from pypolychord.settings import PolyChordSettings
  27 |     from pypolychord.priors import UniformPrior
  28 |     from pypolychord.priors import GaussianPrior
  29 | except:
  30 |     pass
  31 | try:
  32 |     from anesthetic import NestedSamples
  33 | except ImportError:
  34 |     pass
  35 | try:
  36 |     from anesthetic.weighted_pandas import WeightedDataFrame
  37 | except ImportError:
  38 |     pass
  39 | from scipy.optimize import minimize
  40 | import gen_EDGES_beams as gen_beams_and_T_vs_LST_v2
  41 | from line_profiler import LineProfiler
  42 | import matplotlib
  43 | from astropy.modeling.powerlaws import LogParabola1D
  44 | from csv import writer
  45 | 
  46 | #RUN PARAMS
  47 | #================================================================================================
  48 | chunks = [[0,10000],[10000,20000],[20000,30000],[30000,40000],[40000,50000],[50000,60000],[60000,70000],[70000,80000],[80000,90000]]#[[0,5000],[5000,10000],[10000,15000],[15000,20000],[20000,25000],[25000,30000],[30000,35000],[35000,40000],[40000,45000],[45000,50000],[50000,55000],[55000,60000],[60000,65000],[65000,70000],[70000,75000]]
  49 | 
  50 | chunks_to_use_for_run = [0,1,2]
  51 | #chunks_to_use_for_run = [6,7,8]
  52 | #chunks_to_use_for_run = [4,5,7]
  53 | #chunks_to_use_for_run = [3,4,5]
  54 | #chunks_to_use_for_run = [6,7]#,8]
  55 | 
  56 | #MODEL PARAMS
  57 | #================================================================================================
  58 | #declare a random seed
  59 | np.random.seed(0)
  60 | use_perturbed_dataset = True #do we want the input dataset to have calibration errors
  61 | 
  62 | #===================================================================#
  63 | #| SET PARAMS FOR THE SIMULATED DATA AND NESTED SAMPLING
  64 | Max_Nside=32 #the Nside at which to generate the set of maps
  65 | Max_m = hp.nside2npix(Max_Nside)
  66 | no_of_comps = 2
  67 | fit_curved=[True,True]
  68 | 
  69 | no_to_fit_curve = np.sum(fit_curved)
  70 | rezero_prior_std=2000
  71 | 
  72 | #params for the prior on the spectra
  73 | spec_min, spec_max = -3.5,1 #the range for the prior on the spectral indexes
  74 | curvature_mean, curvature_std = 0,3 #the range for the prior on the spectral index curvature parameter
  75 | 
  76 | #params for fitting the reference frequency 
  77 | 
  78 | fixed_f0 = 150 #if you dont fit a seperate reference freq for each comp then we fix f0 to this value
  79 | 
  80 | 
  81 | #params for the prior on the true maps
  82 | map_prior_variance_spec_index = -2.6
  83 | map_prior_variance_f0 = 408#fixed_f0
  84 | map_prior_std=300
  85 | 
  86 | calibrate = True
  87 | calibrate_all_but_45_150 = False#True #calibrate all maps in the dataset but the 45 and 150 MHz maps
  88 | calibrate_all = True#False #calibrate every map in the dataset
  89 | 
  90 | 
  91 | use_equal_spaced_LSTs = True
  92 | fit_haslam_noise = False
  93 | subtract_CMB = 0#-2.726
  94 | print_vals_as_calc = False
  95 | 
  96 | 
  97 | 
  98 | reject_criterion = 1e-3#None #how close can two spectral idexes be in value before being rejected
  99 | cond_no_threshold =1e+9
 100 | 
 101 | 
 102 | 
 103 | 
 104 | 
 105 | unobs_marker = -32768
 106 | 
 107 | nlive = 500#2500*no_of_comps
 108 | 
 109 | precision_criterion = 1e-3
 110 | 
 111 | 
 112 | f0=150 #ref freq used for some fitting of spec indexes for plots (not used during any model fitting)
 113 | 
 114 | freqs = np.array([45.0,50.0,60.0,70.0,74.0,80.0,150.0,159.0,408.0]) #the frequencies in MHz of maps used to generate the model
 115 | #specify what to fit for each map
 116 | if calibrate_all == True:
 117 |     freqs_to_calibrate = np.array([True,True,True,True,True,True,True,True,True]) #calibrate all the maps except the 45 and 150 MHz
 118 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 119 | 
 120 | if calibrate_all_but_45_150 == True:
 121 |     freqs_to_calibrate = np.array([False,True,True,True,True,True,False,True,True]) #calibrate all the maps except the 45 and 150 MHz
 122 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 123 | if calibrate==False:
 124 |     freqs_to_calibrate = np.array([False,False,False,False,False,False,False,False,False])#np.array([True,True,True,True,True,True,False,False])
 125 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 126 | 
 127 | 
 128 | 
 129 | 
 130 | main_label = "_petur:"+str(use_perturbed_dataset)+"_"+str(no_of_comps)+"_comp_cal:"+str(calibrate)+"_rezro_pri_std:"+str(rezero_prior_std)+"_CMB="+str(subtract_CMB)+"_map_pri_std:"+str(map_prior_std)+"_mu:0_map_pri_std_spec_ind="+str(map_prior_variance_spec_index)+"_map_pri_f0="+str(map_prior_variance_f0)+"_cond_no_thres="+str(np.round(np.log10(cond_no_threshold),1))+"_crv_N_std="+str(curvature_std)+"_spec="+str(spec_min)+"_to:"+str(spec_max)#+"_rej_crit="+str(reject_criterion)#+"_nlive="+str(nlive)+"_nrept="+str(nrepeat)+"_precision_criterion="+str(precision_criterion)
 131 | 
 132 | if use_equal_spaced_LSTs==True:
 133 |     #LSTs_for_comparison = np.array([2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])#np.array([0,2,4,6,8,10,12,14,16,18,20,22])#np.array([2.5,18]) #the LSTs in hours at which we will make comparison between the mean sky and EDGES for likelihood calls
 134 |     LSTs_for_comparison = np.linspace(0,24,73)[:-1]
 135 |     print (LSTs_for_comparison)
 136 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 137 |     if calibrate_all==True:
 138 |         #root="uni_EDGES_v4_data_mk24_no_of_curved:"+str(no_to_fit_curve)+"_cal_all_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 139 |         root="uni_EDGES_v4_data_mk24_no_of_curved:"+str(no_to_fit_curve)+"_cal_all_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 140 |     
 141 |     else:
 142 |         if calibrate_all_but_45_150==True:
 143 |             root="uni_EDGES_v4_dat_mk24_no_curve:"+str(no_to_fit_curve)+"_no_cal_45_150_f0="+str(fixed_f0)+main_label
 144 | else:
 145 |     #LSTs_for_comparison = np.array([0,1,2,3,4,5,6,7,8,9,10,11,11.25,11.5,11.75,12,12.25,12.5,12.75,13,13.25,13.5,13.75,14,14.25,14.5,14.75,15,15.25,15.5,15.75,16,16.25,16.5,16.75,17,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18,18.25,18.5,18.75,19,19.25,19.5,19.75,20,20.25,20.5,20.75,21,21.25,21.5,21.75,22,22.25,22.5,22.75,23,23.25,23.5,23.75])
 146 |     LSTs_for_comparison = np.array([0,2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])
 147 |     print (LSTs_for_comparison)
 148 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 149 |     
 150 |     root="mk24_extra_uneq_LSTs:"+str(len(LSTs_for_comparison))+"_fixed_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 151 |     
 152 | 
 153 | #specify what to fit for each map
 154 | 
 155 | 
 156 | #set the prior for the noise (on the Haslam map)
 157 | noise_prior_lower, noise_prior_upper = np.array([0.01]),np.array([50])
 158 | 
 159 | print ("noise prior is from:",noise_prior_lower,"to",noise_prior_upper,"Kelvin")
 160 | 
 161 | n_spec_pars = 3*no_of_comps #the number of parameters for the spectra (for each comp we have: break_freq, spec_index1, spec_index2)
 162 | 
 163 | 
 164 | nv=len(freqs) #the number of freqs that have maps
 165 | no_of_fitted_noise = np.sum(freqs_to_fit_noise) #the number of freqs at which we fit noise level
 166 | no_of_calibrated = np.sum(freqs_to_calibrate) #the number of freqs at which we fit the calibration
 167 | 
 168 | if np.sum(freqs_to_calibrate)!=0:
 169 |     zero_lev_prior_std = rezero_prior_std*np.ones(np.sum(freqs_to_calibrate))#200*((np.array(freqs)[freqs_to_calibrate]/100)**-2.5)
 170 |     zero_lev_prior_means = np.zeros(np.sum(freqs_to_calibrate))
 171 | 
 172 | 
 173 |     #set the prior params for the scale corrections
 174 |     scale_prior_lower=0.85
 175 |     scale_prior_upper=1.25
 176 |     
 177 | 
 178 |     print ("zero level prior is gauss with mean 0K, stds (Kelvin):")
 179 |     print (zero_lev_prior_std)
 180 |     print ("temp scale prior is uniform from:",scale_prior_lower,"to",scale_prior_upper)
 181 | 
 182 | 
 183 | 
 184 | 
 185 | 
 186 | 
 187 | #CREATE A DIR TO STORE RESULTS
 188 | #====================================================================#
 189 | 
 190 | 
 191 | #make a dir to store the results
 192 | p=os.getcwd()+"/"
 193 | path = p+root+"/"
 194 | try:
 195 |     os.mkdir(root)
 196 | except:
 197 |     pass
 198 | #make a dir to store the results as we run
 199 | root2 = path+"running_results/"
 200 | try:
 201 |     os.mkdir(root2)
 202 | except:
 203 |     pass
 204 | 
 205 | #LOAD THE DATASET AND THE ERROR MAPS
 206 | #====================================================================#
 207 | 
 208 | obs_maps = []
 209 | inv_err_maps = []
 210 | data_err_maps = []
 211 | load_path = p+"mock_data_file/mock_dataset_v4/"
 212 | for i in range(len(freqs)):
 213 |     f=freqs[i]
 214 | 
 215 |     if use_perturbed_dataset==True:
 216 |         fname1 = "noisy_perturbed_sky_"+str(f)
 217 |         err_fname = "perturbed_err_map_"+str(f)
 218 |     else:
 219 |         fname1 = "noisy_sky_"+str(f)
 220 |         err_fname = "err_map_"+str(f)
 221 |     
 222 |     #fname2 = "noise_"+str(f)
 223 | 
 224 |     if freqs_to_fit_noise[i]==False:
 225 |         try:
 226 |             #m1, err_m = np.loadtxt(load_path+fname1), np.loadtxt(load_path+fname2)
 227 |             m1 = np.loadtxt(load_path+fname1)
 228 | 
 229 |             err_m = np.loadtxt(load_path+err_fname)
 230 |         except:
 231 |             print ("cant find the files for freq:",f)
 232 | 
 233 |         
 234 |         #mask out any pixels with negative temps
 235 |         bool_arr = m1<=0
 236 |         err_m[bool_arr] = unobs_marker
 237 |         m1[bool_arr] = unobs_marker
 238 | 
 239 |         
 240 | 
 241 |         inv_err_m = 1/err_m
 242 |         inv_err_m[(err_m==unobs_marker)] = 0
 243 | 
 244 |         m1[m1!=unobs_marker] = m1[m1!=unobs_marker]+subtract_CMB
 245 |         obs_maps.append(m1)
 246 |     
 247 |         inv_err_maps.append(inv_err_m)
 248 | 
 249 |         data_err_maps.append(err_m)
 250 |     else:
 251 |         m1 = np.loadtxt(load_path+fname1)
 252 |         m1[m1!=unobs_marker] = m1[m1!=unobs_marker]+subtract_CMB
 253 |         obs_maps.append(m1)
 254 |     
 255 | 
 256 | obs_maps=np.array(obs_maps)
 257 | inv_err_maps=np.array(inv_err_maps)
 258 | print ("dataset loaded")
 259 | 
 260 | #CREATE THE INVERSE NOISE MATRICES
 261 | #====================================================================#
 262 | #generate the inverse noise covariance matrix for each pixel
 263 | inv_noise_mats = np.empty(shape=(Max_m,len(freqs),len(freqs)))
 264 | for p in range(0,Max_m):
 265 |     inv_stds_for_pixel = np.zeros(len(freqs))
 266 |     inv_stds_for_pixel[~freqs_to_fit_noise] = inv_err_maps[:,p]
 267 |     #print (inv_stds_for_pixel)
 268 | 
 269 |     Np_inv = np.diag(inv_stds_for_pixel**2)
 270 |     #print (Np_inv)
 271 |     inv_noise_mats[p,:,:] = Np_inv
 272 | 
 273 | print ("max and min for the inv noise mats: ",np.max(inv_noise_mats),np.min(inv_noise_mats[(inv_noise_mats!=0)]))
 274 | 
 275 | print ("inverse noise matrices created")
 276 | #PLOT THE DATASET
 277 | #====================================================================#
 278 | fig = plt.figure(figsize=(12,16))
 279 | for i in range(len(freqs)):
 280 |     map_i = np.copy(obs_maps[i])
 281 |     ax = plt.subplot(5,3,int(i+1))
 282 |     map_i[(map_i==unobs_marker)]=float("NaN")
 283 |     plt.axes(ax)
 284 |     hp.mollview(map_i,title="Synthetic Data Freq="+str(freqs[i]),hold=True,notext=True,norm="log")
 285 | 
 286 | 
 287 | 
 288 | plt.savefig(path+"sky_maps_for_dataset_for_plt")
 289 | #plt.show()
 290 | fig = plt.figure(figsize=(12,16))
 291 | for i in range(len(freqs)):
 292 |     map_i = np.copy(inv_err_maps[i])
 293 |     ax = plt.subplot(5,3,int(i+1))
 294 |     map_i[(map_i==0)]=float("NaN")
 295 |     plt.axes(ax)
 296 |     hp.mollview(1/map_i,title="input errs freq="+str(freqs[i]),hold=True,notext=True,norm="log")
 297 | 
 298 | 
 299 | 
 300 | plt.savefig(path+"err_maps_for_dataset_for_plt")
 301 | freqs=np.array(freqs)
 302 | #plot the priors and the data
 303 | #====================================================================
 304 | log_mean_temps = []
 305 | mean_temps = []
 306 | for i in range(len(freqs)):
 307 |     the_map = obs_maps[i]
 308 |     mean = np.mean(the_map[(the_map!=unobs_marker)])
 309 |     log_mean_temps.append(np.log(mean))
 310 |     mean_temps.append(mean)
 311 | log_mean_temps = np.array(log_mean_temps)
 312 | mean_temps = np.array(mean_temps)
 313 | 
 314 | log_freqs = np.log(freqs/f0)
 315 | fun = lambda x: np.nansum((log_mean_temps - x[0]*log_freqs -x[1])**2)
 316 | res = minimize(fun,[-2.15,np.log(16)])
 317 | print (res)
 318 | 
 319 | fig = plt.figure(figsize=(6,6))#figsize=(12,16))
 320 | #the fitted powerlaw
 321 | targ = np.exp(res.x[1])*((freqs/f0)**res.x[0])
 322 | 
 323 | ax1 = plt.subplot(1,1,1)
 324 | ax1.plot(freqs,targ,c="red",label="fitted power law")
 325 | ax1.scatter(freqs,mean_temps,label="data")
 326 | ax1.set_title("map mean temps")
 327 | ax1.legend(loc="upper right")
 328 | ax1.set_xscale("log")
 329 | ax1.set_yscale("log")
 330 | plt.savefig(path+"/input_map_for_plt_mean_temps.png")
 331 | 
 332 | 
 333 | #NOTE: we set the terms relating to the EDGES beams to none. The component amplitude generation dosen't requier these terms
 334 | #we still have to pass None as the pix_by_pix_mk24_final2 code needs us to pass something for these terms
 335 | freqs_for_T_v_LST_comp = None
 336 | #generate a set of pre rotated EDGES beams at each of the frequencies that we want to compare the model to EDGES for
 337 | #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 338 | print ("generating EDGES beams")
 339 | 
 340 | 
 341 | EDGES_beams = None
 342 | #Load the mock EDGES T vs LST plots (these are produced by convolving the mock sky map (before any pertubation) with the model beam at that freq)
 343 | #===============================================================================================
 344 | 
 345 | 
 346 | EDGES_temps_at_calib_LSTs_and_freqs = None#np.array(EDGES_temps_at_calib_LSTs_and_freqs)
 347 | EDGES_errs = None#np.array(EDGES_errs)
 348 | #generate the EDGES noise covar mats for each freq assuming noise for each LST is independent of the other LSTs
 349 | 
 350 | EDGES_inv_noise_mats = None#np.array(EDGES_inv_noise_mats)
 351 | #=====================================================================================================
 352 | #set up the likelihood function
 353 | #bayes_eval = likelihood.bayes_mod(obs_maps=obs_maps,obs_freqs=freqs,inv_noise_mats=inv_noise_mats,gaussian_prior_covar_mat=gaussian_prior_covar_matrix,gaussian_prior_mean=gaussian_prior_mean,no_of_comps=no_of_comps,f0=f0,un_obs_marker=unobs_marker)
 354 | 
 355 | #set up the likelihood function
 356 | bayes_eval = likelihood.bayes_mod(obs_maps=obs_maps,obs_freqs=freqs,inv_noise_mats=inv_noise_mats,EDGES_beams=EDGES_beams,EDGES_temps_at_calib_LSTs_and_freqs=EDGES_temps_at_calib_LSTs_and_freqs,EDGES_errs=EDGES_errs,EDGES_inv_noise_mats=EDGES_inv_noise_mats,freqs_for_T_v_LST_comp=None,LSTs_for_comparison=None,no_of_comps=no_of_comps,save_root=root2,un_obs_marker=unobs_marker,map_prior_std=map_prior_std,map_prior_spec_index=map_prior_variance_spec_index,map_prior_f0=map_prior_variance_f0)
 357 | 
 358 | #test
 359 | #print ("===================================")
 360 | #print ("testing")
 361 | #sample_map = bayes_eval.gen_comp_map_sample([150,-2.50,0,150,-0.5,-2.5,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1],freqs_to_fit_noise,freqs_to_calibrate)
 362 | #for c in range(no_of_comps):
 363 | #    m = sample_map[:,c,0]
 364 | #    hp.mollview(m)
 365 | #    plt.savefig("test_c="+str(c)+".png")
 366 | 
 367 | 
 368 | no_of_params_for_spec_mod = no_of_comps + no_to_fit_curve 
 369 | print ("we arn't fitting f0: f0=",fixed_f0," no of comps with curved spectra is",no_to_fit_curve," no of params for spectral model is",no_of_params_for_spec_mod)
 370 | 
 371 | 
 372 | #-------------NESTED SAMPLING PARAMS-------------
 373 | nDims =  int(no_of_params_for_spec_mod + no_of_fitted_noise + 2*no_of_calibrated)
 374 | nrepeat = 5*nDims #the nrepeat is set to 5 times the total number of pars that we fit
 375 | print ("no of dimensions for sampling region is:",nDims)
 376 | nDerived = 0 #we don't derive any parameters 
 377 | settings = PolyChordSettings(nDims, nDerived)
 378 | settings.file_root = root
 379 | settings.nlive = nlive
 380 | settings.nrepeats = nrepeat
 381 | settings.do_clustering = True
 382 | settings.read_resume = True
 383 | settings.write_resume = True
 384 | settings.maximise = False #find the maximum of the poseterior
 385 | settings.precision_criterion = precision_criterion
 386 | #------------------------------------------------
 387 | 
 388 | #find the mean spectral params and their standard deviations
 389 | samples = NestedSamples(root= settings.base_dir + '/' + settings.file_root)
 390 | 
 391 | mean_logZ = samples.logZ()#mean
 392 | std_logZ = samples.logZ(100).std()#100 posterior samples from estimate of log Z (use these to calculate standard deviation in Z)
 393 | 
 394 | print ("log Z",mean_logZ,"std",std_logZ)
 395 | 
 396 | 
 397 | column_names_in_dataframe = list(samples.columns.values)
 398 | print ("column names for dataframe:", column_names_in_dataframe)
 399 | nested_samples = []
 400 | for i in range(nDims):
 401 |     print (i,column_names_in_dataframe[i])
 402 |     nested_samples.append(list(samples.loc[:,column_names_in_dataframe[i]]))
 403 | 
 404 |     
 405 | nested_samples = np.array(nested_samples).T#samples.loc[:,par_llamo].to_numpy()
 406 | weights = samples.weight
 407 | samps_with_weights = WeightedDataFrame(nested_samples,weight=weights)
 408 | mean_params = samps_with_weights.mean().to_numpy()
 409 | std_params = samps_with_weights.std().to_numpy()
 410 | print ("posterior mean for specs:")
 411 | print (mean_params)
 412 | print ("std")
 413 | print (std_params)
 414 | 
 415 | 
 416 | spec_pars_mean = mean_params[:no_of_params_for_spec_mod]
 417 | spec_pars_std = std_params[:no_of_params_for_spec_mod]
 418 | everything_else_mean = mean_params[no_of_params_for_spec_mod:]
 419 | everything_else_std = std_params[no_of_params_for_spec_mod:]
 420 | 
 421 | 
 422 | scale_corrections = np.ones(nv)
 423 | scale_correction_errs = np.zeros(nv)
 424 | scale_corrections[freqs_to_calibrate] = everything_else_mean[np.sum(freqs_to_fit_noise)+np.sum(freqs_to_calibrate):]
 425 | scale_correction_errs[freqs_to_calibrate] = everything_else_std[np.sum(freqs_to_fit_noise)+np.sum(freqs_to_calibrate):]
 426 | zero_corrections = np.zeros(nv)
 427 | zero_correction_errs = np.zeros(nv)
 428 | zero_corrections[freqs_to_calibrate] = everything_else_mean[np.sum(freqs_to_fit_noise):np.sum(freqs_to_fit_noise)+np.sum(freqs_to_calibrate)]
 429 | zero_correction_errs[freqs_to_calibrate] = everything_else_std[np.sum(freqs_to_fit_noise):np.sum(freqs_to_fit_noise)+np.sum(freqs_to_calibrate)]
 430 | 
 431 | print ("params for calibration:")
 432 | print (scale_corrections)
 433 | print (zero_corrections)
 434 | no_of_samples=nested_samples.shape[0]
 435 | print ("no of samples drawn from posterior is:",no_of_samples)
 436 | 
 437 | fig, axes = samples.plot_2d(['p%i' %i for i in range(1,nDims+1)])
 438 | 
 439 | fig.set_size_inches(16, 16)
 440 | 
 441 | matplotlib.rc('xtick', labelsize=25) 
 442 | matplotlib.rc('ytick', labelsize=25) 
 443 | fig.savefig(path+'sampled_posterior_update.png')
 444 | 
 445 | #write a text file containing the key results
 446 | with open(path+"post_run_sumary_stats.txt","w") as f:
 447 |     f.write("\n")
 448 |     f.write("a sumary file for the results of the nested sampling run.\n")
 449 |     f.write("freqs with maps in input dataset \n")
 450 |     f.write(str(freqs)+"\n")
 451 |     f.write("=========================================================\n")
 452 |     f.write("log(Z)="+str(mean_logZ)+" +/- "+str(std_logZ)+"\n")
 453 |     f.write("fixed f0 set as: "+str(fixed_f0)+"\n")
 454 |     f.write("spec indexes\n")
 455 |     f.write(str(spec_pars_mean)+"\n")
 456 |     f.write("+/-\n")
 457 |     f.write(str(spec_pars_std)+"\n")
 458 |     f.write("zero level corrections\n")
 459 |     f.write(str(zero_corrections)+"\n")
 460 |     f.write("+/-\n")
 461 |     f.write(str(zero_correction_errs)+"\n")
 462 |     f.write("temp scale correction factors\n")
 463 |     f.write(str(scale_corrections)+"\n")
 464 |     f.write("+/-\n")
 465 |     f.write(str(scale_correction_errs)+"\n")
 466 | 
 467 |     f.close()
 468 | 
 469 | #create a bool array of which params need updates
 470 | #select the spec indexes 
 471 | print ("=======================================================")
 472 | print ("Seting up spectral params selecter")
 473 | if no_to_fit_curve==no_of_comps:
 474 |     print ("all comps are curved spec")
 475 |     spec_indexes_select = np.tile(np.array([False,True,False]),no_of_comps)
 476 |     spec_curvature_select =  np.tile(np.array([False,False,True]),no_of_comps)
 477 |     spec_f0_select =  np.tile(np.array([True,False,False]),no_of_comps)
 478 |     spec_not_curve_select = np.tile(np.array([False,False,False]),no_of_comps)
 479 |     print ("spec_curvature_select =",spec_curvature_select)
 480 | else:
 481 |     print ("not all comps are curved spec")
 482 |     spec_indexes_select = np.tile(np.array([False,True,False]),no_of_comps)
 483 |     spec_f0_select =  np.tile(np.array([True,False,False]),no_of_comps)
 484 |     spec_curvature_select = np.array([np.array([False,False,True])*fit_curve_for_comp for fit_curve_for_comp in fit_curved]).flatten() #the indexes for spectral curvature in the final param array
 485 |     print ("spec_curvature_select =",spec_curvature_select)
 486 |     spec_not_curve_select = np.array([np.array([False,False,not_fit_curve_for_comp]) for not_fit_curve_for_comp in ~np.array(fit_curved)]).flatten() #the indexes with no curvature in the final param array
 487 |     print ("spec_not_curve select =",spec_not_curve_select)
 488 | 
 489 | #| Use the samples of the spectras to generate samples of component maps
 490 | ##########################################################################
 491 | import gzip
 492 | class generate_component_map_and_sky_map_samples():
 493 |     """this loads the nested samples and generates the associated set of sample component maps and sample sky maps"""
 494 |     
 495 |     def __init__(self):
 496 |         #load the samples of the poseterior for the spectra
 497 |         samples = NestedSamples(root= settings.base_dir + '/' + settings.file_root)
 498 | 
 499 |     
 500 | 
 501 |         column_names_in_dataframe = list(samples.columns.values)
 502 |         print ("column names for dataframe:", column_names_in_dataframe)
 503 |         nested_samples = []
 504 |         
 505 |         for i in range(nDims):
 506 |             print (i,column_names_in_dataframe[i])
 507 |             nested_samples.append(list(samples.loc[:,column_names_in_dataframe[i]]))
 508 |     
 509 |         nested_samples = np.array(nested_samples).T#samples.loc[:,par_llamo].to_numpy()
 510 |         weights = samples.weight
 511 | 
 512 |         #put the params into the correct positions
 513 |         ns_samps_arr = np.empty((len(weights),int(3*no_of_comps + no_of_fitted_noise + 2*no_of_calibrated)))
 514 |         #put all the calibration params into the sample array
 515 |         ns_samps_arr[:,n_spec_pars:] = nested_samples[:,no_of_params_for_spec_mod:]
 516 |         spec_pars_in_correct_order = np.empty((len(weights),int(3*no_of_comps)))
 517 |         #put the spectral slope params into the correct places:
 518 |         spec_pars_in_correct_order[:,spec_indexes_select] = nested_samples[:,:no_of_comps] #all compontents have a spec slope param these are the first no_of_comp params
 519 |         #put the spectral curvature params into the correct places
 520 |         spec_pars_in_correct_order[:,spec_curvature_select] = nested_samples[:,no_of_comps:no_of_params_for_spec_mod] #the curvature params for any with fitted curvature
 521 |         #for any comps that we didn't fit the curvature, we set the curvature to 0
 522 |         spec_pars_in_correct_order[:,spec_not_curve_select] = 0
 523 |         #put the reference freqs into the correct posititions    
 524 |         spec_pars_in_correct_order[:,spec_f0_select] = fixed_f0*np.ones((len(weights),no_of_comps))
 525 |         #finaly we put the spec pars into the samples array
 526 |         ns_samps_arr[:,:n_spec_pars] = spec_pars_in_correct_order
 527 |         print ("parameters in the correct order")
 528 |         print (ns_samps_arr[0,:])
 529 |         print (ns_samps_arr[100,:])
 530 |         self.weights = weights
 531 |         print ("no of samps with non zero weight: ",np.count_nonzero(weights))
 532 |         self.nested_samples = ns_samps_arr
 533 |         print ("saving the marginal posterior samples")
 534 |         np.savetxt(path+"post_samples_marginal.csv",self.nested_samples,delimiter=",")
 535 |         np.savetxt(path+"post_samples_marginal_weights.csv",self.weights,delimiter=",")
 536 |     
 537 |         
 538 |     def calc_posterior_mean_and_std_for_comp_maps(self,chunks_to_use):
 539 |         #| Use the samples of the spectras to generate samples of component maps
 540 |         no_of_samples=self.nested_samples.shape[0]
 541 |         print ("no of samples is:",no_of_samples)
 542 |         print ("generating componet map samples")
 543 |         #sample_maps=np.empty((no_of_samples,bayes_eval.no_of_pixels*no_of_comps))
 544 | 
 545 |         
 546 |         #no_of_samples = 5000
 547 |         
 548 |         
 549 |         
 550 | 
 551 |         for i in chunks_to_use:
 552 |             chunk = chunks[i]
 553 | 
 554 |             samp_chunk = self.nested_samples[chunk[0]:chunk[1],:]
 555 |             print ("=======================================")
 556 |             save_name = path+"post_samples_chunk_"+str(int(i+1))+"_comp_maps.npy.gz"
 557 |             try:
 558 |                 os.remove(save_name)
 559 |             except:
 560 |                 print("No pre exisiting comp maps")
 561 | 
 562 |             print (samp_chunk)
 563 |             print ("running for chunk: ",i+1," covering indexes",chunk)
 564 |             
 565 |             
 566 |             sample_maps_ret_block = self.process_chunk(samp_chunk)
 567 | 
 568 |             print (sample_maps_ret_block)
 569 | 
 570 |             print ("saving these maps:","post_samples_chunk_"+str(int(i+1))+"_comp_maps.npy.gz")
 571 |             f = gzip.GzipFile(save_name,"w")
 572 |             np.save(file=f,arr=sample_maps_ret_block)
 573 |             f.close()
 574 |             
 575 |         
 576 |                 
 577 | 
 578 |     def process_chunk(self,sample_chunk):
 579 |         sample_maps_block=[]
 580 |         for i in range(sample_chunk.shape[0]):
 581 |             samp = sample_chunk[i]
 582 |             #print (samp)
 583 |             sample_map = bayes_eval.gen_comp_map_sample(samp,freqs_to_fit_noise,freqs_to_calibrate)
 584 |             #with open(path+"post_samples_comp_maps.csv","a") as f_object:
 585 |             #    writer_object = writer(f_object)
 586 |             #    writer_object.writerow(list(sample_map.flatten()))
 587 |             sample_maps_block.append(sample_map.flatten())
 588 |         return np.array(sample_maps_block)
 589 | 
 590 | #plot the comp maps and their errors
 591 | #######################################################################################################
 592 | samples_generator = generate_component_map_and_sky_map_samples()
 593 | 
 594 | samples_generator.calc_posterior_mean_and_std_for_comp_maps(chunks_to_use_for_run)
 595 | 

```

`Bayesian-Global-Sky-Model-B-GSM-Paper-1/mock_data_v4_mk24_mixed_models_nested_sampling_perfect_noise.py`:

```py
   1 | from random import sample
   2 | import matplotlib as mpl
   3 | mpl.use('Agg')
   4 | from itertools import combinations_with_replacement
   5 | from itertools import product
   6 | import numpy as np
   7 | import healpy as hp
   8 | import matplotlib.pyplot as plt
   9 | from mpl_toolkits.axes_grid1 import make_axes_locatable
  10 | from matplotlib.colors import SymLogNorm
  11 | from matplotlib.colors import LogNorm
  12 | import os
  13 | import pandas
  14 | from numpy import pi, log, sqrt
  15 | import pix_by_pix_mk24_final2 as likelihood
  16 | 
  17 | import scipy.optimize as so
  18 | import math
  19 | import scipy
  20 | 
  21 | import matplotlib.cm as cm
  22 | from scipy.optimize import minimize
  23 | 
  24 | try:
  25 |     import pypolychord
  26 |     from pypolychord.settings import PolyChordSettings
  27 |     from pypolychord.priors import UniformPrior
  28 |     from pypolychord.priors import GaussianPrior
  29 | except:
  30 |     pass
  31 | try:
  32 |     from anesthetic import NestedSamples
  33 | except ImportError:
  34 |     pass
  35 | try:
  36 |     from anesthetic.weighted_pandas import WeightedDataFrame
  37 | except ImportError:
  38 |     pass
  39 | from scipy.optimize import minimize
  40 | import gen_EDGES_beams as gen_beams_and_T_vs_LST_v2
  41 | 
  42 | #declare a random seed
  43 | np.random.seed(0)
  44 | 
  45 | use_perturbed_dataset = True #do we want the input dataset to have calibration errors
  46 | 
  47 | #===================================================================#
  48 | #| SET PARAMS FOR THE SIMULATED DATA AND NESTED SAMPLING
  49 | Max_Nside=32 #the Nside at which to generate the set of maps
  50 | Max_m = hp.nside2npix(Max_Nside)
  51 | no_of_comps = 2
  52 | #which components will have curved powerlaws
  53 | fit_curved=[True,True]#[False,False,False]
  54 | #fit_curved=[False,False,False]
  55 | 
  56 | 
  57 | no_to_fit_curve = np.sum(fit_curved)
  58 | rezero_prior_std=2000
  59 | 
  60 | #params for the prior on the spectra
  61 | spec_min, spec_max = -3.5,1 #the range for the prior on the spectral indexes
  62 | curvature_mean, curvature_std = 0,3 #the range for the prior on the spectral index curvature parameter
  63 | 
  64 | fixed_f0 = 150 #if you dont fit a seperate reference freq for each comp then we fix f0 to this value
  65 | 
  66 | 
  67 | #params for the prior on the true maps
  68 | map_prior_variance_spec_index = -2.6#-2.8#-4
  69 | map_prior_variance_f0 = 408#fixed_f0
  70 | map_prior_std=300#295#1000
  71 | 
  72 | #misalanious params
  73 | calibrate = True
  74 | calibrate_all_but_45_150 = False#True #calibrate all maps in the dataset but the 45 and 150 MHz maps
  75 | calibrate_all = True#False #calibrate every map in the dataset
  76 | 
  77 | use_equal_spaced_LSTs = True
  78 | fit_haslam_noise = False
  79 | subtract_CMB = 0#-2.726
  80 | print_vals_as_calc = False
  81 | 
  82 | 
  83 | 
  84 | reject_criterion = 1e-3#None #how close can two spectral idexes be in value before being rejected
  85 | cond_no_threshold =1e+9
  86 | 
  87 | 
  88 | 
  89 | 
  90 | unobs_marker = -32768
  91 | 
  92 | 
  93 | 
  94 | nlive = 500#2500*no_of_comps
  95 | 
  96 | precision_criterion = 1e-3
  97 | 
  98 | 
  99 | f0=150 #ref freq used for some fitting of spec indexes for plots (not used during any model fitting)
 100 | 
 101 | 
 102 | freqs = np.array([45.0,50.0,60.0,70.0,74.0,80.0,150.0,159.0,408.0]) #the frequencies in MHz of maps used to generate the model
 103 | #specify what to fit for each map
 104 | if calibrate == True:
 105 |     if calibrate_all == True:
 106 |         freqs_to_calibrate = np.array([True,True,True,True,True,True,True,True,True]) #calibrate all the maps
 107 |         freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 108 |     if calibrate_all_but_45_150 == True:
 109 |         freqs_to_calibrate = np.array([False,True,True,True,True,True,False,True,True]) #calibrate all the maps except 45 and 150
 110 |         freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 111 | 
 112 | else:
 113 |     freqs_to_calibrate = np.array([False,False,False,False,False,False,False,False,False])#np.array([True,True,True,True,True,True,False,False])
 114 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 115 | 
 116 | #freqs_for_T_v_LST_comp = np.linspace(40,190,76,dtype="int")
 117 | #freqs_for_T_v_LST_comp = np.array([40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,122.5,125,127.5,130,132.5,135,137.5,140,142.5,145,147.5,150,152.5,155,157.5,160,162.5,165,167.5,170,172.5,175,177.5,180,182.5,185,187.5,190,192.5,195,197.5,200])
 118 | freqs_for_T_v_LST_comp = np.array([40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.0,115.0,120.0,125.0,130.0,135.0,140.0,145.0,150.0,155.0,160.0,165.0,170.0,175.0,180.0,185.0,190.0,195.0,200.0])
 119 | n_spec_pars = 3*no_of_comps #the number of parameters for the spectra (for each comp we have: break_freq, spec_index1, spec_index2)
 120 | 
 121 | 
 122 | nv=len(freqs) #the number of freqs that have maps
 123 | no_of_fitted_noise = np.sum(freqs_to_fit_noise) #the number of freqs at which we fit noise level
 124 | no_of_calibrated = np.sum(freqs_to_calibrate) #the number of freqs at which we fit the calibration
 125 | 
 126 | 
 127 | 
 128 | #set the prior for the noise (on the Haslam map)
 129 | noise_prior_lower, noise_prior_upper = np.array([0.01]),np.array([50])
 130 | 
 131 | print ("noise prior is from:",noise_prior_lower,"to",noise_prior_upper,"Kelvin")
 132 | 
 133 | 
 134 | 
 135 | 
 136 | if np.sum(freqs_to_calibrate)!=0:
 137 |     #set the prior params for the zero levels
 138 |     zero_lev_prior_std = rezero_prior_std*np.ones(np.sum(freqs_to_calibrate))#200*((np.array(freqs)[freqs_to_calibrate]/100)**-2.5)
 139 |     zero_lev_prior_means = np.zeros(np.sum(freqs_to_calibrate))
 140 | 
 141 | 
 142 |     #set the prior params for the scale corrections
 143 |     scale_prior_lower=0.85
 144 |     scale_prior_upper=1.25
 145 |     
 146 | 
 147 |     print ("zero level prior is gauss with mean 0K, stds (Kelvin):")
 148 |     print (zero_lev_prior_std)
 149 |     print ("temp scale prior is uniform from:",scale_prior_lower,"to",scale_prior_upper)
 150 | 
 151 | 
 152 | 
 153 | 
 154 | 
 155 | #true_spec_indexes = [-2.13,-1.86,-1.46]
 156 | #CREATE A DIR TO STORE RESULTS
 157 | #====================================================================#
 158 | #set the file root
 159 | 
 160 | main_label = "_petur:"+str(use_perturbed_dataset)+"_"+str(no_of_comps)+"_comp_cal:"+str(calibrate)+"_rezro_pri_std:"+str(rezero_prior_std)+"_CMB="+str(subtract_CMB)+"_map_pri_std:"+str(map_prior_std)+"_mu:0_map_pri_std_spec_ind="+str(map_prior_variance_spec_index)+"_map_pri_f0="+str(map_prior_variance_f0)+"_cond_no_thres="+str(np.round(np.log10(cond_no_threshold),1))+"_crv_N_std="+str(curvature_std)+"_spec="+str(spec_min)+"_to:"+str(spec_max)#+"_rej_crit="+str(reject_criterion)#+"_nlive="+str(nlive)+"_nrept="+str(nrepeat)+"_precision_criterion="+str(precision_criterion)
 161 | 
 162 | if use_equal_spaced_LSTs==True:
 163 |     #LSTs_for_comparison = np.array([2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])#np.array([0,2,4,6,8,10,12,14,16,18,20,22])#np.array([2.5,18]) #the LSTs in hours at which we will make comparison between the mean sky and EDGES for likelihood calls
 164 |     LSTs_for_comparison = np.linspace(0,24,73)[:-1]
 165 |     print (LSTs_for_comparison)
 166 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 167 |     if calibrate_all==True:
 168 |         root="uni_EDGES_v4_data_mk24_no_of_curved:"+str(no_to_fit_curve)+"_cal_all_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 169 |     else:
 170 |         if calibrate_all_but_45_150==True:
 171 |             root="uni_EDGES_v4_dat_mk24_no_curve:"+str(no_to_fit_curve)+"_no_cal_45_150_f0="+str(fixed_f0)+main_label
 172 | else:
 173 |     #LSTs_for_comparison = np.array([0,1,2,3,4,5,6,7,8,9,10,11,11.25,11.5,11.75,12,12.25,12.5,12.75,13,13.25,13.5,13.75,14,14.25,14.5,14.75,15,15.25,15.5,15.75,16,16.25,16.5,16.75,17,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18,18.25,18.5,18.75,19,19.25,19.5,19.75,20,20.25,20.5,20.75,21,21.25,21.5,21.75,22,22.25,22.5,22.75,23,23.25,23.5,23.75])
 174 |     LSTs_for_comparison = np.array([0,2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])
 175 |     print (LSTs_for_comparison)
 176 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 177 |     
 178 |     root="mk24_extra_uneq_LSTs:"+str(len(LSTs_for_comparison))+"_fixed_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 179 |     
 180 | 
 181 | #make a dir to store the results
 182 | p=os.getcwd()+"/"
 183 | path = p+root+"/"
 184 | try:
 185 |     os.mkdir(root)
 186 | except:
 187 |     pass
 188 | #make a dir to store the results as we run
 189 | root2 = path+"running_results/"
 190 | try:
 191 |     os.mkdir(root2)
 192 | except:
 193 |     pass
 194 | 
 195 | #write a text file containing the key results
 196 | with open(path+"run_details.txt","w") as f:
 197 |     f.write("\n")
 198 |     f.write("a sumary file for the results of the nested sampling run.\n")
 199 |     f.write("=========================================================\n")
 200 |     f.write("f0 = "+str(fixed_f0)+"\n")
 201 |     f.write("map prior std is: "+str(map_prior_std)+" Kelvin\n")
 202 |     f.write("map prior spec index is: "+str(map_prior_variance_spec_index)+"\n")
 203 |     f.write("freqs with maps in input dataset \n")
 204 |     f.write(str(freqs)+"\n")
 205 |     f.write("calibration freqs\n")
 206 |     f.write(str(freqs_to_calibrate)+"\n")
 207 |     f.write("fit noise at freqs\n")
 208 |     f.write(str(freqs_to_fit_noise)+"\n")
 209 |     f.write("freqs for comparison of T vs LST: \n")
 210 |     f.write(str(freqs_for_T_v_LST_comp)+"\n")
 211 |     f.write("LSTs for comparison: \n")
 212 |     f.write(str(LSTs_for_comparison)+"\n")
 213 |     f.close()
 214 | 
 215 | #LOAD THE DATASET AND THE ERROR MAPS
 216 | #====================================================================#
 217 | 
 218 | obs_maps = []
 219 | inv_err_maps = []
 220 | data_err_maps = []
 221 | load_path = p+"mock_data_file/mock_dataset_v4/"
 222 | for i in range(len(freqs)):
 223 |     f=freqs[i]
 224 | 
 225 |     if use_perturbed_dataset==True:
 226 |         fname1 = "noisy_perturbed_sky_"+str(f)
 227 |         err_fname = "perturbed_err_map_"+str(f)
 228 |     else:
 229 |         fname1 = "noisy_sky_"+str(f)
 230 |         err_fname = "err_map_"+str(f)
 231 |     
 232 |     #fname2 = "noise_"+str(f)
 233 | 
 234 |     if freqs_to_fit_noise[i]==False:
 235 |         try:
 236 |             #m1, err_m = np.loadtxt(load_path+fname1), np.loadtxt(load_path+fname2)
 237 |             m1 = np.loadtxt(load_path+fname1)
 238 | 
 239 |             err_m = np.loadtxt(load_path+err_fname)
 240 |         except:
 241 |             print ("cant find the files for freq:",f)
 242 | 
 243 |         
 244 |         #mask out any pixels with negative temps
 245 |         bool_arr = m1<=0
 246 |         err_m[bool_arr] = unobs_marker
 247 |         m1[bool_arr] = unobs_marker
 248 | 
 249 |         
 250 | 
 251 |         inv_err_m = 1/err_m
 252 |         inv_err_m[(err_m==unobs_marker)] = 0
 253 | 
 254 |         m1[m1!=unobs_marker] = m1[m1!=unobs_marker]+subtract_CMB
 255 |         obs_maps.append(m1)
 256 |     
 257 |         inv_err_maps.append(inv_err_m)
 258 | 
 259 |         data_err_maps.append(err_m)
 260 |     else:
 261 |         m1 = np.loadtxt(load_path+fname1)
 262 |         m1[m1!=unobs_marker] = m1[m1!=unobs_marker]+subtract_CMB
 263 |         obs_maps.append(m1)
 264 |     
 265 | 
 266 | obs_maps=np.array(obs_maps)
 267 | inv_err_maps=np.array(inv_err_maps)
 268 | print ("dataset loaded")
 269 | 
 270 | #CREATE THE INVERSE NOISE MATRICES
 271 | #====================================================================#
 272 | #generate the inverse noise covariance matrix for each pixel
 273 | inv_noise_mats = np.empty(shape=(Max_m,len(freqs),len(freqs)))
 274 | for p in range(0,Max_m):
 275 |     inv_stds_for_pixel = np.zeros(len(freqs))
 276 |     inv_stds_for_pixel[~freqs_to_fit_noise] = inv_err_maps[:,p]
 277 |     #print (inv_stds_for_pixel)
 278 | 
 279 |     Np_inv = np.diag(inv_stds_for_pixel**2)
 280 |     #print (Np_inv)
 281 |     inv_noise_mats[p,:,:] = Np_inv
 282 | print ("==============================================================")
 283 | print ("the mean inv_noise mat has diagonal elements of:")
 284 | print (np.diag(np.nanmean(inv_noise_mats,axis=0))/np.diag(np.nanmean(inv_noise_mats,axis=0))[-1])
 285 | print ("max and min for the inv noise mats: ",np.max(inv_noise_mats),np.min(inv_noise_mats[(inv_noise_mats!=0)]))
 286 | 
 287 | print ("fiting a power law to these matrix diagonal elements, gives:")
 288 | logs = np.log(np.diag(np.nanmean(inv_noise_mats,axis=0)))
 289 | 
 290 | log_freqs = np.log(np.array(freqs)/map_prior_variance_f0)
 291 | fun = lambda x: np.nansum((logs - x[0]*log_freqs -x[1])**2)
 292 | res = minimize(fun,[-2.5,1])
 293 | print ("least squares fit to ",res.x)
 294 | print ("diag elements (using this spec) are",np.exp(res.x[1])*(np.array(freqs)/map_prior_variance_f0)**res.x[0]/(np.exp(res.x[1])*(408/map_prior_variance_f0)**res.x[0]))
 295 | print ("diag elements (using our prior assumption) are",(np.array(freqs)/map_prior_variance_f0)**(-2*map_prior_variance_spec_index))
 296 | 
 297 | print ("fitting a curved power spectrum to the noise covars")
 298 | print ("inverse noise matrices created")
 299 | #PLOT THE DATASET
 300 | #====================================================================#
 301 | fig = plt.figure(figsize=(12,16))
 302 | for i in range(len(freqs)):
 303 |     map_i = np.copy(obs_maps[i])
 304 |     ax = plt.subplot(5,3,int(i+1))
 305 |     map_i[(map_i==unobs_marker)]=float("NaN")
 306 |     plt.axes(ax)
 307 |     hp.mollview(map_i,title="mock data freq="+str(freqs[i]),hold=True,notext=True,norm="log")
 308 | 
 309 | 
 310 | 
 311 | plt.savefig(path+"/mock_sky_maps_for_dataset")
 312 | #plt.show()
 313 | fig = plt.figure(figsize=(12,16))
 314 | for i in range(len(freqs)):
 315 |     map_i = np.copy(inv_err_maps[i])
 316 |     ax = plt.subplot(5,3,int(i+1))
 317 |     map_i[(map_i==0)]=float("NaN")
 318 |     plt.axes(ax)
 319 |     hp.mollview(1/map_i,title="input errs freq="+str(freqs[i]),hold=True,notext=True,norm="log")
 320 | 
 321 | 
 322 | 
 323 | plt.savefig(path+"/err_maps_for_dataset")
 324 | freqs=np.array(freqs)
 325 | #plot the priors and the data
 326 | #====================================================================
 327 | log_mean_temps = []
 328 | mean_temps = []
 329 | for i in range(len(freqs)):
 330 |     the_map = obs_maps[i]
 331 |     mean = np.mean(the_map[(the_map!=unobs_marker)])
 332 |     log_mean_temps.append(np.log(mean))
 333 |     mean_temps.append(mean)
 334 | log_mean_temps = np.array(log_mean_temps)
 335 | mean_temps = np.array(mean_temps)
 336 | 
 337 | log_freqs = np.log(freqs/f0)
 338 | fun = lambda x: np.nansum((log_mean_temps - x[0]*log_freqs -x[1])**2)
 339 | res = minimize(fun,[-2.15,np.log(16)])
 340 | print (res)
 341 | 
 342 | fig = plt.figure(figsize=(6,6))#figsize=(12,16))
 343 | #the fitted powerlaw
 344 | targ = np.exp(res.x[1])*((freqs/f0)**res.x[0])
 345 | 
 346 | ax1 = plt.subplot(1,1,1)
 347 | ax1.plot(freqs,targ,c="red",label="fitted power law")
 348 | ax1.scatter(freqs,mean_temps,label="data")
 349 | ax1.set_title("map mean temps")
 350 | ax1.legend(loc="upper right")
 351 | ax1.set_xscale("log")
 352 | ax1.set_yscale("log")
 353 | plt.savefig(path+"/input_map_mean_temps.png")
 354 | 
 355 | 
 356 | #generate a set of pre rotated EDGES beams at each of the frequencies that we want to compare the model to EDGES for
 357 | #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 358 | print ("generating EDGES beams")
 359 | EDGES_beams = np.empty(shape=(obs_maps.shape[1],len(freqs_for_T_v_LST_comp),len(LSTs_for_comparison)))
 360 | j=0
 361 | for f in freqs_for_T_v_LST_comp:
 362 |     print ("gen EDGES beams for freq:",f)
 363 |     beams_at_LSTS = gen_beams_and_T_vs_LST_v2.gen_EDGES_beams_at_LSTs(f,LSTs_for_comparison,Max_Nside)
 364 |     EDGES_beams[:,j,:] = beams_at_LSTS
 365 |     #for i in range(len(LSTs_for_comparison)):
 366 |     #    beam = beams_at_LSTS[:,i]
 367 |         #print (beam)
 368 |         #print (beam.shape)
 369 |         #hp.mollview(beam.flatten(),title="EDGES beam (galactic coords) for freq: "+str(f)+" for LST: "+str(LSTs_for_comparison[i]))
 370 |         #plt.savefig(path+"EDGES_beam_galactic_coords_for_freq:"+str(f)+"_for_LST:"+str(LSTs_for_comparison[i])+".png")
 371 |         #plt.close("all")
 372 |     j+=1
 373 | 
 374 | #Load the mock EDGES T vs LST plots (these are produced by convolving the mock sky map (before any pertubation) with the model beam at that freq)
 375 | #===============================================================================================
 376 | EDGES_temps_at_calib_LSTs_and_freqs = []
 377 | EDGES_errs = []
 378 | print ("generating the EDGES T vs LST traces")
 379 | #we use uniform errors for EDGES
 380 | 
 381 | for f in freqs_for_T_v_LST_comp:
 382 |    # if f<100:
 383 |    #     T_vs_LST, T_vs_LST_errs = gen_beams_and_T_vs_LST_v2.gen_EDGES_low_T_LST_trace(f,LSTs_for_comparison)
 384 |    #     EDGES_temps_at_calib_LSTs_and_freqs.append(T_vs_LST)
 385 |    #     EDGES_errs.append(T_vs_LST_errs)
 386 |    #     fig=plt.figure()
 387 |    #     plt.errorbar(LSTs_for_comparison,T_vs_LST,T_vs_LST_errs)
 388 |    #     plt.savefig(path+"EDGES_T_vs_LST_freq="+str(f)+".png")
 389 |    # if f>=100:
 390 |    #     T_vs_LST, T_vs_LST_errs = gen_beams_and_T_vs_LST_v2.gen_EDGES_high_T_LST_trace(f,LSTs_for_comparison)
 391 |    #     EDGES_temps_at_calib_LSTs_and_freqs.append(T_vs_LST)
 392 |    #     EDGES_errs.append(T_vs_LST_errs)
 393 |    #     fig=plt.figure()
 394 |    #     plt.errorbar(LSTs_for_comparison,T_vs_LST,T_vs_LST_errs)
 395 |    #     plt.savefig(path+"EDGES_T_vs_LST_freq="+str(f)+".png")
 396 |     T_vs_LST = np.loadtxt(load_path+"true_noisy_TvsLST_"+str(f))
 397 |     T_vs_LST_errs = np.loadtxt(load_path+"TvsLST_errs_"+str(f))
 398 |     print ("TvsLST for freq:",f)
 399 |     print (T_vs_LST)
 400 |     print (T_vs_LST_errs)
 401 |     EDGES_temps_at_calib_LSTs_and_freqs.append(T_vs_LST)
 402 |     EDGES_errs.append(T_vs_LST_errs)
 403 |     fig=plt.figure()
 404 |     plt.errorbar(LSTs_for_comparison,T_vs_LST,T_vs_LST_errs)
 405 |     plt.savefig(path+"EDGES_T_vs_LST_freq="+str(f)+".png")
 406 | plt.close("all")
 407 | 
 408 | EDGES_temps_at_calib_LSTs_and_freqs = np.array(EDGES_temps_at_calib_LSTs_and_freqs)
 409 | EDGES_errs = np.array(EDGES_errs)
 410 | 
 411 | 
 412 | 
 413 | 
 414 | #generate the EDGES noise covar mats for each freq assuming noise for each LST is independent of the other LSTs
 415 | EDGES_inv_noise_mats = []
 416 | for i in range(len(freqs_for_T_v_LST_comp)):
 417 |     inv_cov = np.diag(1/EDGES_errs[i,:]**2)
 418 |     EDGES_inv_noise_mats.append(inv_cov)
 419 | EDGES_inv_noise_mats = np.array(EDGES_inv_noise_mats)
 420 | #=====================================================================================================
 421 | #set up the likelihood function
 422 | bayes_eval = likelihood.bayes_mod(obs_maps=obs_maps,obs_freqs=freqs,inv_noise_mats=inv_noise_mats,EDGES_beams=EDGES_beams,EDGES_temps_at_calib_LSTs_and_freqs=EDGES_temps_at_calib_LSTs_and_freqs,EDGES_errs=EDGES_errs,EDGES_inv_noise_mats=EDGES_inv_noise_mats,freqs_for_T_v_LST_comp=freqs_for_T_v_LST_comp,LSTs_for_comparison=LSTs_for_comparison,no_of_comps=no_of_comps,save_root=root2,un_obs_marker=unobs_marker,map_prior_std=map_prior_std,map_prior_spec_index=map_prior_variance_spec_index,map_prior_f0=map_prior_variance_f0)
 423 | 
 424 | 
 425 | 
 426 | no_of_params_for_spec_mod = no_of_comps + no_to_fit_curve 
 427 | print ("we arn't fitting f0: f0=",fixed_f0," no of comps with curved spectra is",no_to_fit_curve," no of params for spectral model is",no_of_params_for_spec_mod)
 428 | 
 429 | 
 430 | #select the spec indexes 
 431 | if no_to_fit_curve==no_of_comps:
 432 |     print ("all comps are curved spec")
 433 |     spec_indexes_select = np.tile(np.array([False,True,False]),no_of_comps)
 434 |     spec_curvature_select =  np.tile(np.array([False,False,True]),no_of_comps)
 435 |     spec_f0_select =  np.tile(np.array([True,False,False]),no_of_comps)
 436 |     print ("spec_curvature_select =",spec_curvature_select)
 437 | else:
 438 |     print ("not all comps are curved spec")
 439 |     spec_indexes_select = np.tile(np.array([False,True,False]),no_of_comps)
 440 |     spec_f0_select =  np.tile(np.array([True,False,False]),no_of_comps)
 441 |     spec_curvature_select = np.array([np.array([False,False,True])*fit_curve_for_comp for fit_curve_for_comp in fit_curved]).flatten() #the indexes for spectral curvature in the final param array
 442 |     print ("spec_curvature_select =",spec_curvature_select)
 443 |     spec_not_curve_select = np.array([np.array([False,False,not_fit_curve_for_comp]) for not_fit_curve_for_comp in ~np.array(fit_curved)]).flatten() #the indexes with no curvature in the final param array
 444 |     print ("spec_not_curve select =",spec_not_curve_select)
 445 | 
 446 | def likelihood1(x):
 447 | 
 448 |     #print ("==================================================================")
 449 |     #print ("=======================LIKELIHOOD CALL============================")
 450 |     #print ("the pars from the prior is len:",len(x))
 451 |     
 452 |     everything_else = x[no_of_params_for_spec_mod:]
 453 |     spec_params = np.empty(3*no_of_comps)
 454 |     spec_params[spec_f0_select] = fixed_f0
 455 |     spec_params[spec_indexes_select] = x[:no_of_comps]
 456 |     spec_params[spec_curvature_select] = x[no_of_comps:no_of_params_for_spec_mod]
 457 |     if no_to_fit_curve!=no_of_comps:
 458 |         spec_params[spec_not_curve_select] = 0
 459 |     
 460 |     #print ("spec_pars:",spec_params)
 461 |     if no_of_fitted_noise!=0:
 462 |         noise_estimates = everything_else[:no_of_fitted_noise]
 463 |     else:
 464 |         noise_estimates = None
 465 |     #print ("noise_estimates:",noise_estimates)
 466 |     if no_of_calibrated!=0:
 467 |         zero_level_estimates = everything_else[no_of_fitted_noise:no_of_fitted_noise+no_of_calibrated]
 468 |         #print ("zero_level_estimates:",zero_level_estimates)
 469 |         scale_estimates = everything_else[no_of_fitted_noise+no_of_calibrated:]
 470 |         #print ("scale_estimates:",scale_estimates)
 471 |     else:
 472 |         zero_level_estimates = None#x[no_of_comps+nv:no_of_comps+2*nv]#[no_of_comps:no_of_comps+len(freqs)]
 473 |         scale_estimates = None#x[no_of_comps+2*nv:no_of_comps+3*nv]
 474 | 
 475 |     #print ("spec params:",spec_params)
 476 |     #print ("noise estimates:",noise_estimates)
 477 |     #print ("scale estimates:",scale_estimates)
 478 |     #print ("zero level estimates:",zero_level_estimates)
 479 |     
 480 |     log_l = bayes_eval.likelihood(spec_params=spec_params,noise_estimates=noise_estimates,freqs_to_fit_noise=freqs_to_fit_noise,scale_estimates=scale_estimates,zero_level_estimates=zero_level_estimates,freqs_to_calibrate=freqs_to_calibrate,joint_prior_func=None,reject_criterion=reject_criterion,print_vals=print_vals_as_calc,condition_no_threshold=cond_no_threshold)
 481 |     #print (log_l)
 482 |     return log_l[-1],[]#[term_to_plot] #return the log_posterior distribtution value for this set of parameters
 483 | 
 484 | #likelihood1(np.array([-2.5,-2.1,-1.8,20]))
 485 | 
 486 | #-------------NESTED SAMPLING PARAMS-------------
 487 | 
 488 | nDims =  int(no_of_params_for_spec_mod + no_of_fitted_noise + 2*no_of_calibrated)
 489 | nrepeat = 5*nDims #the nrepeat is set to 5 times the total number of pars that we fit
 490 | print ("no of dimensions for sampling region is:",nDims)
 491 | nDerived = 0 #we don't derive any parameters 
 492 | settings = PolyChordSettings(nDims, nDerived)
 493 | settings.file_root = root
 494 | settings.nlive = nlive
 495 | settings.nrepeats = nrepeat
 496 | settings.do_clustering = True
 497 | settings.read_resume = True
 498 | settings.write_resume = True
 499 | settings.maximise = False #find the maximum of the poseterior
 500 | settings.precision_criterion = precision_criterion
 501 | 
 502 | prior_lower = spec_min #I have reduced the size of the prior and have centered it on the values that have previously given the best results
 503 | prior_upper = spec_max
 504 | #------------------------------------------------
 505 | 
 506 | #Define a box uniform prior over the specified range of values
 507 |  
 508 | 
 509 | print ("the indexes with the spec index are:",spec_indexes_select)
 510 | def prior(hypercube):
 511 | 
 512 |     #print ("prior called")
 513 |     #print (hypercube)
 514 |     
 515 |     #=====================================================#
 516 |     #define a uniform prior from spec_min to spec_max for the spectral indexes
 517 |     
 518 |     #generate the spectral params
 519 |     ret_array=np.empty(no_of_params_for_spec_mod)
 520 |     spec_par_inits=hypercube[:no_of_params_for_spec_mod]
 521 | 
 522 |     #select only the spectral index terms
 523 |     spec_inits = spec_par_inits[:no_of_comps]
 524 |     #select only the spectral curvature terms
 525 |     curve_inits = spec_par_inits[no_of_comps:]
 526 |         
 527 |     #define the prior for the spectral indexes
 528 |     ret_array[:no_of_comps] = UniformPrior(spec_min, spec_max)(spec_inits)
 529 |     ret_array[no_of_comps:] = GaussianPrior(curvature_mean, curvature_std)(curve_inits)
 530 |         
 531 |     
 532 |     #define a uniform prior from prior_lower to prior_upper
 533 |     uniform_prior_func = lambda vars: UniformPrior(vars[0], vars[1])(vars[2])
 534 | 
 535 |     #define a gaussian prior around a mean with std
 536 |     gauss_prior_func = lambda vars: GaussianPrior(vars[0],vars[1])(vars[2])
 537 |     
 538 |     if no_of_fitted_noise!=0:
 539 |         noise_inits = hypercube[len(spec_par_inits):len(spec_par_inits)+no_of_fitted_noise]
 540 |         for i in range(len(noise_inits)):
 541 |             init_vars = [noise_prior_lower[i],noise_prior_upper[i],noise_inits[i]]
 542 |             #print (init_vars)
 543 |             ret_array= np.append(ret_array, uniform_prior_func(init_vars))
 544 |     else:
 545 |         pass
 546 | 
 547 |     if no_of_calibrated!=0:
 548 |         #print (no_of_calibrated)
 549 |         zero_inits = hypercube[len(spec_par_inits)+no_of_fitted_noise:len(spec_par_inits)+no_of_fitted_noise+no_of_calibrated]
 550 |         scale_inits = hypercube[len(spec_par_inits)+no_of_fitted_noise+no_of_calibrated:]
 551 |         #print (zero_inits)
 552 |         #print (scale_inits)
 553 |         for i in range(len(zero_inits)):
 554 |             init_vars = [zero_lev_prior_means[i],zero_lev_prior_std[i],zero_inits[i]]
 555 |             ret_array= np.append(ret_array, gauss_prior_func(init_vars))
 556 |         #for i in range(len(scale_inits)):
 557 |         #    init_vars = [scale_prior_mean[i],scale_prior_std[i],scale_inits[i]]
 558 |         #    ret_array= np.append(ret_array, gauss_prior_func(init_vars))
 559 |         scales = UniformPrior(scale_prior_lower,scale_prior_upper)(scale_inits)
 560 |         ret_array = np.append(ret_array,scales)
 561 | 
 562 |     else:
 563 |         pass
 564 |     #print ("array returned")
 565 |     #print (ret_array)
 566 |     return ret_array
 567 | 
 568 | if print_vals_as_calc == True:
 569 |     if no_of_comps==2:
 570 |             print ("testing the prior")
 571 |             if no_to_fit_curve==no_of_comps:
 572 |                 init_vals_for_prior = np.append(np.array([0.49,0.5,0.4,0.55]),0.9*np.ones(nDims-no_of_params_for_spec_mod))
 573 |                 
 574 |             else:
 575 |                 init_vals_for_prior = np.concatenate((np.array([0.49,0.5]),0.4*np.ones(no_to_fit_curve),0.9*np.ones(nDims-no_of_params_for_spec_mod)))
 576 |             
 577 |             prior(init_vals_for_prior)
 578 |             print ("testing likelihood")
 579 |             print (likelihood1(prior(init_vals_for_prior)))
 580 |     if no_of_comps==3:
 581 |             print ("testing the prior")
 582 |             if no_to_fit_curve==no_of_comps:
 583 |                 init_vals_for_prior = np.append(np.array([0.49,0.5,0.51,0.4,0.55,0.6]),0.9*np.ones(nDims-no_of_params_for_spec_mod))
 584 |                 
 585 |             else:
 586 |                 init_vals_for_prior = np.concatenate((np.array([0.49,0.5,0.51]),0.4*np.ones(no_to_fit_curve),0.9*np.ones(nDims-no_of_params_for_spec_mod)))
 587 |             
 588 |             prior(init_vals_for_prior)
 589 |             print ("testing likelihood")
 590 |             print (likelihood1(prior(init_vals_for_prior)))
 591 | 
 592 |     
 593 | def dumper(live, dead, logweights, logZ, logZerr):
 594 |     print("Last dead point:", dead[-1])
 595 | 
 596 | 
 597 | print ("running polychord nlive=",settings.nlive)
 598 | #| Run PolyChord
 599 | 
 600 | #| PERFORM NESTED SAMPLING TO DETERMINE LOCATIONS OF CLUSTERS IN THE LOG LIKELIHOOD SURFACE
 601 | output = pypolychord.run_polychord(likelihood1, nDims, nDerived, settings, prior, dumper)
 602 | 
 603 | #| Create a paramnames file
 604 | 
 605 | paramnames = [('p%i' % i, r'n_%i' % i) for i in range(1,nDims+1)]
 606 | output.make_paramnames_files(paramnames)
 607 | 
 608 | #Plot the Nested Sampling results
 609 | 
 610 | mpl.rc("axes", titlesize=20,labelsize=16)
 611 | print ("loading samples")
 612 | try:
 613 |     print (settings.base_dir + '/' + settings.file_root)
 614 |     samples = NestedSamples(root= settings.base_dir + '/' + settings.file_root)
 615 | 
 616 |     print ("samples loaded")
 617 |     #plot the samples
 618 |     print (['p%i' %i for i in range(nDims)])
 619 |     print ("ploting samples")
 620 |     fig, axes = samples.plot_2d(['p%i' %i for i in range(1,nDims+1)])
 621 |     
 622 |     #fig.xticks(fontsize=15)
 623 |     #fig.yticks(fontsize=15)
 624 |     fig.savefig(path+'sampled_posterior.png')
 625 |     #plt.show()
 626 | except:
 627 |     print ("======================ERROR=======================================")
 628 |     print ("could not load files, there is an error check the chains files for this run")

```

`Bayesian-Global-Sky-Model-B-GSM-Paper-1/pix_by_pix_mk24_final2.py`:

```py
   1 | import numpy as np
   2 | import matplotlib.pyplot as plt
   3 | import healpy as hp
   4 | import scipy
   5 | from itertools import combinations
   6 | from astropy.modeling.powerlaws import LogParabola1D
   7 | 
   8 | import traceback
   9 | 
  10 | class bayes_mod():
  11 |     def __init__(self,obs_maps,obs_freqs,inv_noise_mats,EDGES_beams,EDGES_temps_at_calib_LSTs_and_freqs,EDGES_errs,EDGES_inv_noise_mats,freqs_for_T_v_LST_comp,LSTs_for_comparison,save_root,un_obs_marker=-32768,no_of_comps=2,map_prior_std=1000,map_prior_spec_index=-2.6,map_prior_f0=200):
  12 |         
  13 |         """variables are:
  14 |         
  15 |         obs_maps: a np array containing the observed maps in healpix form (entry [v,:] should be the map at freq v). NOTE: maps must contain the same no of pixels
  16 |         obs_freqs: a 1d np array containing the frequencies of the observed maps
  17 |         inv_noise_mats: the inverse noise covariance matrices for each pixel shape is (npix,nv,nv)
  18 |         un_obs_marker: the value used to indicate that a pixel is unobserved
  19 |         EDGES_beams: a set of beams at all freqs and LSTs that we want to use for calibrating the model
  20 |         EDGES_temps_at_calib_LSTs_and_freqs: the temperature of the sky as seen by EDGES for all the LSTs and freqs we use for calibration
  21 |         EDGES_errs: the uncertainty on the temps observed by EDGES
  22 |         freqs_for_T_v_LST_comp: the freqs at whice we will compare our model predictions to the EDGES T v LST traces
  23 |         no_of_comps: an integer number of components to use
  24 |         map_prior_std: the amplitude range (in Kelvin) over which we integrated the component maps during marginalisation default is 1000.
  25 |         
  26 |         NOTE: for this version of the code component spectra are single-index power-laws. Future versions will need to use spectra objects, as this will allow us to have arbitary parametarisation"""
  27 |         
  28 |         self.freqs = obs_freqs#/f0
  29 | 
  30 |         self.data_vecs = np.empty((obs_maps.shape[1],obs_maps.shape[0],1))
  31 |         
  32 |         self.Ninvs=inv_noise_mats
  33 | 
  34 |         
  35 |         #misalainous params
  36 |         self.no_of_comps = no_of_comps
  37 |         self.no_of_pixels = obs_maps.shape[1]
  38 |         self.no_of_freqs = len(obs_freqs)
  39 | 
  40 |         self.D_T_inv = np.diag((1/((obs_freqs/map_prior_f0)**(2*map_prior_spec_index)))) # the assumed spectral behaviour of the prior inverse covariance for the true sky maps 
  41 |         
  42 |         self.map_prior_var = map_prior_std**2 #the assumed prior variance for the true sky map at the reference frequency
  43 | 
  44 |         print ("===============================================")
  45 |         print ("true map prior variance at ref freq of: ",map_prior_f0," is: ",self.map_prior_var," i.e. std of: ",map_prior_std)
  46 |         print ("the assumed spectral behaviour of the prior inverse covariance for the true sky maps (as a matrix) is:")
  47 |         print (np.diag(self.D_T_inv))
  48 |         print ("spec index for this prior inv covar is")
  49 |         print (np.log(np.diag(self.D_T_inv)[:-1])/(np.log(self.freqs[:-1])-np.log(408)))
  50 |         
  51 |         
  52 | 
  53 |         unobs_pixs = obs_maps==un_obs_marker #a boolian array that indicates which pixels have been observed 
  54 |         self.no_of_obs_map = np.count_nonzero(~unobs_pixs,axis=0)
  55 |         self.pix_to_disregard = self.no_of_obs_map<self.no_of_comps #for any pixel that has fewer observations than comps to fit, mask out predictions
  56 |         self.pix_to_regard = self.no_of_obs_map>=1 #only include pixels with at least 1 observed frequency in the likelihood calculations
  57 | 
  58 |         self.fully_obs_pix = self.no_of_obs_map==self.no_of_freqs
  59 |         self.partial_obs_pix = self.pix_to_regard * (self.no_of_obs_map<self.no_of_freqs)
  60 | 
  61 | 
  62 |         self.no_of_defined_pix = np.sum(self.pix_to_regard)
  63 |         self.pix_for_map_gen = self.no_of_obs_map>=self.no_of_comps
  64 |         self.no_of_pix_for_map_gen = np.sum(self.pix_for_map_gen)
  65 |         self.no_of_partial_obs_pix = np.sum(self.partial_obs_pix)
  66 |         self.no_of_fully_obs_pix = np.sum(self.fully_obs_pix)
  67 |         print ("===============================================")
  68 |         print ("no of pix with >=1 observation")
  69 |         print (self.no_of_defined_pix)
  70 |         print ("no of pix with < n_v but >=1 observations")
  71 |         print (self.no_of_partial_obs_pix)
  72 |         print ("no of fully observed pixels")
  73 |         print (self.no_of_fully_obs_pix)
  74 |     
  75 |         
  76 | 
  77 |         self.mean_Ninv = np.mean(inv_noise_mats[self.pix_to_regard],axis=0)
  78 |         #generate the inverse noise matrix for each pixel
  79 |         #self.Ninvs = np.empty((self.no_of_pixels,Ninv_temp.shape[0],Ninv_temp.shape[1]))
  80 |         for p in range(0,self.no_of_pixels):
  81 | 
  82 |             #get the data vectors into the right order for numpy array broadcasting rules
  83 |             self.data_vecs[p,:,0] = obs_maps[:,p]
  84 |             self.data_vecs_trans = np.reshape(self.data_vecs,(self.no_of_pixels,1,obs_maps.shape[0]))
  85 |         
  86 |         
  87 | 
  88 |         
  89 |        
  90 |         try:
  91 |             self.freqs_for_T_v_LST_comp = freqs_for_T_v_LST_comp
  92 |             self.EDGES_beams = (1/(4*np.pi))*EDGES_beams*hp.nside2pixarea(nside=hp.npix2nside(self.no_of_pixels))
  93 |             self.EDGES_noise_det_term = -np.sum(np.log(2*np.pi*EDGES_errs))
  94 |         except:
  95 |             print ("freqs cant be divided setting to None")
  96 |             self.freqs_for_T_v_LST_comp = None
  97 |             self.EDGES_beams = None
  98 |             self.EDGES_noise_det_term = None
  99 |             
 100 |         self.EDGES_temps = EDGES_temps_at_calib_LSTs_and_freqs
 101 |         self.EDGES_inv_noise_mats = EDGES_inv_noise_mats
 102 |         self.EDGES_errs = EDGES_errs
 103 |         self.LSTs_for_comparison = LSTs_for_comparison
 104 |         print ("EDGES freqs to check are:")
 105 |         print (self.freqs_for_T_v_LST_comp)
 106 |         print ("LSTs to test:")
 107 |         print (self.LSTs_for_comparison)
 108 |         #print ("EDGES Low Band powerlaw temps at LSTs")
 109 |         #print (self.EDGES_temps)
 110 |         #print ("EDGES Low Band errors")
 111 |         #print (EDGES_errs)
 112 |         print ("EDGES det term:",self.EDGES_noise_det_term)
 113 |         self.save_root = save_root
 114 |         self.select_index1 = np.tile(np.array([False,True,False]),no_of_comps)
 115 |         self.select_index2 = np.tile(np.array([False,False,True]),no_of_comps)
 116 | 
 117 |     def gen_A(self,spec_params,freqs_for_mix_mat,print_vals=False):
 118 | 
 119 |         
 120 | 
 121 |         A = np.empty(shape=(len(freqs_for_mix_mat),self.no_of_comps))
 122 |         for c in range(self.no_of_comps):
 123 |             pars = spec_params[3*c:3*(c+1)]
 124 |             #print (pars)
 125 |             spec = LogParabola1D(1,pars[0],-pars[1],-pars[2])(freqs_for_mix_mat)
 126 |             A[:,c] = spec
 127 | 
 128 |         #print (A)
 129 |         if print_vals==True:
 130 |             print ("=======================")
 131 |             print ("printing the A matrix and plotting comp spectra for pars:")
 132 |             print (spec_params)
 133 |             print ("A:")
 134 |             print (A)
 135 |             fig = plt.figure(figsize=(6,12))
 136 |             for i in range(self.no_of_comps):
 137 |                 plt.subplot(self.no_of_comps,1,i+1)
 138 |                 plt.plot(freqs_for_mix_mat,A[:,i],label="spec_for_comp="+str(i+1))
 139 |                 plt.legend(loc="upper right")
 140 |             plt.suptitle("spectral indexes are: "+str(spec_params))
 141 |             plt.savefig(self.save_root+"spec="+str(spec_params)+"_A_specs.png")
 142 |             plt.close("all")
 143 |         return A
 144 |     def likelihood(self,spec_params,noise_estimates,freqs_to_fit_noise,scale_estimates,zero_level_estimates,freqs_to_calibrate,joint_prior_func,reject_criterion=None,print_vals=False,condition_no_threshold=1e+7,height_threshold=0.1):
 145 |         """params are:
 146 |         spec_params: the spectral indexes to use for components
 147 |         reject_criterion: if two or more spec indexes are too close in value return -inf
 148 |         """
 149 |         
 150 | 
 151 |         fail_ret = -1e+30#float("NaN")#-np.inf # the value to return if the spec_params fail one of our tests
 152 | 
 153 |         #print_vals=True
 154 | 
 155 |         #check if any of the spectral indexes are repeated or are within the rejection crition of another spectral index  
 156 |         if reject_criterion==None:
 157 |             pass
 158 |         else:
 159 |             #check for spec indexes (we dont check for the spectral curvature or the ref freqs these can be identicle between spectra)
 160 |             sps = spec_params[self.select_index1]
 161 | 
 162 |             within_reject = np.any([np.count_nonzero(sps[:i] >= (sps[i]-reject_criterion)) for i in range(1,len(sps))])
 163 |             if within_reject==True:
 164 |                 #print (spec_params,"  rejected for 1")
 165 |                 #return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]#[-np.inf,-np.inf,-np.inf,-np.inf]
 166 |                 return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]
 167 |             
 168 |             
 169 |         
 170 |         
 171 |         #*****************************************SET UP THE VARIOUS MATRICES**************************************************
 172 | 
 173 |         #generate the mixing matrix
 174 |         A = self.gen_A(spec_params,self.freqs)#self.freqs[:,np.newaxis] ** spec_params[np.newaxis,:] #this models the component spectra as single-index power-laws
 175 | 
 176 |         
 177 |        
 178 | 
 179 | 
 180 |         #generate the data calibration matrix and vector
 181 |         if np.sum(freqs_to_calibrate)==0:
 182 |             #print ("nothing to calibrate")
 183 | 
 184 |             a=np.identity(self.no_of_freqs)
 185 |             a_inv = np.identity(self.no_of_freqs)
 186 |             b=np.zeros((self.no_of_freqs,1))
 187 |         else:
 188 |             a_estimates = np.ones(self.no_of_freqs)
 189 |             a_estimates[freqs_to_calibrate] = scale_estimates
 190 |             #print (a_estimates)
 191 |             b_estimates = np.zeros(self.no_of_freqs)
 192 |             b_estimates[freqs_to_calibrate] = zero_level_estimates
 193 |             a = np.diag(a_estimates)
 194 |             a_inv = np.diag(1/a_estimates)
 195 |             b = np.reshape(b_estimates,(self.no_of_freqs,1))
 196 | 
 197 |         #print (a)
 198 |         #print (a_inv)
 199 |         #print (b)
 200 |         #generate the template noise covariance matrix
 201 |         if np.sum(freqs_to_fit_noise)==0:
 202 |             inv_noise_mats = self.Ninvs[self.pix_to_regard]
 203 |             mean_noise_mat = self.mean_Ninv
 204 |         else:
 205 |             #fill in the diagonals corresponding to the freqs we are fitting the noise
 206 |             temp_diag = np.zeros(self.no_of_freqs)
 207 |             temp_diag[freqs_to_fit_noise] = 1/noise_estimates**2
 208 |             temp = np.diag(temp_diag)
 209 | 
 210 |             inv_noise_mats = self.Ninvs[self.pix_to_regard] + temp
 211 |             mean_noise_mat = self.mean_Ninv + temp
 212 |             #print (mean_noise_mat)
 213 |         
 214 |         #rescale the noise to account for the calibration
 215 |         if np.sum(freqs_to_calibrate)==0:
 216 |             calibrated_inv_noise_mats = inv_noise_mats
 217 |             test_mat = A.T @ mean_noise_mat @ A 
 218 |         else:
 219 |             calibrated_inv_noise_mats = a_inv.T @ inv_noise_mats @ a_inv
 220 |             test_mat = A.T @ a_inv @ mean_noise_mat @ a_inv @ A
 221 | 
 222 | 
 223 | 
 224 |         #compute the Lambda_p matrices for all pixels with at least 1 observation
 225 |         Lambda_ps = A.T @ calibrated_inv_noise_mats @ A
 226 | 
 227 |         #compute S^T D^-1 S
 228 |         try:
 229 |             STDS = A.T @ self.D_T_inv @ A
 230 | 
 231 |             STDS_inv = np.linalg.inv(STDS)
 232 |         except:
 233 |             #print ("couldnt compute STS_inv")
 234 |             return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]
 235 |         #compute the matrix for the determinant 
 236 |         mat_for_det = self.map_prior_var*(STDS_inv @ Lambda_ps) + np.identity(self.no_of_comps)
 237 | 
 238 |         #compute the matrix for the inverse 
 239 |         mat_for_inv = Lambda_ps + (1/self.map_prior_var)*STDS
 240 | 
 241 |         #test that the mean condition number doesn't exceed our threshold
 242 |         mean_cond_no = np.mean(np.linalg.cond(mat_for_inv))
 243 |         overall_test_val = mean_cond_no<condition_no_threshold
 244 |         if overall_test_val == False:
 245 |             #print (mean_cond_no)
 246 |             #print (overall_test_val)
 247 |             #print ("mean condition number exceeded threshold")
 248 |             #return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]#[-np.inf,-np.inf,-np.inf,-np.inf] #if the matrix is singular return -inf
 249 |             return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]
 250 | 
 251 |         
 252 |         #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 253 |         #compute the determinant term
 254 |         try:
 255 |             #calculate the term 2 (-ln(det(Lambda_p)))
 256 |             t2s = -1 * np.linalg.slogdet(mat_for_det)[1] #only calculate for the pixels that have >=self.no_of_comps observations
 257 |         except np.linalg.LinAlgError:
 258 |             #print ("failed to calculate log determinant")
 259 |             #return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]#[-np.inf,-np.inf,-np.inf,-np.inf] #if the matrix is singular return -inf
 260 |             return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]
 261 | 
 262 |         t2 = np.sum(t2s)
 263 |         #t2 = 0
 264 | 
 265 |         #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 266 |         #generate the calibrated data vectors
 267 |         if np.sum(freqs_to_calibrate)==0:
 268 |             cal_dat_vecs = self.data_vecs[self.pix_to_regard]
 269 |             cal_dat_vecs_trans = np.reshape(cal_dat_vecs,(self.no_of_defined_pix,1,self.no_of_freqs))
 270 |         else:
 271 |             cal_dat_vecs = a @ self.data_vecs[self.pix_to_regard] + b
 272 |             cal_dat_vecs_trans = np.reshape(cal_dat_vecs,(self.no_of_defined_pix,1,self.no_of_freqs))
 273 |         
 274 |         #generate the bp vectors
 275 |         bps = A.T @ calibrated_inv_noise_mats @ cal_dat_vecs
 276 |         bps_trans = np.reshape(bps,(self.no_of_defined_pix,1,self.no_of_comps))
 277 | 
 278 |         #add on the offsets caussed by the non zero mean of the map prior
 279 |         shifted_bps = bps #+ self.C0_inv_mu
 280 |         shifted_bps_trans = bps_trans #+ self.C0_inv_mu_trans
 281 |         #print ("============================================================")
 282 |         #print ("============================================================")
 283 |         #print ("============================================================")
 284 |         #print (bps)
 285 |         #print ("============================================================")
 286 |         #print (shifted_bps)
 287 |         #print ("============================================================")
 288 |         #print (bps_trans)
 289 |         #print ("============================================================")
 290 |         #print (shifted_bps_trans)
 291 | 
 292 |         #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 293 |         try:
 294 |             #calculate term 1 (bp.T @ Cp_inv @ bp) for each pixel
 295 |             mean_comp_maps_for_obs_pix = np.linalg.solve(mat_for_inv, shifted_bps)
 296 |             t1s = shifted_bps_trans @ mean_comp_maps_for_obs_pix 
 297 | 
 298 |         except np.linalg.LinAlgError:
 299 |             #print ("failed to calc mean comp maps")
 300 |             return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]
 301 | 
 302 |         t1 = np.sum(t1s)
 303 | 
 304 |         #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 305 |         try:
 306 |             t3as = cal_dat_vecs_trans @ calibrated_inv_noise_mats @ cal_dat_vecs
 307 | 
 308 |             t3bs = 2*np.pi*inv_noise_mats[inv_noise_mats!=0]
 309 | 
 310 |             t3 = np.sum(np.log(t3bs)) - np.sum(t3as)
 311 |             #print ("new value for t3 is:",t3)
 312 |             #print ("old version gave:",np.sum(t3bs) - np.sum(t3as))
 313 |         except:
 314 |             return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]
 315 | 
 316 |         #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 317 |         #compute the mean of the conditional distribution of sky maps for this set of spectra and calibration pars
 318 |         try:
 319 |             #print ("step 1a")
 320 |             mean_comp_maps = np.zeros(shape=(self.no_of_pixels,self.no_of_comps,1))
 321 |             mean_comp_maps[self.pix_to_regard] = mean_comp_maps_for_obs_pix
 322 |             #print ("mean comp maps")
 323 |             #print (mean_comp_maps)
 324 |             #compute the mean of the posterior of the skys for this set of parameters at all observed freqs
 325 |             #print ("step 2a")
 326 |             spectral_mix_mat_for_test = self.gen_A(spec_params,self.freqs_for_T_v_LST_comp,print_vals=print_vals)#self.freqs_for_T_v_LST_comp[:,np.newaxis] ** spec_params[np.newaxis,:]
 327 |             #print ("spectral mix mat")
 328 |             #print (spectral_mix_mat_for_test)
 329 |             mean_sky_preds = spectral_mix_mat_for_test @ mean_comp_maps
 330 |             mean_sky_preds[~self.pix_to_regard] = float("NaN")
 331 |             #print ("mean sky preds")
 332 |             #print (mean_sky_preds)
 333 |             #print (mean_sky_preds.shape)
 334 |             #print (self.no_of_pixels)
 335 | 
 336 |             #multiply the skys by the EDGES beam at these freqs
 337 |             #print ("step 3a")
 338 |             convolved_sky_preds = mean_sky_preds * self.EDGES_beams
 339 | 
 340 |             #compute the integrated sky temp for each freq and LST in the freqs to compare
 341 |             integrated_skys = np.nansum(convolved_sky_preds,axis=0)
 342 | 
 343 |             #diff = integrated_skys - self.EDGES_temps #the difference between the model and the EDGES obs at all freqs and LSTs
 344 |             #diff_trans = np.reshape(diff,(diff.shape[0],1,diff.shape[1]))
 345 |             #diff = np.reshape(diff,(diff.shape[0],diff.shape[1],1))
 346 |             #print (diff[0,:,:])
 347 |             #print (diff_trans[0,:,:])
 348 |             #print ("step 4a")
 349 |             EDGES_likelihood_t2s = ((integrated_skys-self.EDGES_temps)/self.EDGES_errs)**2
 350 |             EDGES_likelihood_t2 = -1*np.sum(EDGES_likelihood_t2s)#-1 * np.sum(diff_trans @ self.EDGES_inv_noise_mats @ diff)
 351 |             EDGES_int_sky_temps_log_likelihood = EDGES_likelihood_t2 + self.EDGES_noise_det_term#(self.no_of_defined_pix/len(self.LSTs_for_comparison))*(EDGES_likelihood_t2 + self.EDGES_noise_det_term)
 352 |             
 353 |             #print (sky_pred_at_45)
 354 |             val=np.random.uniform(0,1)
 355 |             #if val<=0.0001:
 356 |             #    print_vals=True
 357 | 
 358 |             if print_vals==True:
 359 |                 #print ("step 1")
 360 |                 #plot the component maps for these params
 361 |                 for i in range(self.no_of_comps):
 362 |                     cm = mean_comp_maps[:,i]
 363 | 
 364 |                     plt.subplot(1,self.no_of_comps,i+1)
 365 |                     hp.mollview(cm.flatten(),title="mean comp "+str(i+1),hold=True)
 366 |                 plt.savefig(self.save_root+"spec="+str(spec_params)+"_zeros="+str(zero_level_estimates.round(3))+"_scales="+str(scale_estimates.round(3))+"_comps.png")
 367 |                 
 368 |                 plt.close("all")
 369 |                 #plt.show()
 370 | 
 371 |                 
 372 | 
 373 |             
 374 |             
 375 |                 fig = plt.figure(figsize=(16,12))
 376 |                 for i in range(len(self.freqs_for_T_v_LST_comp)):
 377 |                     true = self.EDGES_temps[i,:]
 378 |                     errs = self.EDGES_errs[i,:]
 379 |                     pred = integrated_skys[i,:]
 380 |                     
 381 |                     plt.subplot(15,10,i+1)
 382 |                     hp.mollview(mean_sky_preds[:,i].flatten(),title="sky pred at "+str(np.round(self.freqs_for_T_v_LST_comp[i],1)),hold=True)
 383 | 
 384 |                 plt.suptitle("spectral indexes are: "+str(spec_params)+"\nCalibration pars, Zeros: "+str(zero_level_estimates)+" scales: "+str(scale_estimates))
 385 |                 fig.subplots_adjust(hspace=0.5)
 386 |                 plt.savefig(self.save_root+"spec="+str(spec_params)+"_zeros="+str(zero_level_estimates.round(3))+"_scales="+str(scale_estimates.round(3))+"_skys.png")
 387 |                 #plt.show()
 388 |                 plt.close("all")
 389 | 
 390 |                 #plot the T vs LST for each of the test freqs
 391 |                 fig = plt.figure(figsize=(16,12))
 392 |                 for i in range(len(self.freqs_for_T_v_LST_comp)):
 393 |                     true = self.EDGES_temps[i,:]
 394 |                     errs = self.EDGES_errs[i,:]
 395 |                     pred = integrated_skys[i,:]
 396 |                     
 397 |                     plt.subplot(15,10,i+1)
 398 |                     plt.errorbar(self.LSTs_for_comparison,true,errs,label="EDGES obs")
 399 |                     plt.plot(self.LSTs_for_comparison,pred,label="model pred")
 400 |                     plt.title("freq="+str(np.round(self.freqs_for_T_v_LST_comp[i],1))+"MHz")
 401 |                     plt.legend(loc="upper left")
 402 |                     #print ("==============================")
 403 |                     #print ("true")
 404 |                     #print (true)
 405 |                     #print ("pred")
 406 |                     #print (pred)
 407 |                 plt.suptitle("spectral indexes are: "+str(spec_params)+"\nCalibration pars, Zeros: "+str(zero_level_estimates)+" scales: "+str(scale_estimates))
 408 |                 fig.subplots_adjust(hspace=0.5)
 409 |                 print ("saving the T vs LST plots")
 410 |                 plt.savefig(self.save_root+"spec="+str(spec_params)+"_zeros="+str(zero_level_estimates.round(3))+"_scales="+str(scale_estimates.round(3))+"_T_vs_LST.png")
 411 |                 print ("saved")
 412 |                 #plt.show()
 413 |                 plt.close("all")
 414 |                     #plt.show()
 415 |             #compute the likelihood
 416 |             #print (integrated_skys.shape)
 417 |             #print (self.EDGES_temps.shape)
 418 |             #print (self.EDGES_inv_noise_mats.shape)
 419 |             #diff = integrated_skys - self.EDGES_temps #the difference between the model and the EDGES obs at all freqs and LSTs
 420 |             #diff_trans = np.reshape(diff,(diff.shape[0],1,diff.shape[1]))
 421 |             #print (diff[0,:])
 422 |             #print (diff.bps_trans[0,:,:])
 423 | 
 424 | 
 425 |         except:
 426 |             traceback.print_exc()
 427 |             #print ("failed to comp the simulated EDGES obs")
 428 |             return [fail_ret,fail_ret,fail_ret,fail_ret,fail_ret]
 429 |         
 430 |         prior_term = 0
 431 |         tot_log_likelihood = t1 + t2 + t3 + EDGES_int_sky_temps_log_likelihood
 432 |         posterior = tot_log_likelihood + prior_term
 433 | 
 434 |         #print_vals=False
 435 |         if print_vals==True:
 436 |             print ("===================================================")
 437 |             print ("spec params:",spec_params)
 438 |             
 439 |             print ("term 2 (det) is",t2)
 440 |             print ("det term + prior:",prior_term + t2)
 441 |             print ("term 1 (bp.T@(Lambda_p+C0_inv)^-1 @bp is",t1)
 442 |             print ("term 3 (gp) is",t3)
 443 |             
 444 |             print ("log likelihood of observing the EDGES vals")
 445 |             print ("term from EDGES Low band temps =",EDGES_int_sky_temps_log_likelihood)
 446 |             print ("likelihood:",tot_log_likelihood)
 447 |             print ("posterior: ",posterior)
 448 | 
 449 |         
 450 | 
 451 |         return [t1, t2 + prior_term, tot_log_likelihood, posterior]
 452 |         
 453 | 
 454 | 
 455 |     def gen_comp_map_sample(self,params,freqs_to_fit_noise,freqs_to_calibrate):
 456 |         """calculate the mean and variance of the conditional distribution of component maps P(Mp|dp,S) for a specified set of spectral params.
 457 |         
 458 |         we then return a sample component map set drawn from this conditional distribution"""
 459 | 
 460 |         spec_params = params[:3*self.no_of_comps]
 461 |         #generate the mixing matrix
 462 |         A = self.gen_A(spec_params,self.freqs)#self.freqs[:,np.newaxis] ** spec_params[np.newaxis,:] #this models the component spectra as single-index power-laws
 463 | 
 464 |         #generate the data calibration matrix and vector
 465 |         if np.sum(freqs_to_calibrate)==0:
 466 |             #print ("nothing to calibrate")
 467 | 
 468 |             a=np.identity(self.no_of_freqs)
 469 |             a_inv = np.identity(self.no_of_freqs)
 470 |             b=np.zeros((self.no_of_freqs,1))
 471 |         else:
 472 |             zero_level_estimates = params[3*self.no_of_comps+np.sum(freqs_to_fit_noise):3*self.no_of_comps+np.sum(freqs_to_fit_noise)+np.sum(freqs_to_calibrate)]
 473 |             scale_estimates = params[3*self.no_of_comps+np.sum(freqs_to_fit_noise)+np.sum(freqs_to_calibrate):]
 474 |             
 475 |             a_estimates = np.ones(self.no_of_freqs)
 476 |             a_estimates[freqs_to_calibrate] = scale_estimates
 477 |             #print (a_estimates)
 478 |             b_estimates = np.zeros(self.no_of_freqs)
 479 |             b_estimates[freqs_to_calibrate] = zero_level_estimates
 480 |             a = np.diag(a_estimates)
 481 |             a_inv = np.diag(1/a_estimates)
 482 |             b = np.reshape(b_estimates,(self.no_of_freqs,1))
 483 | 
 484 |         #print (a)
 485 |         #print (a_inv)
 486 |         #print (b)
 487 |         #generate the template noise covariance matrix
 488 |         if np.sum(freqs_to_fit_noise)==0:
 489 |             inv_noise_mats = self.Ninvs[self.pix_to_regard]
 490 |             mean_noise_mat = self.mean_Ninv
 491 |         else:
 492 |             #fill in the diagonals corresponding to the freqs we are fitting the noise
 493 |             temp_diag = np.zeros(self.no_of_freqs)
 494 |             noise_estimates = params[3*self.no_of_comps:3*self.no_of_comps+np.sum(freqs_to_fit_noise)]
 495 |             temp_diag[freqs_to_fit_noise] = 1/noise_estimates**2
 496 |             temp = np.diag(temp_diag)
 497 | 
 498 |             inv_noise_mats = self.Ninvs[self.pix_to_regard] + temp
 499 |             mean_noise_mat = self.mean_Ninv + temp
 500 |             #print (mean_noise_mat)
 501 |         
 502 |         #rescale the noise to account for the calibration
 503 |         if np.sum(freqs_to_calibrate)==0:
 504 |             calibrated_inv_noise_mats = inv_noise_mats
 505 |             test_mat = A.T @ mean_noise_mat @ A 
 506 |         else:
 507 |             calibrated_inv_noise_mats = a_inv.T @ inv_noise_mats @ a_inv
 508 |             test_mat = A.T @ a_inv @ mean_noise_mat @ a_inv @ A
 509 | 
 510 |         #generate the set of Lambda_p matrixes
 511 |         Lambda_ps = A.T @ calibrated_inv_noise_mats @ A
 512 | 
 513 |         #compute S^T S
 514 |         try:
 515 |             STDS = A.T @ self.D_T_inv @ A
 516 | 
 517 |             STDS_inv = np.linalg.inv(STDS)
 518 |         except:
 519 |             print ("couldnt compute STS_inv")
 520 |             return np.zeros(shape=y.shape)
 521 | 
 522 |         #compute the matrix for the inverse 
 523 |         mat_for_inv = Lambda_ps + (1/self.map_prior_var)*STDS
 524 | 
 525 |         #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 526 |         #generate the calibrated data vectors
 527 |         if np.sum(freqs_to_calibrate)==0:
 528 |             cal_dat_vecs = self.data_vecs[self.pix_to_regard]
 529 |             cal_dat_vecs_trans = np.reshape(cal_dat_vecs,(self.no_of_defined_pix,1,self.no_of_freqs))
 530 |         else:
 531 |             cal_dat_vecs = a @ self.data_vecs[self.pix_to_regard] + b
 532 |             cal_dat_vecs_trans = np.reshape(cal_dat_vecs,(self.no_of_defined_pix,1,self.no_of_freqs))
 533 |         
 534 |         #generate the bp vectors
 535 |         bps = A.T @ calibrated_inv_noise_mats @ cal_dat_vecs
 536 |         bps_trans = np.reshape(bps,(self.no_of_defined_pix,1,self.no_of_comps))
 537 | 
 538 | 
 539 |         #add on the offsets caussed by the non zero mean of the map prior
 540 |         shifted_bps = bps #+ self.C0_inv_mu
 541 | 
 542 |         y = np.zeros(shape=(self.no_of_pixels,self.no_of_comps,1))
 543 |     
 544 | 
 545 |         
 546 |         
 547 |         
 548 |         
 549 |         
 550 |         mean_maps = np.empty(shape=y.shape)
 551 |         #calculate the mean of the distribution for each pixel of this component-map-set
 552 |         try:#print (spec_params)
 553 |             
 554 |             mean_maps_for_defined_pix = np.linalg.solve(mat_for_inv, shifted_bps)#Lambda_p_invs @ bps
 555 |             
 556 | 
 557 |             mean_maps[self.pix_to_regard]=mean_maps_for_defined_pix
 558 |             
 559 |             
 560 |         except:
 561 |         
 562 |             print ("failed to gen mean comp maps")
 563 |             return np.zeros(shape=y.shape)
 564 |         
 565 | 
 566 |         #calculate cholesky decomp of Lambda_p_ivns matrices
 567 |         #NOTE: we only do this for pixels with at least k observations (i.e. where Lambda_p should be non singular).
 568 |         #For the other pixels we dont add any noise. This is fine as we mask out the predictions for these anyway.
 569 |         try:
 570 |             #invert all of the Lambda_p matrices
 571 |             C_ps =  np.linalg.inv(mat_for_inv) #the covariance matrix for the conditional distribution of comp maps for each pixel
 572 |             
 573 |             C_ps_chole = np.linalg.cholesky(C_ps) #take the cholesky decomp of the covar for each pixel
 574 |         except:
 575 |             print ("cholesky decomp failed")
 576 |             return mean_maps
 577 |         #generate noise
 578 |         #noise drawn from a unit variance gaussian (covar matraix is unit variance and is no_of_comps by no_of_comps, we generate no_of_pixels samples from this)
 579 |         standard_noise_for_defined_pix = np.random.standard_normal(size=(self.no_of_defined_pix,self.no_of_comps,1))
 580 |         
 581 |         #multiply by cholesky decomposition of the covariance matrix for each pixel
 582 |         noise_for_pix_to_regard = C_ps_chole @ standard_noise_for_defined_pix
 583 |         
 584 |         #asign the noise to the defined pixels (unobs pixels are all set to zero)
 585 |         noise = np.zeros(shape=y.shape)
 586 |         noise[self.pix_to_regard] = noise_for_pix_to_regard
 587 |         #print (noise.shape)
 588 | 
 589 |         samps = mean_maps + noise
 590 |         
 591 |         samps[~self.pix_to_regard]=float("NaN") #mask out regions of the sky with <k obs
 592 |         #print (samps)
 593 |         return samps
 594 | 

```

`Bayesian-Global-Sky-Model-B-GSM-Paper-1/plot_mean_comps_and_skys_v2.py`:

```py
   1 | from random import sample
   2 | import matplotlib as mpl
   3 | #mpl.use('Agg')
   4 | from itertools import combinations_with_replacement
   5 | from itertools import product
   6 | import numpy as np
   7 | import healpy as hp
   8 | import matplotlib.pyplot as plt
   9 | from mpl_toolkits.axes_grid1 import make_axes_locatable
  10 | from matplotlib.colors import SymLogNorm
  11 | from matplotlib.colors import LogNorm
  12 | import os
  13 | import pandas
  14 | from numpy import pi, log, sqrt
  15 | import pix_by_pix_mk24_final2 as likelihood
  16 | 
  17 | import scipy.optimize as so
  18 | import math
  19 | import scipy
  20 | 
  21 | import matplotlib.cm as cm
  22 | from scipy.optimize import minimize
  23 | 
  24 | try:
  25 |     import pypolychord
  26 |     from pypolychord.settings import PolyChordSettings
  27 |     from pypolychord.priors import UniformPrior
  28 |     from pypolychord.priors import GaussianPrior
  29 | except:
  30 |     pass
  31 | try:
  32 |     from anesthetic import NestedSamples
  33 | except ImportError:
  34 |     pass
  35 | try:
  36 |     from anesthetic.weighted_pandas import WeightedDataFrame
  37 | except ImportError:
  38 |     pass
  39 | from scipy.optimize import minimize
  40 | #import gen_beams_and_T_vs_LST as gen_beams_and_T_vs_LST_v2
  41 | from line_profiler import LineProfiler
  42 | import matplotlib
  43 | from astropy.modeling.powerlaws import LogParabola1D
  44 | from csv import writer
  45 | 
  46 | #========================================================================================#
  47 | #compute mean posterior component amplitude maps
  48 | calc_mean_comps = False
  49 | #compute mean posterior sky predictions
  50 | calc_skys = True
  51 | 
  52 | #========================================================================================#
  53 | #a very large number (~60 thousand) posterior samples have very low weighting <10^-50
  54 | #setting this to True means we only use the more highly weighted posterior samples
  55 | use_limited_post = True
  56 | 
  57 | 
  58 | #specify the frequencies we will generate posterior predictions at
  59 | #========================================================================================#
  60 | #test_freqs = np.array([45.0,50,60])
  61 | #test_freqs = np.array([70.0,74,80])
  62 | #test_freqs = np.array([150.0,159,408])
  63 | #test_freqs = np.array([47.5,75,100,140])
  64 | #test_freqs = np.array([200.0,250.0])
  65 | test_freqs = np.array([300.0,350.0])
  66 | #test_freqs = np.array([400.0])
  67 | 
  68 | print ("test freqs are:",test_freqs)
  69 | #declare a random seed
  70 | np.random.seed(0)
  71 | use_perturbed_dataset = True #do we want the input dataset to have calibration errors
  72 | 
  73 | #===================================================================#
  74 | #| SET PARAMS FOR THE SIMULATED DATA AND NESTED SAMPLING
  75 | Max_Nside=32 #the Nside at which to generate the set of maps
  76 | Max_m = hp.nside2npix(Max_Nside)
  77 | no_of_comps = 2
  78 | fit_curved=[True,True]
  79 | 
  80 | 
  81 | no_to_fit_curve = np.sum(fit_curved)
  82 | rezero_prior_std=2000
  83 | 
  84 | #params for the prior on the spectra
  85 | spec_min, spec_max = -3.5,1 #the range for the prior on the spectral indexes
  86 | curvature_mean, curvature_std = 0,3 #the range for the prior on the spectral index curvature parameter
  87 | 
  88 | #params for fitting the reference frequency 
  89 | 
  90 | fixed_f0 = 150 #if you dont fit a seperate reference freq for each comp then we fix f0 to this value
  91 | 
  92 | 
  93 | #params for the prior on the true maps
  94 | map_prior_variance_spec_index = -2.6
  95 | map_prior_variance_f0 = 408#fixed_f0
  96 | map_prior_std=300
  97 | 
  98 | calibrate = True
  99 | calibrate_all_but_45_150 = False#True #calibrate all maps in the dataset but the 45 and 150 MHz maps
 100 | calibrate_all = True#False #calibrate every map in the dataset
 101 | 
 102 | 
 103 | use_equal_spaced_LSTs = True
 104 | fit_haslam_noise = False
 105 | subtract_CMB = 0#-2.726
 106 | print_vals_as_calc = False
 107 | 
 108 | 
 109 | 
 110 | reject_criterion = 1e-3#None #how close can two spectral idexes be in value before being rejected
 111 | cond_no_threshold =1e+9
 112 | 
 113 | 
 114 | 
 115 | 
 116 | test_LSTs = np.linspace(0,24,73)[:-1]#np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])#np.array([0,2,4,6,8,10,12,14,16,18,20,22]) #the LSTs in hours at which we will make comparison between the mean sky and EDGES for likelihood calls
 117 | #test_freqs = [47.5,75]#,250,300,350]#[45,50.005,59.985,70.007,73.931,79.960,100,125,150,159,200,408]
 118 | #test_freqs = [100,140,200]
 119 | #test_freqs = [250,300,350]
 120 | 
 121 | #np.array([45.0,50.0,60.0,70,74,80,150,159,408])
 122 | #test_freqs = [70.0,74.0,80.0]
 123 | #test_freqs = [150.0,159.0,408.0]
 124 | 
 125 | 
 126 | 
 127 | unobs_marker = -32768
 128 | 
 129 | nlive = 500#2500*no_of_comps
 130 | 
 131 | precision_criterion = 1e-3
 132 | 
 133 | 
 134 | #specify what to fit for each map
 135 | if calibrate_all_but_45_150 == True:
 136 |     freqs_to_calibrate = np.array([False,True,True,True,True,True,False,True,True]) #calibrate all the maps except the 45 and 150 MHz
 137 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 138 | 
 139 | if calibrate_all == True:
 140 |     freqs_to_calibrate = np.array([True,True,True,True,True,True,True,True,True]) #calibrate all the maps except the 45 and 150 MHz
 141 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 142 | 
 143 | if calibrate == False:
 144 |     freqs_to_calibrate = np.array([False,False,False,False,False,False,False,False,False])#np.array([True,True,True,True,True,True,False,False])
 145 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 146 | 
 147 | 
 148 | main_label = "_petur:"+str(use_perturbed_dataset)+"_"+str(no_of_comps)+"_comp_cal:"+str(calibrate)+"_rezro_pri_std:"+str(rezero_prior_std)+"_CMB="+str(subtract_CMB)+"_map_pri_std:"+str(map_prior_std)+"_mu:0_map_pri_std_spec_ind="+str(map_prior_variance_spec_index)+"_map_pri_f0="+str(map_prior_variance_f0)+"_cond_no_thres="+str(np.round(np.log10(cond_no_threshold),1))+"_crv_N_std="+str(curvature_std)+"_spec="+str(spec_min)+"_to:"+str(spec_max)#+"_rej_crit="+str(reject_criterion)#+"_nlive="+str(nlive)+"_nrept="+str(nrepeat)+"_precision_criterion="+str(precision_criterion)
 149 | 
 150 | if use_equal_spaced_LSTs==True:
 151 |     #LSTs_for_comparison = np.array([2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])#np.array([0,2,4,6,8,10,12,14,16,18,20,22])#np.array([2.5,18]) #the LSTs in hours at which we will make comparison between the mean sky and EDGES for likelihood calls
 152 |     LSTs_for_comparison = np.linspace(0,24,73)[:-1]
 153 |     print (LSTs_for_comparison)
 154 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 155 |     if calibrate_all==True:
 156 |         root="uni_EDGES_v4_data_mk24_no_of_curved:"+str(no_to_fit_curve)+"_cal_all_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 157 |     else:
 158 |         if calibrate_all_but_45_150==True:
 159 |             root="uni_EDGES_v4_dat_mk24_no_curve:"+str(no_to_fit_curve)+"_no_cal_45_150_f0="+str(fixed_f0)+main_label
 160 | else:
 161 |     #LSTs_for_comparison = np.array([0,1,2,3,4,5,6,7,8,9,10,11,11.25,11.5,11.75,12,12.25,12.5,12.75,13,13.25,13.5,13.75,14,14.25,14.5,14.75,15,15.25,15.5,15.75,16,16.25,16.5,16.75,17,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18,18.25,18.5,18.75,19,19.25,19.5,19.75,20,20.25,20.5,20.75,21,21.25,21.5,21.75,22,22.25,22.5,22.75,23,23.25,23.5,23.75])
 162 |     LSTs_for_comparison = np.array([0,2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])
 163 |     print (LSTs_for_comparison)
 164 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 165 |     
 166 |     root="mk24_extra_uneq_LSTs:"+str(len(LSTs_for_comparison))+"_fixed_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 167 | 
 168 | 
 169 | p=os.getcwd()+"/"
 170 | path = p+root+"/"
 171 | 
 172 | #load the marginal samples and comp map samples
 173 | print ("loading marginal samples and weights")
 174 | marg_samps = np.loadtxt(path+"post_samples_marginal.csv",delimiter=",")
 175 | weights = np.loadtxt(path+"post_samples_marginal_weights.csv",delimiter=",")
 176 | 
 177 | import gzip
 178 | 
 179 | if use_limited_post == True:
 180 |     #we throw away the 60000 lowest weighted samples
 181 |     #this significantly reduces memory requierments 
 182 |     weights = weights[60000:]
 183 |     marg_samps = marg_samps[60000:,:]
 184 | 
 185 |     print ("loading the component amplitude samples")
 186 |     start = 7
 187 | else:
 188 |     #we use all posterior samples 
 189 |     start = 1
 190 | for i in range(start,10):
 191 |     name = path+"post_samples_chunk_"+str(int(i))+"_comp_maps.npy.gz"
 192 | 
 193 |     try:
 194 |         print ("opening chunk:",i)
 195 |         f = gzip.GzipFile(name,"r")
 196 | 
 197 |         chunk = np.load(f)
 198 |         f.close()
 199 | 
 200 |         if i ==start:
 201 |             maps = chunk
 202 |         else:
 203 |             maps = np.append(maps,chunk,axis=0)
 204 |     except:
 205 |         print ("no file for this chunk")
 206 |     print (maps.shape)
 207 | 
 208 | no_of_pixels = int(maps.shape[1]/no_of_comps)
 209 | no_of_samps = maps.shape[0]
 210 | 
 211 | if calc_mean_comps==True:
 212 |     samples_of_maps_DF = WeightedDataFrame(maps,weight=weights)#np.array(samples_of_maps)
 213 |     print (samples_of_maps_DF)
 214 | 
 215 |     posterior_mean = samples_of_maps_DF.mean() #,axis=0)
 216 |     posterior_std = samples_of_maps_DF.std()#np.nanstd(samples_of_maps,axis=0)#/np.sqrt(nested_samples.shape[0])
 217 | 
 218 | 
 219 |     #convert dataframe to np array
 220 |     posterior_mean = posterior_mean.to_numpy()
 221 |     posterior_std = posterior_std.to_numpy()
 222 |     #print (posterior_mean.shape)
 223 | 
 224 |     recovered_map = np.empty((no_of_pixels,no_of_comps))
 225 |     err_map = np.empty((no_of_pixels,no_of_comps))
 226 |     bool_arr_template = np.zeros(no_of_comps,dtype=bool)
 227 |     #comps_arr = np.empty((no_of_comps,no_of_comps,no_of_pixels))
 228 |     for c in range(no_of_comps):
 229 |         bool_arr_temp = np.copy(bool_arr_template)
 230 |         bool_arr_temp[c] = 1
 231 | 
 232 |         bool_arr = np.tile(bool_arr_temp,no_of_pixels)
 233 | 
 234 |         recovered_map[:,c]=posterior_mean[bool_arr]
 235 |         err_map[:,c]=posterior_std[bool_arr]#[c*no_of_pixels:(c+1)*no_of_pixels]
 236 |     #    print (recovered_map)
 237 | 
 238 | #    comps_arr[:,c,:] = maps[:,bool_arr]
 239 | 
 240 |     print (recovered_map)
 241 |     #save the mean componetent maps and their errors
 242 |     for i in range(no_of_comps):
 243 |         rec_comp, rec_comp_errs = recovered_map[:,i], err_map[:,i]
 244 | 
 245 |         np.savetxt(path+'posterior_mean_for_comp_'+str(i+1),rec_comp,delimiter=",")
 246 |         np.savetxt(path+'posterior_std_for_comp_'+str(i+1),rec_comp_errs,delimiter=",")
 247 |     #plot the component maps
 248 |     fig = plt.figure(figsize=(15,10))
 249 | 
 250 |     for i in range(no_of_comps):
 251 |         ax = plt.subplot2grid((2,no_of_comps),(0,i))
 252 | 
 253 |         plt.axes(ax)
 254 |         hp.mollview(recovered_map[:,i],title="recovered mean comp "+str(i+1),hold=True,notext=True)#,min=-10,max=10)#np.max(c_map_1))
 255 | 
 256 |     for i in range(no_of_comps):
 257 |         ax = plt.subplot2grid((2,no_of_comps),(1,i))
 258 | 
 259 |         plt.axes(ax)
 260 |         hp.mollview(err_map[:,i],title="errs mean comp "+str(i+1),hold=True)#,notext=True,min=-10,max=10)#np.max(c_map_1))
 261 | 
 262 |     plt.savefig(path+"recovered_comp_comparison.png")
 263 | 
 264 | #gen the skys
 265 | bool_arr_template = np.zeros(no_of_comps,dtype=bool)
 266 | comps_arr = np.empty((no_of_samps,no_of_comps,no_of_pixels))
 267 | for c in range(no_of_comps):
 268 |     bool_arr_temp = np.copy(bool_arr_template)
 269 |     bool_arr_temp[c] = 1
 270 | 
 271 |     bool_arr = np.tile(bool_arr_temp,no_of_pixels)
 272 |     comps_arr[:,c,:] = maps[:,bool_arr]
 273 | 
 274 | def gen_samp_sky(index):
 275 | 
 276 |     spec_params = marg_samps[index,:3*no_of_comps]
 277 | 
 278 |     comp_maps_samp = comps_arr[index]
 279 | 
 280 |     specs = np.empty((len(test_freqs),no_of_comps))
 281 |     for i in range(no_of_comps):
 282 |         sp = spec_params[3*i:3*(i+1)]
 283 |     #    print (sp)
 284 |         spec = LogParabola1D(1,sp[0],-sp[1],-sp[2])(test_freqs)
 285 |     #    print (spec)
 286 |         specs[:,i]=spec
 287 |     specs = np.array(specs)
 288 |     #print (specs)
 289 |     s1 = specs @ comp_maps_samp
 290 | 
 291 |     #print (s1)
 292 |     #print (s1.flatten())
 293 |     return s1.flatten()
 294 | 
 295 | if calc_skys == True:
 296 |     print ("generating sky samples")
 297 |     sky_samps = []
 298 |     for i in range(no_of_samps):
 299 |         sky_samps.append(gen_samp_sky(i))
 300 | 
 301 |     sky_samps = np.array(sky_samps)
 302 | 
 303 |     samples_of_maps_DF = WeightedDataFrame(sky_samps,weight=weights)#np.array(samples_of_maps)
 304 |     print (samples_of_maps_DF)
 305 | 
 306 |     posterior_mean = samples_of_maps_DF.mean() #,axis=0)
 307 |     posterior_std = samples_of_maps_DF.std()#np.nanstd(samples_of_maps,axis=0)#/np.sqrt(nested_samples.shape[0])
 308 | 
 309 | 
 310 |     #convert dataframe to np array
 311 |     posterior_mean = posterior_mean.to_numpy()
 312 |     posterior_std = posterior_std.to_numpy()
 313 | 
 314 |     for i in range(len(test_freqs)):
 315 |         f=test_freqs[i]
 316 |         sky = posterior_mean[i*no_of_pixels:(i+1)*no_of_pixels]
 317 |         sky_err = posterior_std[i*no_of_pixels:(i+1)*no_of_pixels]
 318 | 
 319 |         np.savetxt(path+"bayesian_pred_"+str(f)+"MHz",sky,delimiter=",")
 320 |         np.savetxt(path+"bayesian_errs_"+str(f)+"MHz",sky_err,delimiter=",")
 321 | 
 322 |         hp.mollview(sky,title=str(f))
 323 |         plt.savefig(path+str(f)+".png")
 324 | 
 325 |         hp.mollview(sky_err,title=str(f))
 326 |         plt.savefig(path+str(f)+"_err.png")
 327 | 

```

`Bayesian-Global-Sky-Model-B-GSM-Paper-1/plot_posterior_results_for_paper_mixed_model_with_crosshair_v2.py`:

```py
   1 | from random import sample
   2 | import matplotlib as mpl
   3 | #mpl.use('Agg')
   4 | from itertools import combinations_with_replacement
   5 | from itertools import product
   6 | import numpy as np
   7 | import healpy as hp
   8 | import matplotlib.pyplot as plt
   9 | from mpl_toolkits.axes_grid1 import make_axes_locatable
  10 | from matplotlib.colors import SymLogNorm
  11 | from matplotlib.colors import LogNorm
  12 | import os
  13 | import pandas
  14 | from numpy import pi, log, sqrt
  15 | #import pix_by_pix_mk19b as likelihood
  16 | 
  17 | import scipy.optimize as so
  18 | import math
  19 | import scipy
  20 | 
  21 | import matplotlib.cm as cm
  22 | from scipy.optimize import minimize
  23 | 
  24 | try:
  25 |     import pypolychord
  26 |     from pypolychord.settings import PolyChordSettings
  27 |     from pypolychord.priors import UniformPrior
  28 |     from pypolychord.priors import GaussianPrior
  29 | except:
  30 |     pass
  31 | try:
  32 |     from anesthetic import NestedSamples
  33 | except ImportError:
  34 |     pass
  35 | try:
  36 |     from anesthetic.weighted_pandas import WeightedDataFrame
  37 | except ImportError:
  38 |     pass
  39 | from scipy.optimize import minimize
  40 | import gen_EDGES_beams as gen_beams_and_T_vs_LST_v2
  41 | from line_profiler import LineProfiler
  42 | import matplotlib
  43 | from fgivenx import plot_contours, plot_lines
  44 | from fgivenx import samples_from_getdist_chains
  45 | from astropy.modeling.powerlaws import LogParabola1D
  46 | unobs_marker = -32768
  47 | p=os.getcwd()+"/"
  48 | Nside = 32
  49 | 
  50 | no_of_comps = 2
  51 | #=======================================================================
  52 | fit_curved=[True,True]
  53 | calibrate_all = True#False
  54 | calibrate_all_but_45_150 = False#True
  55 | 
  56 | rezero_prior_std=2000
  57 | subtract_CMB = 0
  58 | #params for the prior on the spectra
  59 | spec_min, spec_max = -3.5,1 #the range for the prior on the spectral indexes
  60 | curvature_mean, curvature_std = 0,3 #the range for the prior on the spectral index curvature parameter
  61 | 
  62 | #params for fitting the reference frequency 
  63 | 
  64 | fixed_f0 = 150 #if you dont fit a seperate reference freq for each comp then we fix f0 to this value
  65 | 
  66 | 
  67 | #params for the prior on the true maps
  68 | map_prior_variance_spec_index = -2.6
  69 | map_prior_variance_f0 = 408#fixed_f0
  70 | map_prior_std=300
  71 | 
  72 | no_to_fit_curve = np.sum(fit_curved)
  73 | 
  74 | #=======================================================================
  75 | #true spectral parameters
  76 | true_specs=[[-2.6,0],[-2.1,-0.5]]
  77 | #true components
  78 | true_c1, true_c2 = np.loadtxt(p+"mock_data_file/true_comp1"), np.loadtxt(p+"mock_data_file/true_comp2")
  79 | true_c1 += 50-np.nanmin(true_c1)
  80 | true_c2 += 10-np.nanmin(true_c2)
  81 | 
  82 | #true calibration values
  83 | true_alphas = [1,0.95,1.05,1.06,1.1,1.15,1,1.18,0.95]
  84 | true_betas = [0,600,400,300,-200,-300,0,-30,5]
  85 | 
  86 | #reference frequency for the spectra
  87 | f0=150
  88 | #freqs for the map comparison
  89 | fmaps_for_test = np.array([47.5,75,100,140,200,250,300,350])#[45.0,50.0,60.0,70.0,74.0,80.0,150.0,159.0,408.0]#[47.5,75,100,140,200,250,300,350]
  90 | 
  91 | #freqs for the T_vs_LST comparison
  92 | fT_vs_LST = [55,75,95,120,140,160,180,220,300,380]
  93 | LSTs = np.linspace(0,24,49)
  94 | 
  95 | #freqs for the mean_temps comparisons
  96 | fmeans = [55,75,95,120,140,160,180,220,300,380]
  97 | 
  98 | main_label = "_petur:True_"+str(no_of_comps)+"_comp_cal:True_rezro_pri_std:"+str(rezero_prior_std)+"_CMB="+str(subtract_CMB)+"_map_pri_std:"+str(map_prior_std)+"_mu:0_map_pri_std_spec_ind="+str(map_prior_variance_spec_index)+"_map_pri_f0="+str(map_prior_variance_f0)+"_cond_no_thres=9.0_crv_N_std="+str(curvature_std)+"_spec="+str(spec_min)+"_to:"+str(spec_max)#+"_rej_crit="+str(reject_criterion)#+"_nlive="+str(nlive)+"_nrept="+str(nrepeat)+"_precision_criterion="+str(precision_criterion)
  99 | 
 100 | if calibrate_all==True:
 101 |         file_name="uni_EDGES_v4_data_mk24_no_of_curved:"+str(no_to_fit_curve)+"_cal_all_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 102 | else:
 103 |     if calibrate_all_but_45_150==True:
 104 |         file_name="uni_EDGES_v4_dat_mk24_no_curve:"+str(no_to_fit_curve)+"_no_cal_45_150_f0="+str(fixed_f0)+main_label
 105 | #file_name = "v3b_data_mk24_LSTs:72_f0=150_perturbed_dataset:True_2_comp_cal_all_inc_408:True_rezro_pri_std:2000_CMB=0_map_pri_std:500_mu:0_map_pri_spec_ind=-2.6_map_pri_f0=408_cond_no_thres=9.0_crv_N_std=2_spec=-3.5_to:1"
 106 | 
 107 | root = p+file_name
 108 | print ("================================================================================")
 109 | print (file_name)
 110 | print ("================================================================================")
 111 | 
 112 | #plot the true compontnets and the spectra used to generate the synthetic dataset
 113 | #================================================================================
 114 | plt_freqs = np.logspace(1.653,2.7,200)
 115 | 
 116 | fig=plt.figure(figsize=(7,4.5))
 117 | ax1 = plt.subplot2grid((2,2),(0,0))
 118 | plt.axes(ax1)
 119 | hp.mollview(true_c1,title="Component 1 true",hold=True,notext=True,norm="log")
 120 | ax2 = plt.subplot2grid((2,2),(1,0))
 121 | plt.axes(ax2)
 122 | hp.mollview(true_c2,title="Component 2 true",hold=True,notext=True,norm="log")
 123 | 
 124 | ax3 = plt.subplot2grid((2,2),(0,1))
 125 | ax3.set_title("Spectrum Component 1",fontsize=12)
 126 | ax3.plot(plt_freqs,LogParabola1D(1,150,2.6,0)(plt_freqs))
 127 | ax3.set_ylabel("scaling",fontsize=12)
 128 | ax3.set_yticks([0.01,0.1,1,10,100],[0.01,0.1,1,10,100])
 129 | ax3.set_ylim(0.01,100)
 130 | ax3.set_xscale('log')
 131 | ax3.set_yscale('log')
 132 | 
 133 | ax4 = plt.subplot2grid((2,2),(1,1))
 134 | ax4.set_title("Spectrum Component 2",fontsize=12)
 135 | ax4.plot(plt_freqs,LogParabola1D(1,150,2.1,0.5)(plt_freqs))
 136 | ax4.set_ylabel("scaling",fontsize=12)
 137 | ax4.set_xlabel("Freq (MHz)",fontsize=12)
 138 | ax4.set_yticks([0.01,0.1,1,10,100],[0.01,0.1,1,10,100])
 139 | ax4.set_ylim(0.01,100)
 140 | ax4.set_xscale('log')
 141 | ax4.set_yscale('log')
 142 | plt.rcParams.update({'font.size':12})
 143 | plt.subplots_adjust(wspace=0.5)
 144 | plt.subplots_adjust(hspace=0.5)
 145 | #plt.tight_layout()
 146 | #plt.rcParams['font.size'] = 15
 147 | #plt.show()
 148 | plt.savefig(root+"/true_comps_and_specs.pdf")
 149 | 
 150 | #produce a plot showing the production of the synthtic EDGES data
 151 | #=====================================================================================
 152 | true_sky45 = np.loadtxt(p+"/mock_data_file/mock_dataset_v4/true_sky_45.0")
 153 | true_TvsLST45 = np.loadtxt(p+"/mock_data_file/mock_dataset_v4/true_noisy_TvsLST_45.0")
 154 | beam45_at_18h = gen_beams_and_T_vs_LST_v2.gen_EDGES_beams_at_LSTs(45,[18],32).flatten()
 155 | beam45_at_18h[(beam45_at_18h<0)]=0
 156 | true_sky45[(true_sky45<=0)]=float("NaN")
 157 | 
 158 | 
 159 | fig=plt.figure(figsize=(6,4))
 160 | ax1 = plt.subplot2grid((1,1),(0,0))
 161 | plt.axes(ax1)
 162 | hp.mollview(true_sky45,title="Synthetic Sky 45MHz\n(no calibration uncertainty)",hold=True,notext=True,norm="log")
 163 | plt.savefig(root+"/syth_sky_45.pdf")
 164 | 
 165 | fig=plt.figure(figsize=(6,4))
 166 | ax2 = plt.subplot2grid((1,1),(0,0))
 167 | plt.axes(ax2)
 168 | hp.mollview(beam45_at_18h,title="Beam Model 45MHz\n(for LST=18h)",hold=True,notext=True)#,norm="log")
 169 | plt.savefig(root+"/syth_beam_45.pdf")
 170 | 
 171 | fig=plt.figure(figsize=(6.5,4.5))
 172 | ax3 = plt.subplot2grid((1,1),(0,0))
 173 | ax3.set_title("Synthetic EDGES data for 45MHz",fontsize=15)
 174 | ax3.plot(np.linspace(0,24,73)[:-1],true_TvsLST45)
 175 | ax3.set_ylabel("Antenna Temp (K)",fontsize=15)
 176 | ax3.set_xlabel("LST (hours)")
 177 | plt.savefig(root+"/syth_EDGES_45.pdf")
 178 | plt.show()
 179 | #========================================================
 180 | #define the nested samping labels and no of dimensions
 181 | if calibrate_all==True:
 182 |     nDims = no_of_comps + no_to_fit_curve + 18
 183 |     if no_to_fit_curve == 2:
 184 |         param_names = [r"$\beta_1$",r"$\beta_2$",r"$\gamma_1$",r"$\gamma_2$",r"$b_{45}$",r"$b_{50}$",r"$b_{60}$",r"$b_{70}$",r"$b_{74}$",r"$b_{80}$",r"$b_{150}$",r"$b_{159}$",r"$b_{408}$",r"$a_{45}$",r"$a_{50}$",r"$a_{60}$",r"$a_{70}$",r"$a_{74}$",r"$a_{80}$",r"$a_{150}$",r"$a_{159}$",r"$a_{408}$"]
 185 |         true_param_vals_for_all = [-2.6,-2.1,0,-0.5,0,600,400,300,-200,-300,0,-30,5,1,0.95,1.05,1.06,1.1,1.15,1,1.18,0.95]
 186 |     if no_to_fit_curve == 1:
 187 |         param_names = [r"$\beta_1$",r"$\beta_2$",r"$\gamma_2$",r"$b_{45}$",r"$b_{50}$",r"$b_{60}$",r"$b_{70}$",r"$b_{74}$",r"$b_{80}$",r"$b_{150}$",r"$b_{159}$",r"$b_{408}$",r"$a_{45}$",r"$a_{50}$",r"$a_{60}$",r"$a_{70}$",r"$a_{74}$",r"$a_{80}$",r"$a_{150}$",r"$a_{159}$",r"$a_{408}$"]
 188 |         true_param_vals_for_all = [-2.6,-2.1,-0.5,0,600,400,300,-200,-300,0,-30,5,1,0.95,1.05,1.06,1.1,1.15,1,1.18,0.95]
 189 | if calibrate_all_but_45_150==True:
 190 |     nDims = no_of_comps + no_to_fit_curve + 14
 191 | 
 192 |     if no_to_fit_curve == 2:
 193 |         param_names = [r"$\beta_1$",r"$\beta_2$",r"$\gamma_1$",r"$\gamma_2$",r"$b_{50}$",r"$b_{60}$",r"$b_{70}$",r"$b_{74}$",r"$b_{80}$",r"$b_{159}$",r"$b_{408}$",r"$a_{50}$",r"$a_{60}$",r"$a_{70}$",r"$a_{74}$",r"$a_{80}$",r"$a_{159}$",r"$a_{408}$"]
 194 |         true_param_vals_for_all = [-2.6,-2.1,0,-0.5,600,400,300,-200,-300,-30,5,0.95,1.05,1.06,1.1,1.15,1.18,0.95]
 195 |     if no_to_fit_curve == 1:
 196 |         param_names = [r"$\beta_1$",r"$\beta_2$",r"$\gamma_2$",r"$b_{50}$",r"$b_{60}$",r"$b_{70}$",r"$b_{74}$",r"$b_{80}$",r"$b_{159}$",r"$b_{408}$",r"$a_{50}$",r"$a_{60}$",r"$a_{70}$",r"$a_{74}$",r"$a_{80}$",r"$a_{159}$",r"$a_{408}$"]
 197 |         true_param_vals_for_all = [-2.6,-2.1,-0.5,600,400,300,-200,-300,-30,5,0.95,1.05,1.06,1.1,1.15,1.18,0.95]
 198 |         upper_lower_bounds=[[-2.62,-2.598],[-2.11,-2.09],[-0.51,-0.48]]
 199 | #load the posterior samples
 200 | samples = NestedSamples(root= p+"chains" + '/' + file_name)
 201 | column_names_in_dataframe = list(samples.columns.values)
 202 | print ("column names for dataframe:", column_names_in_dataframe)
 203 | nested_samples = []
 204 | for i in range(nDims):
 205 |     print (i,column_names_in_dataframe[i])
 206 |     nested_samples.append(list(samples.loc[:,column_names_in_dataframe[i]]))
 207 |     
 208 | nested_samples = np.array(nested_samples).T#samples.loc[:,par_llamo].to_numpy()
 209 | 
 210 | weights = samples.weight
 211 | samps_with_weights = WeightedDataFrame(nested_samples,weight=weights)
 212 | mean_params = samps_with_weights.mean().to_numpy()
 213 | std_params = samps_with_weights.std().to_numpy()
 214 | no_of_samples = nested_samples.shape[0]
 215 | print ("=========================================================")
 216 | print ("mean param solultions")
 217 | print (mean_params)
 218 | print ("stds")
 219 | print (std_params)
 220 | print ("mean param solutions (with no weighting of samples)")
 221 | print (np.mean(nested_samples,axis=0))
 222 | print ("=========================================================")
 223 | 
 224 | 
 225 | #===========================================================
 226 | #plot a histogram of the Log evidence values for the posterior
 227 | samples.gui()
 228 | plt.show()
 229 | logZ = samples.logZ() #average log evidence for the model
 230 | print (type(logZ),logZ)
 231 | logL = samples.logL.to_numpy()
 232 | 
 233 | plt.hist(logL,bins=np.linspace(-1000000,-800000,100),weights=weights)
 234 | 
 235 | #samples.logL.plot(kind="hist",weights=weights)
 236 | plt.show()
 237 | plt.hist(logL,bins=np.linspace(-883000,-881000,1000),weights=weights)
 238 | plt.xlabel("log(L)")
 239 | plt.ylabel("occurrence")
 240 | plt.show()
 241 | 
 242 | #plt.hist(logL,bins=np.linspace(-809000,-808900,1000))#,weights=weights)
 243 | #plt.show()
 244 | #===========================================================
 245 | #produce a corner plot for the posterior samples
 246 | 
 247 | #define the axis labels for the posterior plot
 248 | tick_labels = []
 249 | upper_lower_bounds = []
 250 | for i in range(len(mean_params)):
 251 |     #define the axis labels to be at the posterior mean value and +/- 5 sigma from this
 252 |     lower_lab, mean_lab ,upper_lab = mean_params[i] - 2.5*std_params[i], mean_params[i], mean_params[i] + 2.5*std_params[i]
 253 | 
 254 |     lower_lab, mean_lab, upper_lab = np.round(lower_lab,3), np.round(mean_lab,3), np.round(upper_lab,3)
 255 |     lower_b,upper_b = mean_params[i] - 5*std_params[i], mean_params[i] + 5*std_params[i]
 256 | 
 257 |     lower_b, upper_b = np.round(lower_b,3), np.round(upper_b,3)
 258 | 
 259 |     if std_params[i]>=1:
 260 |         lower_lab, mean_lab, upper_lab = np.round(lower_lab,0), np.round(mean_lab,0), np.round(upper_lab,0)
 261 |     
 262 |     if no_to_fit_curve==2:
 263 |         if i == 1:
 264 |             upper_lower_bounds.append([-2.105,-2.085])
 265 |             lower_lab, upper_lab = -2.1,-2.09
 266 |         if i == 2:
 267 |             upper_lower_bounds.append([-0.012,0.002])
 268 |             lower_lab, upper_lab = -0.009, -0.001
 269 |         if i!=1:
 270 |             if i!=2:
 271 |                 upper_lower_bounds.append([lower_b,upper_b])
 272 |     if no_to_fit_curve==1:
 273 |         if i==2:
 274 |             upper_lower_bounds.append([-0.504,-0.487])
 275 |         else:
 276 |             upper_lower_bounds.append([lower_b,upper_b])
 277 |     tick_labels.append([lower_lab,upper_lab])
 278 | print (tick_labels)
 279 | matplotlib.rc('xtick', labelsize=10) 
 280 | matplotlib.rc('ytick', labelsize=10) 
 281 | 
 282 | fig, axes = samples.plot_2d(['p%i' %i for i in range(1,nDims+1)])
 283 | 
 284 | fig.set_size_inches(8, 8)
 285 | for i in range(nDims):
 286 |     #set the axis labels for the y and x axis to the correct variable name
 287 |     axes.iloc[i,0].set_ylabel(param_names[i], fontsize=15,rotation=0,labelpad=22,horizontalalignment="center", verticalalignment="center",rotation_mode="anchor")
 288 |     axes.iloc[int(nDims-1),i].set_xlabel(param_names[i], fontsize=15,rotation=90)
 289 |     #set the tick labels for the subplots
 290 |     axes.iloc[i,0].set_yticks(tick_labels[i])
 291 |     axes.iloc[int(nDims-1),i].set_xticks(tick_labels[i])
 292 |     axes.iloc[i,0].set_yticklabels(tick_labels[i])
 293 |     axes.iloc[int(nDims-1),i].set_xticklabels(tick_labels[i],rotation=90)
 294 |     #set axis limits
 295 |     axes.iloc[int(nDims-1),i].set_xlim(upper_lower_bounds[i][0],upper_lower_bounds[i][1])
 296 |     axes.iloc[i,0].set_ylim(upper_lower_bounds[i][0],upper_lower_bounds[i][1])
 297 | 
 298 | #plot the true values for the params as a crossed set of lines
 299 | for i in range(nDims):
 300 |     for j in range(nDims):
 301 |         par_i_true = true_param_vals_for_all[i]
 302 |         par_j_true = true_param_vals_for_all[j]
 303 | 
 304 |         if j>i:
 305 |             pass
 306 |         else:
 307 |             if j==i:
 308 |                 axes.iloc[i,j].plot(par_j_true*np.ones(2),[-400,700],c="red",linestyle="--")
 309 |             else:
 310 |                 axes.iloc[i,j].plot([-400,700],par_i_true*np.ones(2),c="red",linestyle="--")
 311 |                 axes.iloc[i,j].plot(par_j_true*np.ones(2),[-400,700],c="red",linestyle="--")
 312 | plt.tight_layout()
 313 | #plt.show()
 314 | fig.savefig(root+'/sampled_posterior_with_cross.pdf')
 315 | 
 316 | 
 317 | 
 318 | #============================================================================
 319 | #plot the component map posterior
 320 | #load the posterior component maps
 321 | post_c1_mean, post_c2_mean = np.loadtxt(root+"/posterior_mean_for_comp_1"), np.loadtxt(root+"/posterior_mean_for_comp_2")
 322 | post_c1_std, post_c2_std = np.loadtxt(root+"/posterior_std_for_comp_1"),np.loadtxt(root+"/posterior_std_for_comp_2")
 323 | post_c1_mean[(post_c1_mean<=0)]=float("NaN")
 324 | post_c2_mean[(post_c2_mean<=0)]=float("NaN")
 325 | fig=plt.figure(figsize=(6,8.5))
 326 | 
 327 | ax1 = plt.subplot2grid((4,2),(0,0))
 328 | #ax1.suptitle("Component 1")
 329 | plt.axes(ax1)
 330 | hp.mollview(true_c1,title="Component 1\ntrue",hold=True,notext=True,norm="log")#,min=-10,max=10)#np.max(c_map_1))
 331 | 
 332 | ax2 = plt.subplot2grid((4,2),(1,0))
 333 | plt.axes(ax2)
 334 | try:
 335 |     hp.mollview(post_c1_mean,title="posterior mean",hold=True,notext=True,norm="log")
 336 | except:
 337 |     hp.mollview(post_c1_mean,title="posterior mean",hold=True,notext=True)
 338 | ax3 = plt.subplot2grid((4,2),(2,0))
 339 | plt.axes(ax3)
 340 | hp.mollview(post_c1_std,title="posterior std",hold=True,notext=True,norm="log")
 341 | 
 342 | ax4 = plt.subplot2grid((4,2),(3,0))
 343 | plt.axes(ax4)
 344 | nm1=(post_c1_mean-true_c1)/post_c1_std
 345 | hp.mollview(nm1,title="norm resid",hold=True,notext=True,min=-5,max=5)
 346 | 
 347 | ax1b = plt.subplot2grid((4,2),(0,1))
 348 | #ax1b.suptitle("Component 2")
 349 | plt.axes(ax1b)
 350 | hp.mollview(true_c2,title="Component 2\ntrue",hold=True,notext=True,norm="log")#,min=-10,max=10)#np.max(c_map_1))
 351 | 
 352 | ax2b = plt.subplot2grid((4,2),(1,1))
 353 | plt.axes(ax2b)
 354 | try:
 355 |     hp.mollview(post_c2_mean,title="posterior mean",hold=True,notext=True,norm="log")
 356 | except:
 357 |     hp.mollview(post_c2_mean,title="posterior mean",hold=True,notext=True)
 358 | ax3b = plt.subplot2grid((4,2),(2,1))
 359 | plt.axes(ax3b)
 360 | hp.mollview(post_c2_std,title="posterior std",hold=True,notext=True,norm="log")
 361 | 
 362 | ax4b = plt.subplot2grid((4,2),(3,1))
 363 | plt.axes(ax4b)
 364 | nm2=(post_c2_mean-true_c2)/post_c2_std
 365 | hp.mollview(nm2,title="norm resid",hold=True,notext=True,min=-5,max=5)
 366 | 
 367 | plt.subplots_adjust(hspace=0.3)
 368 | plt.tight_layout()
 369 | plt.rcParams.update({'font.size':12})
 370 | plt.savefig(root+"/posterior_component_comparison.pdf")
 371 | 
 372 | nm1_for_hist=nm1[~np.isnan(nm1)]
 373 | nm2_for_hist=nm2[~np.isnan(nm2)]
 374 | print ("for component 1 norm resids: mean",np.nanmean(nm1),"std",np.nanstd(nm1))
 375 | print ("for component 2 norm resids: mean",np.nanmean(nm2),"std",np.nanstd(nm2))
 376 | 
 377 | fig = plt.figure(figsize=(8,5))
 378 | plt.subplot(1,2,1)
 379 | plt.hist(nm1_for_hist,bins=40)
 380 | plt.title("Component 1\n"+r"$\mu=$"+" "+str(np.round(np.nanmean(nm1),2))+" "+r"$\sigma=$"+" "+str(np.round(np.nanstd(nm1),2)))
 381 | plt.xlabel("norm resid")
 382 | plt.ylabel("occurrence")
 383 | plt.xlim(-4,4)
 384 | plt.subplot(1,2,2)
 385 | plt.hist(nm2_for_hist,bins=40)
 386 | plt.title("Component 2\n"+r"$\mu=$"+" "+str(np.round(np.nanmean(nm2),2))+" "+r"$\sigma=$"+" "+str(np.round(np.nanstd(nm2),2)))
 387 | plt.xlabel("norm resid")
 388 | plt.xlim(-4,4)
 389 | plt.tight_layout()
 390 | plt.savefig(root+"/posterior_hist_for_comps.pdf")
 391 | np.savetxt(root+"/comp1_norm_resids_for_hist",nm1_for_hist,delimiter=",")
 392 | np.savetxt(root+"/comp2_norm_resids_for_hist",nm2_for_hist,delimiter=",")
 393 | #plt.show()
 394 | 
 395 | #=========================================================================
 396 | #plot the posterior sky predictions
 397 | fig = plt.figure(figsize=(9,13.5))
 398 | test_freqs = fmaps_for_test
 399 | norm_resid_vals = []
 400 | for i in range(len(test_freqs)):
 401 |     f = test_freqs[i]
 402 | 
 403 |     post_sky = np.loadtxt(root+"/bayesian_pred_"+str(f)+"MHz")
 404 |     post_sky_err = np.loadtxt(root+"/bayesian_errs_"+str(f)+"MHz")
 405 | 
 406 |     true_sky = np.loadtxt(p+"/mock_data_file/mock_dataset_v4/true_sky_"+str(float(f)))
 407 | 
 408 |     true_sky[(true_sky==unobs_marker)]=float("NaN")
 409 | 
 410 |     norm_resids = (post_sky-true_sky)/post_sky_err
 411 |     if i == 0:
 412 | 
 413 |         ax1 = plt.subplot2grid((len(test_freqs),4),(i,0))
 414 |         plt.axes(ax1)
 415 |         hp.mollview(true_sky,title="true "+str(f)+" MHz",hold=True,notext=True,norm="log")
 416 | 
 417 |         ax2 = plt.subplot2grid((len(test_freqs),4),(i,1))
 418 |         plt.axes(ax2)
 419 |         hp.mollview(post_sky,title="post mean",hold=True,notext=True,norm="log")
 420 | 
 421 |         ax3 = plt.subplot2grid((len(test_freqs),4),(i,2))
 422 |         plt.axes(ax3)
 423 |         hp.mollview(post_sky_err,title="post std",hold=True,notext=True,norm="log")
 424 | 
 425 |         ax4 = plt.subplot2grid((len(test_freqs),4),(i,3))
 426 |         plt.axes(ax4)
 427 |         hp.mollview(norm_resids,title="norm resid",hold=True,notext=True,min=-5,max=5)
 428 |     else:
 429 | 
 430 |         ax1 = plt.subplot2grid((len(test_freqs),4),(i,0))
 431 |         plt.axes(ax1)
 432 |         hp.mollview(true_sky,title="true "+str(f)+" MHz",hold=True,notext=True,norm="log")
 433 | 
 434 |         ax2 = plt.subplot2grid((len(test_freqs),4),(i,1))
 435 |         plt.axes(ax2)
 436 |         hp.mollview(post_sky,title="post mean",hold=True,notext=True,norm="log")
 437 | 
 438 |         ax3 = plt.subplot2grid((len(test_freqs),4),(i,2))
 439 |         plt.axes(ax3)
 440 |         hp.mollview(post_sky_err,title="post std",hold=True,notext=True,norm="log")
 441 | 
 442 |         ax4 = plt.subplot2grid((len(test_freqs),4),(i,3))
 443 |         plt.axes(ax4)
 444 |         hp.mollview(norm_resids,title="norm resid",hold=True,notext=True,min=-5,max=5)
 445 | 
 446 |     norm_resid_vals.append(norm_resids)
 447 | plt.subplots_adjust(hspace=0.6)
 448 | plt.tight_layout()
 449 | plt.rcParams.update({'font.size':18})
 450 | plt.savefig(root+"/posterior_sky_map_comparison.pdf")
 451 | #plt.show()
 452 | 
 453 | norm_resid_vals=np.array(norm_resid_vals).flatten()
 454 | norm_resids_for_hist = norm_resid_vals[~np.isnan(norm_resid_vals)]
 455 | 
 456 | mean_norm_resid,std_norm_resid= np.nanmean(norm_resids_for_hist), np.nanstd(norm_resids_for_hist)
 457 | print ("mean norm_resid =",mean_norm_resid," std = ",std_norm_resid)
 458 | 
 459 | fig=plt.figure(figsize=(5.5,5))
 460 | plt.hist(norm_resids_for_hist,bins=50)
 461 | plt.ylabel("occurrence")
 462 | plt.xlabel("norm resid")
 463 | plt.title("Sky Posterior (all freqs)\n"+r"$\mu=$"+" "+str(np.round(mean_norm_resid,2))+" "+r"$\sigma=$"+" "+str(np.round(std_norm_resid,2)))
 464 | plt.savefig(root+"/posterior_sky_maps_norm_resid_histo.pdf")
 465 | 
 466 | np.savetxt(root+"/norm_resids_for_hist",norm_resids_for_hist,delimiter=",")
 467 | #=============================================================================
 468 | #plot the posterior of the TvsLST
 469 | test_freqs_for_TvsLST = [45.0,50.0,60.0,70.0,74.0,80.0,150.0,159.0,408.0]#[45.0,50.0,60.0,70.0,74.0,80.0]#,150.0,159.0,408.0]
 470 | test_LSTs = np.linspace(0,24,73)[:-1]
 471 | 
 472 | 
 473 | 
 474 | #fig=plt.figure(figsize=(12,8.5))
 475 | TvsLST_true =[]
 476 | TvsLST_uncal = []
 477 | TvsLST_post =[]
 478 | print ("========================================================")
 479 | for i in range(len(test_freqs_for_TvsLST)):
 480 |     true_TvsLST = np.loadtxt(p+"/mock_data_file/mock_dataset_v4/true_noisy_TvsLST_"+str(test_freqs_for_TvsLST[i]))
 481 |     uncal_dataset_TvsLST = np.loadtxt(p+"/mock_data_file/mock_dataset_v4/peturbed_TvsLST_"+str(test_freqs_for_TvsLST[i]))
 482 | 
 483 |     TvsLST_true.append(true_TvsLST)
 484 |     TvsLST_uncal.append(uncal_dataset_TvsLST)
 485 | 
 486 |     #plt.subplot(3,3,i+1)
 487 |     #plt.title(str(test_freqs_for_TvsLST[i])+" MHz")
 488 |     #plt.plot(test_LSTs,true_TvsLST,c="red",label="true")
 489 |     #plt.plot(test_LSTs,uncal_dataset_TvsLST,c="green",label="uncal dataset")
 490 | 
 491 | 
 492 |     #generate EDGES beams set for this freq
 493 |     beams = gen_beams_and_T_vs_LST_v2.gen_EDGES_beams_at_LSTs(test_freqs_for_TvsLST[i],test_LSTs,Nside)
 494 |     beams = (1/(4*np.pi))*beams*hp.nside2pixarea(nside=Nside)
 495 | 
 496 |     #load the posterior sky
 497 |     post_sky = np.loadtxt(root+"/bayesian_pred_"+str(test_freqs_for_TvsLST[i])+"MHz")
 498 |     #generate a TvsLST for the posterior sky
 499 |     post_mean_TvsLST = np.array([np.nansum(beams[:,i]*post_sky) for i in range(len(test_LSTs))])
 500 |     
 501 |     
 502 |     TvsLST_post.append(post_mean_TvsLST)
 503 |     print ("for freq:",test_freqs_for_TvsLST[i])
 504 |     print ("mean resid (for posterior):",np.mean(post_mean_TvsLST-true_TvsLST))
 505 |     print ("rms resid (uncal):",np.sqrt(np.mean((uncal_dataset_TvsLST-true_TvsLST)**2)))
 506 |     print ("rms resid (posterior):",np.sqrt(np.mean((post_mean_TvsLST-true_TvsLST)**2)))
 507 |     #plt.plot(test_LSTs,post_mean_TvsLST,c="blue",linestyle="--",label="posterior mean")
 508 | 
 509 |     #plt.legend(loc="upper left")
 510 | #plt.tight_layout()
 511 | #plt.savefig(root+"/posterior_TvsLST_comp_all_freqs.png")
 512 | #plt.show()
 513 | plt.rcParams.update({'font.size':14})
 514 | fig=plt.figure(figsize=(8.5,13))
 515 | 
 516 | #row 1
 517 | ax1 = plt.subplot2grid((9,3),(0,0),rowspan=2)
 518 | ax1.plot(test_LSTs,TvsLST_true[0],label="true",c="red",linewidth=2)
 519 | ax1.plot(test_LSTs,TvsLST_uncal[0],label="uncal data",c="green",linewidth=2)
 520 | ax1.plot(test_LSTs,TvsLST_post[0],label="post mean",linewidth=2,linestyle="--",c="blue")
 521 | ax1.legend(loc="upper left")
 522 | ax1.set_title(r"45MHz")
 523 | ax1.set_ylabel(r"Ant Temp (K)")
 524 | 
 525 | ax2 = plt.subplot2grid((9,3),(2,0),rowspan=1,sharex=ax1)
 526 | #ax2.plot(test_LSTs,TvsLST_true[0]-TvsLST_uncal[0],label="true - uncal data",c="green",linewidth=2)
 527 | ax2.plot(test_LSTs,-TvsLST_true[0]+TvsLST_post[0],label="true - post mean",linewidth=2,c="blue")
 528 | #ax2.legend(loc="upper left")
 529 | ax2.set_ylabel(r"$\Delta$ Ant"+"\nTemp (K)")
 530 | 
 531 | ax1b = plt.subplot2grid((9,3),(0,1),rowspan=2)
 532 | ax1b.plot(test_LSTs,TvsLST_true[1],label="true",c="red",linewidth=2)
 533 | ax1b.plot(test_LSTs,TvsLST_uncal[1],label="uncal data",c="green",linewidth=2)
 534 | ax1b.plot(test_LSTs,TvsLST_post[1],label="post mean",linewidth=2,linestyle="--",c="blue")
 535 | #ax1b.legend(loc="upper left")
 536 | ax1b.set_title(r"50MHz")
 537 | #ax1b.set_ylabel(r"Ant Temp (K)")
 538 | 
 539 | ax2b = plt.subplot2grid((9,3),(2,1),rowspan=1,sharex=ax1b)
 540 | #ax2b.plot(test_LSTs,TvsLST_true[1]-TvsLST_uncal[1],label="true - uncal data",c="green",linewidth=2)
 541 | ax2b.plot(test_LSTs,-TvsLST_true[1]+TvsLST_post[1],label="true - post mean",linewidth=2,c="blue")
 542 | #ax2b.legend(loc="upper left")
 543 | #ax2b.set_ylabel(r"$\Delta$ Ant Temp (K)")
 544 | 
 545 | ax1c = plt.subplot2grid((9,3),(0,2),rowspan=2)
 546 | ax1c.plot(test_LSTs,TvsLST_true[2],label="true",c="red",linewidth=2)
 547 | ax1c.plot(test_LSTs,TvsLST_uncal[2],label="uncal data",c="green",linewidth=2)
 548 | ax1c.plot(test_LSTs,TvsLST_post[2],label="post mean",linewidth=2,linestyle="--",c="blue")
 549 | #ax1c.legend(loc="upper left")
 550 | ax1c.set_title(r"60MHz")
 551 | #ax1c.set_ylabel(r"Ant Temp (K)")
 552 | 
 553 | ax2c = plt.subplot2grid((9,3),(2,2),rowspan=1,sharex=ax1c)
 554 | #ax2c.plot(test_LSTs,TvsLST_true[2]-TvsLST_uncal[2],label="true - uncal data",c="green",linewidth=2)
 555 | ax2c.plot(test_LSTs,-TvsLST_true[2]+TvsLST_post[2],label="true - post mean",linewidth=2,c="blue")
 556 | #ax2c.legend(loc="upper left")
 557 | #ax2c.set_ylabel(r"$\Delta$ Ant Temp (K)")
 558 | 
 559 | #row 2
 560 | ax1 = plt.subplot2grid((9,3),(3,0),rowspan=2)
 561 | ax1.plot(test_LSTs,TvsLST_true[3],label="true",c="red",linewidth=2)
 562 | ax1.plot(test_LSTs,TvsLST_uncal[3],label="uncal data",c="green",linewidth=2)
 563 | ax1.plot(test_LSTs,TvsLST_post[3],label="post mean",linewidth=2,linestyle="--",c="blue")
 564 | #ax1.legend(loc="upper left")
 565 | ax1.set_title(r"70MHz")
 566 | ax1.set_ylabel(r"Ant Temp (K)")
 567 | 
 568 | ax2 = plt.subplot2grid((9,3),(5,0),rowspan=1,sharex=ax1)
 569 | #ax2.plot(test_LSTs,TvsLST_true[3]-TvsLST_uncal[3],label="true - uncal data",c="green",linewidth=2)
 570 | ax2.plot(test_LSTs,-TvsLST_true[3]+TvsLST_post[3],label="true - post mean",linewidth=2,c="blue")
 571 | #ax2.legend(loc="upper left")
 572 | ax2.set_ylabel(r"$\Delta$ Ant"+"\nTemp (K)")
 573 | 
 574 | ax1b = plt.subplot2grid((9,3),(3,1),rowspan=2)
 575 | ax1b.plot(test_LSTs,TvsLST_true[4],label="true",c="red",linewidth=2)
 576 | ax1b.plot(test_LSTs,TvsLST_uncal[4],label="uncal data",c="green",linewidth=2)
 577 | ax1b.plot(test_LSTs,TvsLST_post[4],label="post mean",linewidth=2,linestyle="--",c="blue")
 578 | #ax1b.legend(loc="upper left")
 579 | ax1b.set_title(r"74MHz")
 580 | #ax1b.set_ylabel(r"Ant Temp (K)")
 581 | 
 582 | ax2b = plt.subplot2grid((9,3),(5,1),rowspan=1,sharex=ax1b)
 583 | #ax2b.plot(test_LSTs,TvsLST_true[4]-TvsLST_uncal[4],label="true - uncal data",c="green",linewidth=2)
 584 | ax2b.plot(test_LSTs,-TvsLST_true[4]+TvsLST_post[4],label="true - post mean",linewidth=2,c="blue")
 585 | #ax2b.legend(loc="upper left")
 586 | #ax2b.set_ylabel(r"$\Delta$ Ant Temp (K)")
 587 | 
 588 | ax1c = plt.subplot2grid((9,3),(3,2),rowspan=2)
 589 | ax1c.plot(test_LSTs,TvsLST_true[5],label="true",c="red",linewidth=2)
 590 | ax1c.plot(test_LSTs,TvsLST_uncal[5],label="uncal data",c="green",linewidth=2)
 591 | ax1c.plot(test_LSTs,TvsLST_post[5],label="post mean",linewidth=2,linestyle="--",c="blue")
 592 | #ax1c.legend(loc="upper left")
 593 | ax1c.set_title(r"80MHz")
 594 | #ax1c.set_ylabel(r"Ant Temp (K)")
 595 | 
 596 | ax2c = plt.subplot2grid((9,3),(5,2),rowspan=1,sharex=ax1c)
 597 | #ax2c.plot(test_LSTs,TvsLST_true[5]-TvsLST_uncal[5],label="true - uncal data",c="green",linewidth=2)
 598 | ax2c.plot(test_LSTs,-TvsLST_true[5]+TvsLST_post[5],label="true - post mean",linewidth=2,c="blue")
 599 | #ax2c.legend(loc="upper left")
 600 | #ax2c.set_ylabel(r"$\Delta$ Ant Temp (K)")
 601 | 
 602 | #row 3
 603 | ax1 = plt.subplot2grid((9,3),(6,0),rowspan=2)
 604 | ax1.plot(test_LSTs,TvsLST_true[6],label="true",c="red",linewidth=2)
 605 | ax1.plot(test_LSTs,TvsLST_uncal[6],label="uncal data",c="green",linewidth=2)
 606 | ax1.plot(test_LSTs,TvsLST_post[6],label="post mean",linewidth=2,linestyle="--",c="blue")
 607 | #ax1.legend(loc="upper left")
 608 | ax1.set_title(r"150MHz")
 609 | ax1.set_ylabel(r"Ant Temp (K)")
 610 | 
 611 | ax2 = plt.subplot2grid((9,3),(8,0),rowspan=1,sharex=ax1)
 612 | #ax2.plot(test_LSTs,TvsLST_true[6]-TvsLST_uncal[6],label="true - uncal data",c="green",linewidth=2)
 613 | ax2.plot(test_LSTs,-TvsLST_true[6]+TvsLST_post[6],label="true - post mean",linewidth=2,c="blue")
 614 | #ax2.legend(loc="upper left")
 615 | ax2.set_ylabel(r"$\Delta$ Ant"+"\nTemp (K)")
 616 | ax2.set_xlabel(r"LST (hours)")
 617 | 
 618 | ax1b = plt.subplot2grid((9,3),(6,1),rowspan=2)
 619 | ax1b.plot(test_LSTs,TvsLST_true[7],label="true",c="red",linewidth=2)
 620 | ax1b.plot(test_LSTs,TvsLST_uncal[7],label="uncal data",c="green",linewidth=2)
 621 | ax1b.plot(test_LSTs,TvsLST_post[7],label="post mean",linewidth=2,linestyle="--",c="blue")
 622 | #ax1b.legend(loc="upper left")
 623 | ax1b.set_title(r"159MHz")
 624 | #ax1b.set_ylabel(r"Ant Temp (K)")
 625 | 
 626 | ax2b = plt.subplot2grid((9,3),(8,1),rowspan=1,sharex=ax1b)
 627 | #ax2b.plot(test_LSTs,TvsLST_true[7]-TvsLST_uncal[7],label="true - uncal data",c="green",linewidth=2)
 628 | ax2b.plot(test_LSTs,-TvsLST_true[7]+TvsLST_post[7],label="true - post mean",linewidth=2,c="blue")
 629 | #ax2b.legend(loc="upper left")
 630 | #ax2b.set_ylabel(r"$\Delta$ Ant Temp (K)")
 631 | ax2b.set_xlabel(r"LST (hours)")
 632 | 
 633 | ax1c = plt.subplot2grid((9,3),(6,2),rowspan=2)
 634 | ax1c.plot(test_LSTs,TvsLST_true[8],label="true",c="red",linewidth=2)
 635 | ax1c.plot(test_LSTs,TvsLST_uncal[8],label="uncal data",c="green",linewidth=2)
 636 | ax1c.plot(test_LSTs,TvsLST_post[8],label="post mean",linewidth=2,linestyle="--",c="blue")
 637 | #ax1c.legend(loc="upper left")
 638 | ax1c.set_title(r"408MHz")
 639 | #ax1c.set_ylabel(r"Ant Temp (K)")
 640 | 
 641 | ax2c = plt.subplot2grid((9,3),(8,2),rowspan=1,sharex=ax1c)
 642 | #ax2c.plot(test_LSTs,TvsLST_true[8]-TvsLST_uncal[8],label="true - uncal data",c="green",linewidth=2)
 643 | ax2c.plot(test_LSTs,-TvsLST_true[8]+TvsLST_post[8],label="true - post mean",linewidth=2,c="blue")
 644 | #ax2c.legend(loc="upper left")
 645 | #ax2c.set_ylabel(r"$\Delta$ Ant Temp (K)")
 646 | ax2c.set_xlabel(r"LST (hours)")
 647 | 
 648 | plt.subplots_adjust(hspace=0.97)
 649 | plt.subplots_adjust(wspace=0.28)
 650 | plt.savefig(root+"/posterior_TvsLST_comp_all_freqs_v2.pdf")
 651 | 
 652 | #plt.show()
 653 | 
 654 | 
 655 | #========================================================
 656 | #plot the posterior of the spectra
 657 | 
 658 | #define the spectra functional form
 659 | def spec_func(fs,params):
 660 |     #print (params)
 661 |     #t = (fs/f0)**(params[0]+params[1]*np.log(fs/f0))
 662 |     t = LogParabola1D(1,f0,-params[0],-params[1])(fs)
 663 |     #print (t)
 664 |     #plt.plot(fs,t)
 665 |     #plt.show()
 666 |     return t
 667 | 
 668 | #select the samples with non zero weight
 669 | sammps_set2 = nested_samples[(weights!=0),:]
 670 | nested_samples = sammps_set2
 671 | print ("no of samps with non-zero weight is:",nested_samples.shape)
 672 | #plot a posterior of the spectra for each component
 673 | print ("plotting posterior of the spectra")
 674 | plt.rcParams['font.size'] = 14
 675 | fig = plt.figure(figsize=(8,4))
 676 | plt_freqs = np.logspace(1.653,2.7,200)
 677 | for i in range(no_of_comps):
 678 |     spec_pars = np.empty(shape=(nested_samples.shape[0],2))
 679 |     if no_to_fit_curve == no_of_comps:
 680 |         spec_pars[:,0] = nested_samples[:,i]
 681 |         spec_pars[:,1] = nested_samples[:,no_of_comps+i]
 682 |     else:
 683 |         if no_to_fit_curve==1:
 684 |             if fit_curved[i] == True:
 685 |                 spec_pars[:,0] = nested_samples[:,i]
 686 |                 spec_pars[:,1] = nested_samples[:,2]
 687 |             else:
 688 |                 spec_pars[:,0] = nested_samples[:,i]
 689 |                 spec_pars[:,1] = 0
 690 |         if no_to_fit_curve ==0:
 691 |             spec_pars[:,0] = nested_samples[:,i]
 692 |             spec_pars[:,1] = 0
 693 |     print ("unweighted mean of spec pars is:",np.mean(spec_pars,axis=0))
 694 | 
 695 |     plt.subplot(1,no_of_comps,i+1)
 696 |     #plot the posterior samples
 697 |     plot_lines(spec_func, plt_freqs, spec_pars)
 698 |     #plot the true spectra
 699 |     print ("true spectral params are:",true_specs[i])
 700 |     plt.plot(plt_freqs,spec_func(plt_freqs,true_specs[i]),linestyle="--",c="red",linewidth=1.5,label="true:"+r"$\beta=$"+str(true_specs[i][0])+"  "+r"$\gamma=$"+str(true_specs[i][1]))
 701 |     print ("mean post pars:",[mean_params[2*i],mean_params[2*i+1]])
 702 |     #plt.plot(plt_freqs,spec_func(plt_freqs,[mean_params[2*i],mean_params[2*i+1]]),linestyle="--",c="purple",label="posterior mean")
 703 |     
 704 |     plt.legend(loc="upper right")
 705 |     plt.xlabel("freq (MHz)",fontsize=15)
 706 |     if i==0:
 707 |         plt.ylabel("spectral scaling",fontsize=15)
 708 |     plt.yticks(ticks=[0.01,0.1,1,10,100],labels=[0.01,0.1,1,10,100])
 709 |     plt.ylim(0.01,100)
 710 |     plt.xscale('log')
 711 |     plt.yscale("log")
 712 |     
 713 | 
 714 | 
 715 | 
 716 | #plt.ylabel("spectral scaling")
 717 | 
 718 | plt.tight_layout()
 719 | plt.savefig(root+"/spectral_posterior.pdf")
 720 | #plt.show()

```

`Bayesian-Global-Sky-Model-B-GSM-Paper-1/test_approximation2.py`:

```py
   1 | from random import sample
   2 | import matplotlib as mpl
   3 | #mpl.use('Agg')
   4 | from itertools import combinations_with_replacement
   5 | from itertools import product
   6 | import numpy as np
   7 | import healpy as hp
   8 | import matplotlib.pyplot as plt
   9 | from mpl_toolkits.axes_grid1 import make_axes_locatable
  10 | from matplotlib.colors import SymLogNorm
  11 | from matplotlib.colors import LogNorm
  12 | import os
  13 | import pandas
  14 | from numpy import pi, log, sqrt
  15 | import pix_by_pix_mk24_final2 as likelihood
  16 | 
  17 | import scipy.optimize as so
  18 | import math
  19 | import scipy
  20 | 
  21 | import matplotlib.cm as cm
  22 | from scipy.optimize import minimize
  23 | 
  24 | try:
  25 |     import pypolychord
  26 |     from pypolychord.settings import PolyChordSettings
  27 |     from pypolychord.priors import UniformPrior
  28 |     from pypolychord.priors import GaussianPrior
  29 | except:
  30 |     pass
  31 | try:
  32 |     from anesthetic import NestedSamples
  33 | except ImportError:
  34 |     pass
  35 | try:
  36 |     from anesthetic.weighted_pandas import WeightedDataFrame
  37 | except ImportError:
  38 |     pass
  39 | from scipy.optimize import minimize
  40 | import gen_EDGES_beams as gen_beams_and_T_vs_LST_v2
  41 | from line_profiler import LineProfiler
  42 | import matplotlib
  43 | from astropy.modeling.powerlaws import LogParabola1D
  44 | from csv import writer
  45 | 
  46 | #RUN PARAMS
  47 | #================================================================================================
  48 | save_no = 50
  49 | chunksize=100
  50 | if save_no!=1:
  51 |     how_many_post_points = [-save_no*chunksize,-(save_no-1)*chunksize] #which posterior points to use
  52 | else:
  53 |     how_many_post_points = [-save_no*chunksize,None]
  54 | print ("running for the posteior points between indexes:",how_many_post_points)
  55 | how_many_maps = 10 #for each posterior point how many component maps will we use
  56 | 
  57 | #MODEL PARAMS
  58 | #================================================================================================
  59 | #declare a random seed
  60 | 
  61 | use_perturbed_dataset = True #do we want the input dataset to have calibration errors
  62 | 
  63 | #===================================================================#
  64 | #| SET PARAMS FOR THE SIMULATED DATA AND NESTED SAMPLING
  65 | Max_Nside=32 #the Nside at which to generate the set of maps
  66 | Max_m = hp.nside2npix(Max_Nside)
  67 | no_of_comps = 2
  68 | fit_curved=[True,True]
  69 | 
  70 | no_to_fit_curve = np.sum(fit_curved)
  71 | rezero_prior_std=2000
  72 | 
  73 | #params for the prior on the spectra
  74 | spec_min, spec_max = -3.5,1 #the range for the prior on the spectral indexes
  75 | curvature_mean, curvature_std = 0,1 #the range for the prior on the spectral index curvature parameter
  76 | 
  77 | #params for fitting the reference frequency 
  78 | 
  79 | fixed_f0 = 150 #if you dont fit a seperate reference freq for each comp then we fix f0 to this value
  80 | 
  81 | 
  82 | #params for the prior on the true maps
  83 | map_prior_variance_spec_index = -2.8
  84 | map_prior_variance_f0 = 408#fixed_f0
  85 | map_prior_std=300
  86 | 
  87 | calibrate = True
  88 | calibrate_all_but_45_150 = False#True #calibrate all maps in the dataset but the 45 and 150 MHz maps
  89 | calibrate_all = True#False #calibrate every map in the dataset
  90 | 
  91 | 
  92 | use_equal_spaced_LSTs = True
  93 | fit_haslam_noise = False
  94 | subtract_CMB = 0#-2.726
  95 | print_vals_as_calc = False
  96 | 
  97 | 
  98 | 
  99 | reject_criterion = 1e-3#None #how close can two spectral idexes be in value before being rejected
 100 | cond_no_threshold =1e+9
 101 | 
 102 | 
 103 | 
 104 | 
 105 | test_LSTs = np.linspace(0,24,73)[:-1]#np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])#np.array([0,2,4,6,8,10,12,14,16,18,20,22]) #the LSTs in hours at which we will make comparison between the mean sky and EDGES for likelihood calls
 106 | #test_freqs = [47.5,75]#,250,300,350]#[45,50.005,59.985,70.007,73.931,79.960,100,125,150,159,200,408]
 107 | #test_freqs = [100,140,200]
 108 | #test_freqs = [250,300,350]
 109 | 
 110 | test_freqs = [45.0,50.0,60.0,70,74,80,150,159,408]
 111 | #test_freqs = [70.0,74.0,80.0]
 112 | #test_freqs = [150.0,159.0,408.0]
 113 | 
 114 | 
 115 | 
 116 | unobs_marker = -32768
 117 | 
 118 | nlive = 500#2500*no_of_comps
 119 | 
 120 | precision_criterion = 1e-3
 121 | 
 122 | 
 123 | f0=150 #ref freq used for some fitting of spec indexes for plots (not used during any model fitting)
 124 | 
 125 | freqs = np.array([45.0,50.0,60.0,70.0,74.0,80.0,150.0,159.0,408.0]) #the frequencies in MHz of maps used to generate the model
 126 | #specify what to fit for each map
 127 | if calibrate_all == True:
 128 |     freqs_to_calibrate = np.array([True,True,True,True,True,True,True,True,True]) #calibrate all the maps except the 45 and 150 MHz
 129 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 130 | 
 131 | if calibrate_all_but_45_150 == True:
 132 |     freqs_to_calibrate = np.array([False,True,True,True,True,True,False,True,True]) #calibrate all the maps except the 45 and 150 MHz
 133 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 134 | if calibrate==False:
 135 |     freqs_to_calibrate = np.array([False,False,False,False,False,False,False,False,False])#np.array([True,True,True,True,True,True,False,False])
 136 |     freqs_to_fit_noise = np.array([False,False,False,False,False,False,False,fit_haslam_noise,fit_haslam_noise])
 137 | 
 138 | freqs_for_T_v_LST_comp = np.array([40.0,45.0,50.0,55.0,60.0,65.0,70.0,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200.0])
 139 | 
 140 | 
 141 | main_label = "_petur:"+str(use_perturbed_dataset)+"_"+str(no_of_comps)+"_comp_cal:"+str(calibrate)+"_rezro_pri_std:"+str(rezero_prior_std)+"_CMB="+str(subtract_CMB)+"_map_pri_std:"+str(map_prior_std)+"_mu:0_map_pri_std_spec_ind="+str(map_prior_variance_spec_index)+"_map_pri_f0="+str(map_prior_variance_f0)+"_cond_no_thres="+str(np.round(np.log10(cond_no_threshold),1))+"_crv_N_std="+str(curvature_std)+"_spec="+str(spec_min)+"_to:"+str(spec_max)#+"_rej_crit="+str(reject_criterion)#+"_nlive="+str(nlive)+"_nrept="+str(nrepeat)+"_precision_criterion="+str(precision_criterion)
 142 | 
 143 | if use_equal_spaced_LSTs==True:
 144 |     #LSTs_for_comparison = np.array([2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])#np.array([0,2,4,6,8,10,12,14,16,18,20,22])#np.array([2.5,18]) #the LSTs in hours at which we will make comparison between the mean sky and EDGES for likelihood calls
 145 |     LSTs_for_comparison = np.linspace(0,24,73)[:-1]
 146 |     print (LSTs_for_comparison)
 147 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 148 |     if calibrate_all==True:
 149 |         #root="uni_EDGES_v4_data_mk24_no_of_curved:"+str(no_to_fit_curve)+"_cal_all_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 150 |         root="uni_EDGES_v4_data_mk24_no_of_curved:"+str(no_to_fit_curve)+"_cal_all_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 151 |     
 152 |     else:
 153 |         if calibrate_all_but_45_150==True:
 154 |             root="uni_EDGES_v4_dat_mk24_no_curve:"+str(no_to_fit_curve)+"_no_cal_45_150_f0="+str(fixed_f0)+main_label
 155 | else:
 156 |     #LSTs_for_comparison = np.array([0,1,2,3,4,5,6,7,8,9,10,11,11.25,11.5,11.75,12,12.25,12.5,12.75,13,13.25,13.5,13.75,14,14.25,14.5,14.75,15,15.25,15.5,15.75,16,16.25,16.5,16.75,17,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18,18.25,18.5,18.75,19,19.25,19.5,19.75,20,20.25,20.5,20.75,21,21.25,21.5,21.75,22,22.25,22.5,22.75,23,23.25,23.5,23.75])
 157 |     LSTs_for_comparison = np.array([0,2,4,6,8,10,12,14,15,15.5,15.75,16,16.25,16.5,16.75,17,17.25,17.5,17.75,18,18.25,18.5,18.75,19,19.25,19.5,20,21,22])
 158 |     print (LSTs_for_comparison)
 159 |     print ("no of LSTs is:",len(LSTs_for_comparison))
 160 |     
 161 |     root="mk24_extra_uneq_LSTs:"+str(len(LSTs_for_comparison))+"_fixed_f0="+str(fixed_f0)+main_label#"very_unequal_LST_lots_freq_vSTRG_BIAS"+main_label#"real_data_mk19_EDGES_"+main_label
 162 |     
 163 | 
 164 | #specify what to fit for each map
 165 | 
 166 | 
 167 | #set the prior for the noise (on the Haslam map)
 168 | noise_prior_lower, noise_prior_upper = np.array([0.01]),np.array([50])
 169 | 
 170 | print ("noise prior is from:",noise_prior_lower,"to",noise_prior_upper,"Kelvin")
 171 | 
 172 | n_spec_pars = 3*no_of_comps #the number of parameters for the spectra (for each comp we have: break_freq, spec_index1, spec_index2)
 173 | 
 174 | 
 175 | nv=len(freqs) #the number of freqs that have maps
 176 | no_of_fitted_noise = np.sum(freqs_to_fit_noise) #the number of freqs at which we fit noise level
 177 | no_of_calibrated = np.sum(freqs_to_calibrate) #the number of freqs at which we fit the calibration
 178 | 
 179 | if np.sum(freqs_to_calibrate)!=0:
 180 |     zero_lev_prior_std = rezero_prior_std*np.ones(np.sum(freqs_to_calibrate))#200*((np.array(freqs)[freqs_to_calibrate]/100)**-2.5)
 181 |     zero_lev_prior_means = np.zeros(np.sum(freqs_to_calibrate))
 182 | 
 183 | 
 184 |     #set the prior params for the scale corrections
 185 |     scale_prior_lower=0.85
 186 |     scale_prior_upper=1.25
 187 |     
 188 | 
 189 |     print ("zero level prior is gauss with mean 0K, stds (Kelvin):")
 190 |     print (zero_lev_prior_std)
 191 |     print ("temp scale prior is uniform from:",scale_prior_lower,"to",scale_prior_upper)
 192 | 
 193 | 
 194 | 
 195 | 
 196 | 
 197 | 
 198 | #CREATE A DIR TO STORE RESULTS
 199 | #====================================================================#
 200 | 
 201 | 
 202 | #make a dir to store the results
 203 | p=os.getcwd()+"/"
 204 | path = p+root+"/"
 205 | try:
 206 |     os.mkdir(root)
 207 | except:
 208 |     pass
 209 | #make a dir to store the results as we run
 210 | root2 = path+"running_results/"
 211 | try:
 212 |     os.mkdir(root2)
 213 | except:
 214 |     pass
 215 | 
 216 | #LOAD THE DATASET AND THE ERROR MAPS
 217 | #====================================================================#
 218 | 
 219 | obs_maps = []
 220 | inv_err_maps = []
 221 | data_err_maps = []
 222 | load_path = p+"mock_data_file/mock_dataset_v4/"
 223 | for i in range(len(freqs)):
 224 |     f=freqs[i]
 225 | 
 226 |     if use_perturbed_dataset==True:
 227 |         fname1 = "noisy_perturbed_sky_"+str(f)
 228 |         err_fname = "perturbed_err_map_"+str(f)
 229 |     else:
 230 |         fname1 = "noisy_sky_"+str(f)
 231 |         err_fname = "err_map_"+str(f)
 232 |     
 233 |     #fname2 = "noise_"+str(f)
 234 | 
 235 |     if freqs_to_fit_noise[i]==False:
 236 |         try:
 237 |             #m1, err_m = np.loadtxt(load_path+fname1), np.loadtxt(load_path+fname2)
 238 |             m1 = np.loadtxt(load_path+fname1)
 239 | 
 240 |             err_m = np.loadtxt(load_path+err_fname)
 241 |         except:
 242 |             print ("cant find the files for freq:",f)
 243 | 
 244 |         
 245 |         #mask out any pixels with negative temps
 246 |         bool_arr = m1<=0
 247 |         err_m[bool_arr] = unobs_marker
 248 |         m1[bool_arr] = unobs_marker
 249 | 
 250 |         
 251 | 
 252 |         inv_err_m = 1/err_m
 253 |         inv_err_m[(err_m==unobs_marker)] = 0
 254 | 
 255 |         m1[m1!=unobs_marker] = m1[m1!=unobs_marker]+subtract_CMB
 256 |         obs_maps.append(m1)
 257 |     
 258 |         inv_err_maps.append(inv_err_m)
 259 | 
 260 |         data_err_maps.append(err_m)
 261 |     else:
 262 |         m1 = np.loadtxt(load_path+fname1)
 263 |         m1[m1!=unobs_marker] = m1[m1!=unobs_marker]+subtract_CMB
 264 |         obs_maps.append(m1)
 265 |     
 266 | 
 267 | obs_maps=np.array(obs_maps)
 268 | inv_err_maps=np.array(inv_err_maps)
 269 | print ("dataset loaded")
 270 | 
 271 | #CREATE THE INVERSE NOISE MATRICES
 272 | #====================================================================#
 273 | #generate the inverse noise covariance matrix for each pixel
 274 | inv_noise_mats = np.empty(shape=(Max_m,len(freqs),len(freqs)))
 275 | for p in range(0,Max_m):
 276 |     inv_stds_for_pixel = np.zeros(len(freqs))
 277 |     inv_stds_for_pixel[~freqs_to_fit_noise] = inv_err_maps[:,p]
 278 |     #print (inv_stds_for_pixel)
 279 | 
 280 |     Np_inv = np.diag(inv_stds_for_pixel**2)
 281 |     #print (Np_inv)
 282 |     inv_noise_mats[p,:,:] = Np_inv
 283 | 
 284 | print ("max and min for the inv noise mats: ",np.max(inv_noise_mats),np.min(inv_noise_mats[(inv_noise_mats!=0)]))
 285 | 
 286 | print ("inverse noise matrices created")
 287 | #PLOT THE DATASET
 288 | #====================================================================#
 289 | fig = plt.figure(figsize=(12,16))
 290 | for i in range(len(freqs)):
 291 |     map_i = np.copy(obs_maps[i])
 292 |     ax = plt.subplot(5,3,int(i+1))
 293 |     map_i[(map_i==unobs_marker)]=float("NaN")
 294 |     plt.axes(ax)
 295 |     hp.mollview(map_i,title="Synthetic Data Freq="+str(freqs[i]),hold=True,notext=True,norm="log")
 296 | 
 297 | 
 298 | 
 299 | plt.savefig(path+"sky_maps_for_dataset_for_plt")
 300 | #plt.show()
 301 | fig = plt.figure(figsize=(12,16))
 302 | for i in range(len(freqs)):
 303 |     map_i = np.copy(inv_err_maps[i])
 304 |     ax = plt.subplot(5,3,int(i+1))
 305 |     map_i[(map_i==0)]=float("NaN")
 306 |     plt.axes(ax)
 307 |     hp.mollview(1/map_i,title="input errs freq="+str(freqs[i]),hold=True,notext=True,norm="log")
 308 | 
 309 | 
 310 | 
 311 | plt.savefig(path+"err_maps_for_dataset_for_plt")
 312 | freqs=np.array(freqs)
 313 | #plot the priors and the data
 314 | #====================================================================
 315 | log_mean_temps = []
 316 | mean_temps = []
 317 | for i in range(len(freqs)):
 318 |     the_map = obs_maps[i]
 319 |     mean = np.mean(the_map[(the_map!=unobs_marker)])
 320 |     log_mean_temps.append(np.log(mean))
 321 |     mean_temps.append(mean)
 322 | log_mean_temps = np.array(log_mean_temps)
 323 | mean_temps = np.array(mean_temps)
 324 | 
 325 | log_freqs = np.log(freqs/f0)
 326 | fun = lambda x: np.nansum((log_mean_temps - x[0]*log_freqs -x[1])**2)
 327 | res = minimize(fun,[-2.15,np.log(16)])
 328 | print (res)
 329 | 
 330 | fig = plt.figure(figsize=(6,6))#figsize=(12,16))
 331 | #the fitted powerlaw
 332 | targ = np.exp(res.x[1])*((freqs/f0)**res.x[0])
 333 | 
 334 | ax1 = plt.subplot(1,1,1)
 335 | ax1.plot(freqs,targ,c="red",label="fitted power law")
 336 | ax1.scatter(freqs,mean_temps,label="data")
 337 | ax1.set_title("map mean temps")
 338 | ax1.legend(loc="upper right")
 339 | ax1.set_xscale("log")
 340 | ax1.set_yscale("log")
 341 | plt.savefig(path+"/input_map_for_plt_mean_temps.png")
 342 | 
 343 | 
 344 | #generate a set of pre rotated EDGES beams at each of the frequencies that we want to compare the model to EDGES for
 345 | #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 346 | print ("generating EDGES beams")
 347 | EDGES_beams = np.empty(shape=(obs_maps.shape[1],len(freqs_for_T_v_LST_comp),len(LSTs_for_comparison)))
 348 | j=0
 349 | for f in freqs_for_T_v_LST_comp:
 350 |     print ("gen EDGES beams for freq:",f)
 351 |     beams_at_LSTS = gen_beams_and_T_vs_LST_v2.gen_EDGES_beams_at_LSTs(f,LSTs_for_comparison,Max_Nside)
 352 |     EDGES_beams[:,j,:] = beams_at_LSTS
 353 |     #for i in range(len(LSTs_for_comparison)):
 354 |     #    beam = beams_at_LSTS[:,i]
 355 |         #print (beam)
 356 |         #print (beam.shape)
 357 |     #    hp.mollview(beam.flatten(),title="EDGES beam (galactic coords) for freq: "+str(f)+" for LST: "+str(LSTs_for_comparison[i]))
 358 |     #    plt.savefig(path+"EDGES_beam_galactic_coords_for_freq:"+str(f)+"_for_LST:"+str(LSTs_for_comparison[i])+".png")
 359 |     #    plt.close("all")
 360 |     j+=1
 361 | 
 362 | #Load the mock EDGES T vs LST plots (these are produced by convolving the mock sky map (before any pertubation) with the model beam at that freq)
 363 | #===============================================================================================
 364 | EDGES_temps_at_calib_LSTs_and_freqs = []
 365 | EDGES_errs = []
 366 | print ("generating the EDGES T vs LST traces")
 367 | for f in freqs_for_T_v_LST_comp:
 368 |    # if f<100:
 369 |    #     T_vs_LST, T_vs_LST_errs = gen_beams_and_T_vs_LST_v2.gen_EDGES_low_T_LST_trace(f,LSTs_for_comparison)
 370 |    #     EDGES_temps_at_calib_LSTs_and_freqs.append(T_vs_LST)
 371 |    #     EDGES_errs.append(T_vs_LST_errs)
 372 |    #     fig=plt.figure()
 373 |    #     plt.errorbar(LSTs_for_comparison,T_vs_LST,T_vs_LST_errs)
 374 |    #     plt.savefig(path+"EDGES_T_vs_LST_freq="+str(f)+".png")
 375 |    # if f>=100:
 376 |    #     T_vs_LST, T_vs_LST_errs = gen_beams_and_T_vs_LST_v2.gen_EDGES_high_T_LST_trace(f,LSTs_for_comparison)
 377 |    #     EDGES_temps_at_calib_LSTs_and_freqs.append(T_vs_LST)
 378 |    #     EDGES_errs.append(T_vs_LST_errs)
 379 |    #     fig=plt.figure()
 380 |    #     plt.errorbar(LSTs_for_comparison,T_vs_LST,T_vs_LST_errs)
 381 |    #     plt.savefig(path+"EDGES_T_vs_LST_freq="+str(f)+".png")
 382 |     T_vs_LST = np.loadtxt(load_path+"true_noisy_TvsLST_"+str(f))
 383 |     T_vs_LST_errs = np.loadtxt(load_path+"TvsLST_errs_"+str(f))
 384 |     print ("TvsLST for freq:",f)
 385 |     print (T_vs_LST)
 386 |     print (T_vs_LST_errs)
 387 |     EDGES_temps_at_calib_LSTs_and_freqs.append(T_vs_LST)
 388 |     EDGES_errs.append(T_vs_LST_errs)
 389 |     fig=plt.figure()
 390 |     plt.errorbar(LSTs_for_comparison,T_vs_LST,T_vs_LST_errs)
 391 |     plt.savefig(path+"EDGES_T_vs_LST_freq="+str(f)+".png")
 392 | plt.close("all")
 393 | 
 394 | EDGES_temps_at_calib_LSTs_and_freqs = np.array(EDGES_temps_at_calib_LSTs_and_freqs)
 395 | EDGES_errs = np.array(EDGES_errs)
 396 | #generate the EDGES noise covar mats for each freq assuming noise for each LST is independent of the other LSTs
 397 | EDGES_inv_noise_mats = []
 398 | for i in range(len(test_freqs)):
 399 |     inv_cov = np.diag(1/EDGES_errs[i,:]**2)
 400 |     EDGES_inv_noise_mats.append(inv_cov)
 401 | EDGES_inv_noise_mats = np.array(EDGES_inv_noise_mats)
 402 | #=====================================================================================================
 403 | #set up the likelihood function
 404 | #bayes_eval = likelihood.bayes_mod(obs_maps=obs_maps,obs_freqs=freqs,inv_noise_mats=inv_noise_mats,gaussian_prior_covar_mat=gaussian_prior_covar_matrix,gaussian_prior_mean=gaussian_prior_mean,no_of_comps=no_of_comps,f0=f0,un_obs_marker=unobs_marker)
 405 | 
 406 | #set up the likelihood function
 407 | bayes_eval = likelihood.bayes_mod(obs_maps=obs_maps,obs_freqs=freqs,inv_noise_mats=inv_noise_mats,EDGES_beams=EDGES_beams,EDGES_temps_at_calib_LSTs_and_freqs=EDGES_temps_at_calib_LSTs_and_freqs,EDGES_errs=EDGES_errs,EDGES_inv_noise_mats=EDGES_inv_noise_mats,freqs_for_T_v_LST_comp=freqs_for_T_v_LST_comp,LSTs_for_comparison=LSTs_for_comparison,no_of_comps=no_of_comps,save_root=root2,un_obs_marker=unobs_marker,map_prior_std=map_prior_std,map_prior_spec_index=map_prior_variance_spec_index,map_prior_f0=map_prior_variance_f0)
 408 | 
 409 | #test
 410 | print ("===================================")
 411 | print ("testing")
 412 | sample_map = bayes_eval.gen_comp_map_sample([150,-2.50,0,150,-0.5,-2.5,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1],freqs_to_fit_noise,freqs_to_calibrate)
 413 | for c in range(no_of_comps):
 414 |     m = sample_map[:,c,0]
 415 |     hp.mollview(m)
 416 |     plt.savefig("test_c="+str(c)+".png")
 417 | 
 418 | 
 419 | no_of_params_for_spec_mod = no_of_comps + no_to_fit_curve 
 420 | print ("we arn't fitting f0: f0=",fixed_f0," no of comps with curved spectra is",no_to_fit_curve," no of params for spectral model is",no_of_params_for_spec_mod)
 421 | 
 422 | 
 423 | #-------------NESTED SAMPLING PARAMS-------------
 424 | nDims =  int(no_of_params_for_spec_mod + no_of_fitted_noise + 2*no_of_calibrated)
 425 | nrepeat = 5*nDims #the nrepeat is set to 5 times the total number of pars that we fit
 426 | print ("no of dimensions for sampling region is:",nDims)
 427 | nDerived = 0 #we don't derive any parameters 
 428 | settings = PolyChordSettings(nDims, nDerived)
 429 | settings.file_root = root
 430 | settings.nlive = nlive
 431 | settings.nrepeats = nrepeat
 432 | settings.do_clustering = True
 433 | settings.read_resume = True
 434 | settings.write_resume = True
 435 | settings.maximise = False #find the maximum of the poseterior
 436 | settings.precision_criterion = precision_criterion
 437 | #------------------------------------------------
 438 | 
 439 | #find the mean spectral params and their standard deviations
 440 | samples = NestedSamples(root= settings.base_dir + '/' + settings.file_root)
 441 | 
 442 | mean_logZ = samples.logZ()#mean
 443 | std_logZ = samples.logZ(100).std()#100 posterior samples from estimate of log Z (use these to calculate standard deviation in Z)
 444 | 
 445 | print ("log Z",mean_logZ,"std",std_logZ)
 446 | 
 447 | def abs_temp_likelihood(comp_map_sample,spec_params):
 448 |     comp_map_samp[np.isnan(comp_map_samp)]=0
 449 |     comp_maps = np.zeros(shape=(bayes_eval.no_of_pixels,bayes_eval.no_of_comps,1))
 450 |     #print (np.sum(bayes_eval.pix_to_regard))
 451 |     comp_maps= comp_map_sample#[bayes_eval.pix_to_regard] = comp_map_sample
 452 | 
 453 |     spectral_mix_mat_for_test = bayes_eval.gen_A(spec_params,bayes_eval.freqs_for_T_v_LST_comp)#self.freqs_for_T_v_LST_comp[:,np.newaxis] ** spec_params[np.newaxis,:]
 454 |     #print ("spectral mix mat")
 455 |     #print (spectral_mix_mat_for_test)
 456 |     sky_preds = spectral_mix_mat_for_test @ comp_maps
 457 |     sky_preds[~bayes_eval.pix_to_regard] = float("NaN")
 458 | 
 459 |     convolved_sky_preds = sky_preds * bayes_eval.EDGES_beams
 460 | 
 461 |     #compute the integrated sky temp for each freq and LST in the freqs to compare
 462 |     integrated_skys = np.nansum(convolved_sky_preds,axis=0)
 463 | 
 464 |             
 465 |     EDGES_likelihood_t2s = ((integrated_skys-bayes_eval.EDGES_temps)/bayes_eval.EDGES_errs)**2
 466 |     EDGES_likelihood_t2 = -1*np.sum(EDGES_likelihood_t2s)#-1 * np.sum(diff_trans @ self.EDGES_inv_noise_mats @ diff)
 467 |     EDGES_int_sky_temps_log_likelihood = EDGES_likelihood_t2 + bayes_eval.EDGES_noise_det_term#(self.no_of_defined_pix/len(self.LSTs_for_comparison))*(EDGES_likelihood_t2 + self.EDGES_noise_det_term)
 468 |             
 469 | 
 470 |     return EDGES_int_sky_temps_log_likelihood
 471 | 
 472 | 
 473 | #load the samples from the marginal posterior
 474 | marg_post_samps = np.loadtxt(root+"/post_samples_marginal.csv",delimiter=",")
 475 | weights = np.loadtxt(root+"/post_samples_marginal_weights.csv",delimiter=",")
 476 | 
 477 | print (marg_post_samps)
 478 | #loop through the marginal posterior samps
 479 | abs_temp_likelihood_values = []
 480 | count = 0
 481 | for marg_samp in marg_post_samps[how_many_post_points[0]:how_many_post_points[1]]:
 482 |     #print (marg_samp)
 483 |     for i in range(how_many_maps):
 484 |         #generate a map sample
 485 |         comp_map_samp = bayes_eval.gen_comp_map_sample(marg_samp,freqs_to_fit_noise,freqs_to_calibrate)
 486 |         #print (i)
 487 |         #print (comp_map_samp)
 488 | 
 489 |         #compute the EDGES likelihood term for this comp map sample
 490 |         log_likelihood_value = abs_temp_likelihood(comp_map_samp,marg_samp[:6]) 
 491 |         #print ("log likelihood:",log_likelihood_value)
 492 |         abs_temp_likelihood_values.append(log_likelihood_value)
 493 |     count+=1
 494 |     print("percent done:",100*count/abs(chunksize))
 495 | 
 496 | abs_temp_likelihood_values = np.array(abs_temp_likelihood_values)
 497 | print ("mean likelihood val:",np.mean(abs_temp_likelihood_values))
 498 | print ("std:",np.std(abs_temp_likelihood_values))
 499 | np.savetxt(root+"/abs_temp_likelihood_values"+str(save_no),abs_temp_likelihood_values,delimiter=",")
 500 | plt.figure(figsize=(6,6))
 501 | plt.rcParams['font.size'] = 14
 502 | plt.hist(abs_temp_likelihood_values,bins=50)
 503 | plt.xlabel(r"log(P(E|M,S))")
 504 | plt.ylabel("occurrence")
 505 | plt.savefig(root+"/histo_of_abs_temp_likelihood_vals.pdf")
 506 | plt.show()

```