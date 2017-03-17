MOSiC : IDL package for IRIS spectral data reduction and Gaussian fitting


+++++++++++++++++++++
Installation
+++++++++++++++++++++

If you are using SolarSoft, you only need the following folders:

iris_reza, reza_lib, other_lib

put them in your idl path and you are done.

If you are not a SolarSoft user,
you have to unpack the whole package and add it to your idl path.


+++++++++++++++++++++
General remarks
+++++++++++++++++++++
I do not use the header information about the line position or dispersion.
The user interactively select lines and the program estimates the amount of the dispersion.
For our Gaussian fits, we use the MPFIT package.
Please note that at the moment, the programs do work on off-limb data but this is experimental.

Each sub-program checks if the input profile is empty before executing anything, so the program
does not crash in case of a dead column or like that. 

Analyzing NUV and FUV data, or every two wavelength band are completely independent.
A common block is used to pass the spectral dispersion to the functions fitting the
profiles.

Due to limited RAM size, I configured the program such that each set of lines should be called
separately. It means that the program reshape the full array first, and then read each raster spectra separately.
so things like

reorder, 'my_data_path'	 	     	     (this step should be executed first !)
analyze_iris, 'my_data_path', /do_cont, /do_mg, /do_si

in principle works but needs a lot of RAM.
A better approach is to call it like in the following order:

reorder, 'my_data_path'	 	     	     (this step should be executed first !)
analyze_iris, 'my_data_path', /do_cont       (this step should be executed second !)

analyze_iris, 'my_data_path', /do_mg, /do_h, /do_gauss

analyze_iris, 'my_data_path', /do_o1		triple Gaussian fit to the O I + two C I lines
	      		      			/do_o1 should be executed before any other FUV channel
						as we get the systematic velocity form the O I line
						
analyze_iris, 'my_data_path', /do_cl		single Gaussian fit
analyze_iris, 'my_data_path', /do_cii		penta Gaussian fit to the two C II and the Ni II lines
analyze_iris, 'my_data_path', /do_si, /do_1394	hexa Gaussian fit to Si IV 140.3, three O IV and the S IV lines

....


The temporal gradient due to orbital velocity is calculated separately for NUV and FUV spectra.
In NUV, one of the strong unblended photospheric lines is used after analyzing the NUV continuum data.
In FUV, we first take average profile of each slit position and then use the O I line center
to construct a temporal curve. 

The velocity maps, e.g. from O I or Cl I are similar to photospheric velocities.
The line width calculated in any program is the 1/e width which is needed to calculate the
non-thermal width. 

A common feature to all analysis programs is that at the very beginning, it analyze an average profile. 
The user can check the screen output and make sure that everything is alright. Also it is advisable to run a few
slit positions with visualization active to one can make sure that parameters are set correctly.


comments on the Gaussian fits:
--------------------------------------------
a) each parameter has a limited range,
b) output of the fit is checked and the fit is repeated if the result is not acceptable,
c) the degree of freedom increases gradually, i.e., we always start from a single Gaussian fit. 
  a higher SNR is required to run high-order fits.
d) quality of the fit is evaluated using the reduced chi-square. It has the chi-square distribution
   around unity which is an indicator of satisfactory fit. 


----------------------------------------------------------------
1) convert data to standard FITS = level 3
----------------------------------------------------------------
It is easier to have normal fits files rather than Rice-compressed files. 

$ tar xzvf iris_l2_20160509_131630_3684510004_raster.tar.gz

and inside idl:

ww = '/data/obs/iris/spots/mercyry2/'
infile = ['iris_l2_20160509_131630_3684510004_raster_t000_r00000.fits']
iris_make_fits_level3, ww + infile, /all, wdir=ww


----------------------------------------------------------------
for very large files, calibrate one by one so instead of /all:
----------------------------------------------------------------
iris_make_fits_level3, ww + infile, 0, wdir=ww
iris_make_fits_level3, ww + infile, 1, wdir=ww
iris_make_fits_level3, ww + infile, 2, wdir=ww


------------------------------------------------------------------------
2) flatten timeseries (skip this step for normal observations)
------------------------------------------------------------------------

iris_timeseri, '/path/to/iris_level3_file/'


it converts 4D data to 3D data (overwrite).

Please note that program perhaps crashes if you put more than one level 3 file in one directory. 

---------------------
3)  reorder it 
---------------------

reorder, '/path/to/data/file/'

we do so to have the wavelength as the first axis, slit as second, and scan as third.
it overwrites the existing calibrated file.



The core program is called "analyze_iris.pro" .

It calls a few ten subroutines.
Due to varying spectral range and dispersions in the data,
some numbers are hardwired. One need to be VERY CAREFUL to modify the numbers, e.g.
for the initial guess of line center position. 
In general, user does not need to modify this file or any other IDL program.

----------------------------------------------------------------
4) I suggest to reduce the data in the following order
----------------------------------------------------------------

$ idl

activate mosic in your idl PATH


> analyze_iris, '/full/path/to/level3_data/', /do_cont

			         This program first creates a quasi-continuum map in one Jpeg file
			         (by averaging some spectral bands),
			         so it is useful for a quicklook what is inside the data cube.
			         If the fits maps are separated (one file for each wavelength band) due to large file size,
				 then the program should be called like


			       	 The user then selects the spectral range in the NUV and FUV (by mouse click).
				 This selection is done once and then used in the other subroutines. 
				 The measured spectral dispersion will be printed and compared to
				 the official dispersion in the header. 
				 After that, the 280 nm range will be analyzed in detail: the program
				 calculates line parameters for two photospheric lines selected by the user.
				 I recommend to use a line which is deep, isolated, and not at the border of the spectral range.
			         Next it computes a real photospheric continuum map by analyzing the whole 280 nm spectra, 
				 and photospheric velocities. It calculates and applies  correction for
				 the satellite orbital velocity. User should select an order for the polynomial fit. 
 	       		         This step must be done (once) before any other step.
				 We used the line positions to correct for the orbital velocity. The two velocity maps
				 shown at the end of this step must have meaningful rad/blue structure of granulation.
				 In not, re-run this step and select a different photospheric line/range.

				 This program also creates a simple file to note the spectral dispersion of the data. 
				 It will be used in the later fitting procedures and functions.
				 In addition, it generates a Jpeg file of the average profile. 
				 
				 
 analyze_iris, '/full/path/to/level3_data/', /do_o1, /do_quiet
 
 	       		         O I   1355.598 A and Cl I   1354.288/1355.84 A 
 	       		         multi-line fitting of the Cl I and O I lines. These lines are often weak.
				 You should not be surprised that the program by itself find the O I line and make a Gaussian fit:
				 we have selected the spectral range in the previous step (above).
				 The maps of the 140 nm continuum, line position, amplitude, and widths get updated progressively,
				 if one keeps the visualization active.
				 /do_quiet blocks all visualizations. 
				 
				 
 analyze_iris, '/full/path/to/level3_data/', /do_mg, /do_h, /do_gauss, /do_quiet
 
 	       			 to analyze the Mg II h/k lines.
				 it performs all the calculations and returns a structured variable to the main program. 
 	       	       	         Its main task is emission peak analysis and a multi-Gaussian line fitting.
				 The program at first creates a big plot and marks all the important line position and ranges,
				 including K1 minimum range, K2v, K3, K2r, the same for H like
				 and also marks two photospheric lines.
				 The user can correct if the shown line positions are not satisfactory. This is often the case
				 close to the limb where the line is slightly broader than on the disk center.
				 the profile analysis program calculates several parameters for the emission peaks
				 (position, width, amplitude, asymmetry, WB, H1/2/3 intensity/velocity, bands, ..). 

				 Please note that the program always analyzes the two lines separately (first k, then h)
				 so if you want the Gaussian fit for both h and k lines, you have to activate it separately.
				 The Gaussian fit first fits a single Gaussian function to the core of the line (between k1 or h1 minima).
				 Then, it fits a double Gaussian which fits much better to the line and finally a triple Gaussian.
				 The runtime without the Gaussian fit is very short but including the Gaussian fit,
				 it takes several minutes per slit.
				 
				 
 analyze_iris, '/full/path/to/level3_data/', /do_si, /do_1394, /do_quiet
 
 	       			 Si IV  1397 and 1403 +  O IV line pair
				 It first plots the line profiles, and over plot a single and a double Gaussian fit,
				 so you can check that the line positions are correct.
 	       		       	 For each profile, the program runs as following:

				 first a single Gaussian fit to the Si IV line (4 params: continuum, line[center, amplitude,width]),
				 then a double Gaussian fit to the Si IV line (6 params: 2xline[center, amplitude,width]),  
				 after that a penta Gaussian fit to Si IV and four other lines:
				 O IV line at 1399.97 A,
				 O IV line at 1401.16 A,
				 Si IV line at 1404.79 A,
				 O IV line at 1404.82 A.
				 This fit has seven degrees of freedom:
				 Si IV 1403[center, amplitude,width], and four amplitudes for the other four lines.
				 widths of O IV is taken as 1.22 times width of Si IV. 
				 The line positions have a fixed offset w.r.t. Si IV 1403. 
				 The continuum level is fitted only in the first step and then is used as a known parameter.
				 In the next step, the program fits the O IV 1401 line with a single Gaussian fit while other parameters
				 are similar to the penta-Gaussian fit. Finally, the last step is repeated with a double Gaussian fit
				 for the Si IV 1403 line. In the last three steps, the program fits several small spectral lines 
				 to improve the overall fit quality and reduced the chi-square. These lines are fitted with no additional
				 degree of freedom (amplitude/line width are scaled from Si IV and O Iv lines). 
				 
				 Like all other Gaussian fittings in this IRIS reduction package, the fit takes
				 into account the photon noise in the data analysis via ir_error.pro considering
				 the conversion factor between photons and electrons in NUV and FUV channels.
				 
				 This guided fitting approach helps to make sure that the fit does not fail as
				 it has very few degrees of freedom (line positions/widths are hard-wired w.r.t. Si IV).
				 Also note that there is a range for each parameter so in worse case, it cannot return a value
				 outside this range.

				 The continuum in 140 nm range is calculated using PDF of a set of spectral points in two
				 continuum windows, and then calculating the median of a subset excluding the extreme values in the
				 upper and the lower limits. This avoids a lot of bad pixels but in case most of
				 continuum is damaged, it still can return a few hot pixels in the final continuum map.
				 There are variations along the scanning axis between successive columns. 
				 It can be clearly identified in FFT of big maps when it has a repeating wavelength of 4-5 pixels, and due to its
				 low amplitude of about 0.5 count, it is not seen in any other intensity map.
				 
				 
				 There is a hardwired threshold in the program to skip profiles with a maximum
				 amplitude below a certain level. 

				 For Si IV 1397, it is simple and the program runs only a single and double Gaussian fit,
				 like the procedure in C II lines. There, the amplitude is larger by a factor two so we have better
				 SNR but the line is blended.


 analyze_iris, 'spot_20160414/', /do_cii, /do_fast, /do_quiet
 	       			 pair of C II lines and the Ni II 133.52 nm line
				 At first, the program performs single and double Gaussian fits
				 to each C II line in a multi-Gaussian fitting scheme (together).
				 Then it includes the Ni II line as a single Gaussian fit
				 in a penta Gaussian fit with 13 free parameters.
				 
				 the single Gaussian fit: (4 params: continuum, line[center, amplitude,width]),
				 the double Gaussian fit: (5 params: 2xline[center, amplitude], one width),
				 the quad Gaussian fit: 2x double Gaussian fit
				 the penta Gaussian fit: quad Gaussian + line[center, amplitude,width] for Ni II
				 



analyze_iris, 'spot_20160414/', /do_cl
 	       			 only Cl I  1351.657 A  (single line).






in each step, a separate Save file is generated in the directory of the level3 data.


---------------------------------------
5) Evaluate results
---------------------------------------
use check_iris_si_fit.pro  to check quality of the Gaussian fit. There are similar programs for other bands:

   > check_iris_si_fit, '/path/to/save/file/'
   
   check_iris_si1_fit.pro  : Si IV 139 nm
   check_iris_si_fit.pro   : Si IV 140 nm
   check_iris_c_fit.pro    : C II  135 nm
   check_iris_mg_fit.pro   : Mg II 297 nm
   check_iris_o_fit.pro    : O I   134 nm
   check_iris_cl_fit.pro   : Cl I  133 nm



Please note that in case of flares or very strong emission, some of the blends or weak lines simply disappear.
at the same time, emission lines can have large Doppler shift like 100 km/s or more so might be
out of range of the normal fitting procedure. A separate analysis should be performed for such pixels.
Fe XII and Fe XXI lines are occationaly observed in 135 nm range as very broad spectral features.  



---------------------------------------
6) Comments ...
---------------------------------------
If you have comments or question, or you found a bug, please send me an e-mail at reza5pm@gmail.com
Thank you.

Reza Rezaei
