'''
Functions for spectrums and fitting procedure with lmfit.
'''
import numpy as np
from scipy.special import erfc
from scipy.interpolate import interp1d
from lmfit import Minimizer, Parameters


def funSpectrum(elements, eVs, SCF_Si, ct, sl, noise, SF0, TF0, TW0, broad_factor, lowTF0, highTF0, lowTW0, highTW0):
    '''
    Model function for the spectrum.

    Parameters
    ----------
    elements : list of Elements
        List of objects Elements.
    eVs : ndarray
        1D array of float containing the channels converted to eVs.
    SCF_Si : ndarray
        2D array of float containing the energy, f1, f2 from CXRO for Si.
    ct : float
        Constant part of the baseline = sl*eVs+ct.
    sl : float
        Linear part of the baseline = sl*eVs+ct.
    noise : float
        Width of the peak in keV, before contribution from the detector.
    SF0 : float
        Shelf fraction of the peak.
    TF0 : float
        Tail fraction of the peak.
    TW0 : float
        Tail width of the peak.
    broad_factor : float
        Broadening factor of the Compton peak's width as compared to an elastic peak.
    lowTF0 : float
        Tail fraction of the peak, on the low energy side.
    highTF0 : float
        Tail fraction of the peak, on the high energy side.
    lowTW0 : float
        Tail width of the peak, on the low energy side.
    highTW0 : float
        Tail width of the peak, on the high energy side.

    Returns
    -------
    ndarray
        spectrum_model, 1D array of float containing the spectrum fit.
    ndarray
        gaussian_part, 1D array of float containing the gaussian part of a spectrum fit.
    ndarray
        shelf_part, 1D array of float containing the shelf part of a spectrum fit.
    ndarray
        tail_part, 1D array of float containing the tail part of a spectrum fit.
    ndarray
        baseline, 1D array of float containing the baseline part of a spectrum fit.
    ndarray
        compton_part, 1D array of float containing the Compton part of a spectrum fit.
    '''
    # Function for the whole spectrum
    spectrum_model = 0.

    ##########################################################
    # Functions composing the spectrum (for control on plots)
    # Gaussian part of the peaks
    gaussian_part = 0.
    # Shelf part of the peaks
    shelf_part = 0.
    # Tail part of the peaks
    tail_part = 0.
    # Compton peaks
    compton_part = 0.


    for element in elements:

        for line in element.lines:

            # Spectrum of the given line
            spectrum_line = 0.

            # Separate the contributions to the full spectrum
            gaussian_line = 0.
            shelf_line = 0.
            tail_line = 0.
            compton_line = 0.

            # Area of the line is line.area
            # Each peak of the line has an individual area = area.line*peak.relative_strength

            for peak in line.peaks:
                if element.name == 'Compton':
                    peakArr, gaussianArr, lowTailArr, highTailArr =\
                    funComptonPeak(eVs, peak.position, peak.relative_strength, noise, broad_factor,
                                   lowTF0, highTF0, lowTW0, highTW0)
                    spectrum_line += line.area*peakArr
                    compton_line += line.area*peakArr
                    gaussian_line += line.area*gaussianArr
                    tail_line += line.area*(lowTailArr+highTailArr)

                else:
                    peakArr, gaussianArr, shelfArr, tailArr = funPeak(eVs, SCF_Si, peak.position, peak.relative_strength,
                                                                      noise, SF0, TF0, TW0)
                    spectrum_line += line.area*peakArr
                    gaussian_line += line.area*gaussianArr
                    shelf_line += line.area*shelfArr
                    tail_line += line.area*tailArr


            # Store the spectrum of the given line in an attribute
            line.peak_series = spectrum_line

            # Add the contribution of the given line to the final spectrum
            spectrum_model += spectrum_line
            gaussian_part += gaussian_line
            shelf_part += shelf_line
            tail_part += tail_line
            compton_part += compton_line


    # Set to True/False to cut the baseline at the end of the elastic peak
    cut_baseline = True

    if cut_baseline:
        # We add a linear baseline, which cannot be < 0, and stops after the elastic peak (if there is one)
        limit_baseline = eVs[-1]
        for element in elements:
            if element.name == 'Elastic':
                # There is only one elastic peak (incident beam)
                line = element.lines[0]
                peak = line.peaks[0]
                # We stop the baseline contribution at one width after the elastic peak position
                limit_baseline = peak.position+1000.*noise #noise is in keV, position in eV
        baseline_tmp = ct+sl*eVs
        baseline = np.where(eVs<=limit_baseline, baseline_tmp, 0.)
    else:
        baseline = ct+sl*eVs

    baseline = np.where(baseline>0.,baseline,0.)
    spectrum_model+= baseline

    return spectrum_model, gaussian_part, shelf_part, tail_part, baseline, compton_part


def funInterpolateSCF_Si(energy, SCF_Si):
    '''
    Interpolation of the real and imaginary parts of the atomic scattering factor.
    From CXRO, with the atomic scattering factor f = f1 + i f2.
    Used here to take into account absorption from Si within the detector.

    Parameters
    ----------
    energy : float
        energy in eV
    SCF_Si : ndarray
        2D array of float containing the energy, f1, f2 from CXRO for Si.

    Returns
    -------
    float
        f1_interp, interpolated for the given energy.
    float
        f2_interp, interpolated for the given energy.
    '''
    energy_CXRO = SCF_Si[:,0]
    f1_CXRO = SCF_Si[:,1]
    f2_CXRO = SCF_Si[:,2]

    fun_f1_interp = interp1d(energy_CXRO, f1_CXRO)
    fun_f2_interp = interp1d(energy_CXRO, f2_CXRO)

    f1_interp = float(fun_f1_interp(energy))
    f2_interp = float(fun_f2_interp(energy))

    return f1_interp, f2_interp


def funPeak(eVs, SCF_Si, pos, amp, noise, SF0, TF0, TW0):
    '''
    Definition of a peak, taking into account the contribution of the Si-based detector
    to the spectrum.
    Following:
    - M. Van Gysel, P. Lemberge & P. Van Espen, “Implementation of a spectrum fitting
    procedure using a robust peak model”, X-Ray Spectrometry 32 (2003), 434–441
    - J. Campbell & J.-X. Wang, “Improved model for the intensity of low-energy tailing in
    Si (Li) x-ray spectra”, X-Ray Spectrometry 20 (1991), 191–197

    Parameters
    ----------
    eVs : ndarray
        1D array of float containing the channels converted to eVs.
    SCF_Si : ndarray
        2D array of float containing the energy, f1, f2 from CXRO for Si.
    pos : float
        Position of the peak in eV.
    amp : float
        Amplitude of the peak.
    noise : float
        Width of the peak in keV, before contribution from the detector.
    SF0 : float
        Shelf fraction of the peak.
    TF0 : float
        Tail fraction of the peak.
    TW0 : float
        Tail width of the peak.

    Returns
    -------
    ndarray
        peakArr, 1D array containing the intensity of the model peak for each value of eVs.
    ndarray
        gaussianArr, 1D array containing the gaussian contribution to peakArr.
    ndarray
        shelfArr, 1D array containing the shelf contribution to peakArr.
    ndarray
        tailArr, 1D array containing the tail contribution to peakArr.
    '''
    # We work in keV for the peak definition
    pos_keV = pos/1000.
    keVs = eVs/1000.

    # Parameters caracteristic to the detection process
    epsilon = 0.0036
    fano = 0.115

    # Peak width after correction from detector resolution (sigmajk)
    wid = np.sqrt((noise/2.3548)**2.+epsilon*fano*pos_keV)

    # Energy dependent attenuation by Si in the detector
    atwe_Si = 28.086 #atomic weight in g/mol
    rad_el = 2.815e-15 #radius electron in m
    Na = 6.02e23 # in mol-1
    llambda = 12398./pos*1e-10 # in m
    # mass attenuation coefficient of Si in cm^2/g
    musi = 2.*llambda*rad_el*Na*1e4*float(funInterpolateSCF_Si(pos, SCF_Si)[1])/atwe_Si

    # We rescale SF to avoid too little value of the input
    SF0_r = SF0/1e4
    SF = SF0_r*musi

    # In this model we consider that TW and TF are not varying with the energy
    TW = TW0
    TF = TF0

    # Definition of gaussian
    arg = (keVs-pos_keV)**2./(2.*wid**2.)
    farg = (keVs-pos_keV)/wid
    gaussianArr = amp/(np.sqrt(2.*np.pi)*wid)*np.exp(-arg)

    # Function shelf S(i, Ejk)
    shelfArr = amp/(2.*pos_keV)*erfc(farg/np.sqrt(2.))

    # Function tail T(i, Ejk)
    # Full function:
    # tailArr = amp/(2.*wid*TW)*np.exp(farg/TW+1/(2*TW**2))*erfc(farg/np.sqrt(2.)+1./(np.sqrt(2.)*TW))
    # for numerical reasons, exp(u) very big and erfc(x) very small, we use an approximation for small TW

    if TW<1e-3:
        tailArr = amp/(np.sqrt(2.*np.pi)*wid)*np.exp(-farg**2/2.)*(1-farg*TW)
    else:
        x = farg/np.sqrt(2.)+1./(np.sqrt(2.)*TW)
        u = farg/TW+1/(2*TW**2)

        # For low u we use the full expression
        small_u = np.where(u<200, u, 0.)
        tailArr_small_u = amp/(2.*wid*TW)*np.exp(small_u)*erfc(x)

        # For high u we have to use an approximation because exp(u) overflows.
        large_u_mask = np.where(u<200, 0., 1.)
        tailArr_large_u = amp/(np.sqrt(2.*np.pi)*wid)*np.exp(-farg**2/2.)*(1-farg*TW)*large_u_mask

        tailArr = tailArr_large_u+tailArr_small_u

    # Avoid numerical instabilities
    gaussianArr = np.where(gaussianArr>1e-10,gaussianArr, 0.)
    tailArr = np.where(tailArr>1e-10,tailArr, 0.)
    shelfArr = np.where(shelfArr>1e-10,shelfArr, 0.)

    # Normalize by own fractions
    shelfArr = SF*shelfArr
    tailArr = TF*tailArr

    # Function Peak
    peakArr = np.array(gaussianArr+shelfArr+tailArr)

    return peakArr, gaussianArr, shelfArr, tailArr


def funComptonPeak(eVs, pos, amp, noise, broad_factor, lowTF0, highTF0, lowTW0, highTW0):
    '''
    Definition of the compton peak, taking into account the contribution of the Si-based detector
    to the spectrum.
    Following:
    - M. Van Gysel, P. Lemberge & P. Van Espen, “Description of Compton peaks in
    energy-dispersive x-ray fluorescence spectra”, X-Ray Spectrometry 32 (2003), 139–147.

    Parameters
    ----------
    eVs : ndarray
        1D array of float containing the channels converted to eVs.
    pos : float
        Position of the peak in eV.
    amp : float
        Amplitude of the peak.
    noise : float
        Width of the peak in keV, before contribution from the detector.
    broad_factor : float
        Broadening factor of the Compton peak's width as compared to regular peak.
    lowTF0 : float
        Tail fraction of the peak, on the low energy side.
    highTF0 : float
        Tail fraction of the peak, on the high energy side.
    lowTW0 : float
        Tail width of the peak, on the low energy side.
    highTW0 : float
        Tail width of the peak, on the high energy side.

    Returns
    -------
    ndarray
        peakArr, 1D array containing the intensity of the model peak for each value of eVs.
    ndarray
        gaussianArr, 1D array containing the gaussian contribution to peakArr.
    ndarray
        lowTailArr, 1D array containing  the tail contribution to peakArr, on the low energy side.
    ndarray
        highTailArr, 1D array containing  the tail contribution to peakArr, on the high energy side.
    '''
    # We work in keV for the peak definition
    pos_keV = pos/1000.
    keVs = eVs/1000.

    # Parameters caracteristic to the detection process
    epsilon = 0.0036
    fano = 0.115

    # Peak width after correction from detector resolution (sigmajk)
    wid = np.sqrt((noise/2.3548)**2.+epsilon*fano*pos_keV)

    # Definition of gaussian
    arg = (keVs-pos_keV)**2./(2.*(broad_factor*wid)**2.)
    farg = (keVs-pos_keV)/wid
    gaussianArr = amp/(np.sqrt(2.*np.pi)*broad_factor*wid)*np.exp(-arg)

    #Low energy tail TA
    if lowTW0<1e-3:
        lowTailArr = amp/(np.sqrt(2.*np.pi)*wid)*np.exp(-farg**2/2.)*(1-farg*lowTW0)
    else:
        x = farg/np.sqrt(2.)+1./(np.sqrt(2.)*lowTW0)
        u = farg/lowTW0+1/(2*lowTW0**2)

        # For low u we use the full expression
        small_u = np.where(u<200, u, 0.)
        lowTailArr_small_u = amp/(2.*wid*lowTW0)*np.exp(u)*erfc(x)

        # For high u we have to use an approximation because exp(u) overflows.
        large_u_mask = np.where(u<200, 0., 1.)
        lowTailArr_large_u = amp/(np.sqrt(2.*np.pi)*wid)*np.exp(-farg**2/2.)*(1-farg*lowTW0)*large_u_mask

        lowTailArr = lowTailArr_large_u+lowTailArr_small_u


    #High energy tail TB
    if highTW0<1e-3:
        highTailArr = amp/(np.sqrt(2.*np.pi)*wid)*np.exp(-farg**2/2.)*(1+farg*highTW0)
    else:
        x = -farg/np.sqrt(2.)+1./(np.sqrt(2.)*highTW0)
        u = -farg/highTW0+1/(2*highTW0**2)

        # For low u we use the full expression
        small_u = np.where(u<200, u, 0.)
        highTailArr_small_u = amp/(2.*wid*highTW0)*np.exp(u)*erfc(x)

        # For high u we have to use an approximation because exp(u) overflows.
        large_u_mask = np.where(u<200, 0., 1.)
        highTailArr_large_u = amp/(np.sqrt(2.*np.pi)*wid)*np.exp(-farg**2/2.)*(1+farg*highTW0)*large_u_mask

        highTailArr = highTailArr_large_u+highTailArr_small_u


    # Normalize by own fractions
    lowTailArr = lowTF0*lowTailArr
    highTailArr = highTF0*highTailArr

    # Avoid numerical instabilities
    lowTailArr = np.where(lowTailArr>1e-10,lowTailArr, 0.)
    highTailArr = np.where(highTailArr>1e-10,highTailArr, 0.)

    # Function Peak
    peakArr = np.array(gaussianArr+lowTailArr+highTailArr)

    # Avoid numerical instabilities
    peakArr = np.where(peakArr>1e-10,peakArr, 0.)


    return peakArr, gaussianArr, lowTailArr, highTailArr


def funFitSpectrum(spectrum_to_fit, list_isfit, elements, channels, SCF_Si, gain, eV0, ct, sl, noise,
                   SF0, TF0, TW0, broad_factor, lowTF0, highTF0, lowTW0, highTW0, 
                   is_bckg_subset, bckg_eVs_inf, bckg_eVs_sup):
    '''
    Fit procedure using lmfit.

    Parameters
    ----------
    spectrum_to_fit : ndarray
        1D array of float containing the spectrum to fit.
    list_isfit : list of str
        List of params names which are active fit parameters.
    elements : list of Elements
        List of objects Elements.
    channels : ndarray
        1D array of int corresponding to the subset of selected channels.
    SCF_Si : ndarray
        2D array of float containing the energy, f1, f2 from CXRO for Si.
    gain : float
        Parameter in the conversion of channels to eVs, eVs = channels*gain + eV0.
    eV0 : float
        Parameter in the conversion of channels to eVs, eVs = channels*gain + eV0.
    ct : float
        Constant part of the baseline = sl*eVs+ct.
    sl : float
        Linear part of the baseline = sl*eVs+ct.
    noise : float
        Width of the peak in keV, before contribution from the detector.
    SF0 : float
        Shelf fraction of the peak.
    TF0 : float
        Tail fraction of the peak.
    TW0 : float
        Tail width of the peak.
    broad_factor : float
        Broadening factor of the Compton peak's width as compared to an elastic peak.
    lowTF0 : float
        Tail fraction of the peak, on the low energy side.
    highTF0 : float
        Tail fraction of the peak, on the high energy side.
    lowTW0 : float
        Tail width of the peak, on the low energy side.
    highTW0 : float
        Tail width of the peak, on the high energy side.
    is_bckg_subset : bool
        True if fitting the background on a subset of the spectrum.
    bckg_eVs_inf : float
        Lower bound of the fitting range for the background (if is_bckg_subset).
    bckg_eVs_sup : float
        Upper bound of the fitting range for the background (if is_bckg_subset).

    Returns
    -------
    object lmfit.MinimizerResult
        result, result of the fit.
    '''
    ############################################################
    # Least-square fit
    # We fit first the line areas only with a least-square fit
    # To get initial conditions for the LM fit

    # Set the area of each line to 1.
    for element in elements:
        for line in element.lines:
            line.area = 1.

    eVs = gain*channels + eV0

    # Compute the model spectrum
    # The aim is to generate line.peak_series for each line
    # line.peak_series is the array containing the part of the spectrum corresponding to the given line
    _, _, _, _, _, _ = funSpectrum(elements, eVs, SCF_Si, ct, sl, noise,
                                   SF0, TF0, TW0, broad_factor,
                                   lowTF0, highTF0, lowTW0, highTW0)

    # Do the LS fit, returning the best area of each line (not considering the baseline)
    a = []
    for element in elements:
        for line in element.lines:
            a.append(line.peak_series)
    area_LS, _, _, _ = np.linalg.lstsq(np.transpose(a), spectrum_to_fit, 1.e-10)

    # Set the area of each line to the value obtained by the LS fit
    # They will be used as initial guess in the LM fit
    i=0
    for element in elements:
        for line in element.lines:
            line.area = area_LS[i]
            i+=1

    ############################################################
    # LM fit

    def funToBeMinimized(dparams, SCF_Si, spectrum_to_fit):
        '''
        Define the objective function to be minized by the LM fit.

        Parameters
        ----------
        dparams : dict
            The fit parameters.
        SCF_Si : ndarray
            2D array of float containing the energy, f1, f2 from CXRO for Si.
        spectrum_to_fit : ndarray
            The experimental spectrum.

        Returns
        -------
        ndarray
            residuals, the 1D array to be minimized.
        '''
        # Construct the array eVs with the current fit values
        # of gain and eV0
        eVs = dparams['gain']*channels + dparams['eV0']

        # Set the area of each line
        for element in elements:
            for line in element.lines:
                line.area = dparams['area_'+line.name]
                for peak in line.peaks:
                    if peak.is_fitpos:
                        peak.position = dparams['pos_'+peak.name]


        # Compute the model
        spectrum_model, _, _, _, _, _ =\
                        funSpectrum(elements, eVs, SCF_Si, dparams['ct'], dparams['sl'], dparams['noise'],
                                    dparams['SF0'], dparams['TF0'], dparams['TW0'], dparams['broad_factor'],
                                    dparams['lowTF0'], dparams['highTF0'], dparams['lowTW0'], dparams['highTW0'])

        # Compute the residuals
        residuals = spectrum_model - spectrum_to_fit

        return residuals


    ####################################################
    # Create the dictonnary of parameters for the LM fit

    dparams = Parameters()

    for element in elements:
        for line in element.lines:

            # The area of each line is fitted, and cannot be negative
            dparams.add('area_'+line.name, value=line.area, vary=True, min=0.)

            # Add the positions of the peaks which should be fitted
            # We do not allow variation of the peak position of more than 100 eVs
            for peak in line.peaks:
                if peak.is_fitpos:
                    dparams.add('pos_'+peak.name, value=peak.position,
                                vary=True, min=peak.position-100, max=peak.position+100)

    # is_bckg_subset True if the background is fitted only on a subset of the spectrums
    if (is_bckg_subset and 'ct' in list_isfit):  

        if 'sl' in list_isfit:
            # Fit of the background on the subset with a*x+b
            ind_inf = np.argmin(np.abs(eVs-bckg_eVs_inf))
            ind_sup = np.argmin(np.abs(eVs-bckg_eVs_sup))
            p, cov = np.polyfit(eVs[ind_inf:ind_sup], spectrum_to_fit[ind_inf:ind_sup], 1, cov=True)
            dparams.add('ct', value=p[1], vary=False)
            dparams.add('sl', value=p[0], vary=False)

        else:
            # Fit of the background on the subset with a constant
            ind_inf = np.argmin(np.abs(eVs-bckg_eVs_inf))
            ind_sup = np.argmin(np.abs(eVs-bckg_eVs_sup))
            p, cov = np.polyfit(channels[ind_inf:ind_sup], spectrum_to_fit[ind_inf:ind_sup], 0, cov=True)
            dparams.add('ct', value=p[0], vary=False)
            dparams.add('sl', value=sl, vary=False)
 
    else:
        # Fit on the whole spectrum with classical parameters of lmfit
        dparams.add('ct', value=ct, vary='ct' in list_isfit, min=0.)
        dparams.add('sl', value=sl, vary='sl' in list_isfit)  
                    
    dparams.add('gain', value=gain, vary='gain' in list_isfit, min=0.)
    dparams.add('eV0', value=eV0, vary='eV0' in list_isfit)
    dparams.add('noise', value=noise, vary='noise' in list_isfit, min=0.)
    dparams.add('SF0', value=SF0, vary='SF0' in list_isfit, min=0.)
    dparams.add('TF0', value=TF0, vary='TF0' in list_isfit, min=0.)
    dparams.add('TW0', value=TW0, vary='TW0' in list_isfit, min=0.)
    dparams.add('broad_factor', value=broad_factor, vary='broad_factor' in list_isfit, min=0.)
    dparams.add('lowTF0', value=lowTF0, vary='lowTF0' in list_isfit, min=0.)
    dparams.add('highTF0', value=highTF0, vary='highTF0' in list_isfit, min=0.)
    dparams.add('lowTW0', value=lowTW0, vary='lowTW0' in list_isfit, min=0.)
    dparams.add('highTW0', value=highTW0, vary='highTW0' in list_isfit, min=0.)

    # Do the fit, here with leastsq model
    minner = Minimizer(funToBeMinimized, dparams, fcn_args=([SCF_Si, spectrum_to_fit]),
                       xtol=1e-6, ftol=1e-6, max_nfev=1000)

    result = minner.minimize(method='leastsq')
    
    # Redefine the result to account for the fitting of the background
    if (is_bckg_subset and 'ct' in list_isfit):   
        if 'sl' in list_isfit:
            result.params['ct'].vary = True    
            result.params['sl'].vary = True
            if result.errorbars:
                result.params['ct'].stderr = np.sqrt(np.diag(cov))[1]
                result.params['sl'].stderr = np.sqrt(np.diag(cov))[0]
        else:
            result.params['ct'].vary = True    
            result.params['sl'].vary = False
            if result.errorbars:
                result.params['ct'].stderr = np.sqrt(np.diag(cov))[0]           
    
    return result
