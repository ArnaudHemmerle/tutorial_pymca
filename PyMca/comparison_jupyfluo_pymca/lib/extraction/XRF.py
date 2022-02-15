'''
Module for extracting, plotting, and saving an XRF scan.
'''
import os
import sys
import io
from contextlib import redirect_stdout
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from lib.extraction import PyNexus as PN

# Define colors for prints
_RED='\x1b[31;01m'
_RESET="\x1b[0m"
_BLUE='\x1b[34;01m'

def Treat(nxs_filename, recording_dir, list_elems,
          absorbers='', logz=True, first_channel=0, last_channel=2048,
          use_eV=False, gain=10., eV0=0., arr_peaks=None, working_dir='', fast=True,
          show_data_stamps=False, plot_spectrogram=False, plot_sum=False, plot_first_last=False,
          save=False, verbose=False):
    '''
    Call functions for extracting, plotting, and saving an XRF scan.

    Parameters
    ----------
    nxs_filename : str
        Nexus filename.
    recording_dir : str
        Directory where the nexus file is stored.
    list_elems : ndarray
        An array with the elements to extract, for ex. list_elems = [1, 2, 3].
    absorbers : str, optional
        Text to display indicating which absorber was used.
    logz : bool, optional
        Log scale on the plots.
    first_channel : int, optional
        The spectrums will be extracted between first_channel and last_channel.
    last_channel : int, optional
        The spectrums will be extracted between first_channel and last_channel.
    use_eV : bool, optional
        Convert the channels to eVs.
    gain : float, optional
        Channels are converted to eVs following eVs = gain*channel+eV0.
    ev0 : float, optional
        Channels are converted to eVs following eVs = gain*channel+eV0.
    arr_peaks : ndarray, optional
        An array with the peaks to display, for ex. arr_peaks = [('Elastic', '12000.'), ('Compton', '11670.')].
    working_dir : str, optional
        Directory where the treated files will be stored.
    fast : bool, optional
        Trigger fast extract of the nexus file.
    show_data_stamps : bool, optional
        Print the list of sensors from the nexus file.
    plot_spectrogram : bool, optional
        Plot the spectrogram.
    plot_sum : bool, optional
        Plot the sum of the spectrums over time.
    plot_first_last : bool, optional
        Plot the first and last spectrum.
    save : bool, optional
        Save the GIXD.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    ndarray
        channels, an array containing the channels.
    ndarray
        eVs, an array containing the channels converted to eVs.
    ndarray
        spectrums, an array containing the spectrums.

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file is not found.
    SystemExit('ICR not found')
        when no ICR is found (most likely because the wrong elements were given).
    '''

    # Deal with empty arr_peaks
    if arr_peaks is None:
        arr_peaks = [(None,None)]

    channels, eVs, spectrums, first_non_zero_spectrum, last_non_zero_spectrum = \
            Extract(nxs_filename, recording_dir,
                    list_elems, first_channel, last_channel,
                    gain, eV0, fast, show_data_stamps, verbose)

    if plot_spectrogram or plot_first_last or plot_sum:
        Plot(channels, eVs, spectrums, first_non_zero_spectrum, last_non_zero_spectrum,
             use_eV, arr_peaks, absorbers, logz,
             nxs_filename, plot_spectrogram, plot_sum, plot_first_last)

    if save:
        Save(nxs_filename, recording_dir, fast, working_dir, verbose)

    return channels, eVs, spectrums

def Extract(nxs_filename, recording_dir,
            list_elems, first_channel, last_channel,
            gain, eV0, fast, show_data_stamps, verbose):
    '''
    Extract the nexus scan and return useful quantities for XRF.

    Parameters
    ----------
    nxs_filename : str
        Nexus filename.
    recording_dir : str
        Directory where the nexus file is stored.
    list_elems : ndarray
        An array with the elements to extract, for ex. list_elems = [1, 2, 3].
    first_channel : int, optional
        The spectrums will be extracted between first_channel and last_channel.
    last_channel : int
        The spectrums will be extracted between first_channel and last_channel.
    gain : float
        Channels are converted to eVs following eVs = gain*channel+eV0.
    ev0 : float
        Channels are converted to eVs following eVs = gain*channel+eV0.
    fast : bool, optional
        Trigger fast extract of the nexus file.
    show_data_stamps : bool, optional
        Print the list of sensors from the nexus file.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    ndarray
        channels, an array containing the channels.
    ndarray
        eVs, an array containing the channels converted to eVs.
    ndarray
        spectrums, an array containing the spectrums.
    int
        first_non_zero_spectrum, index of the first scan extracted.
    int
        last_non_zero_spectrum, index of the last scan extracted.

    Raises
    ------
    SystemExit('Nexus not found')
        when Nexus file is not found.
    SystemExit('ICR not found')
        when no ICR is found (most likely because the wrong elements were given).
    '''
    nxs_path = recording_dir+nxs_filename

    if not os.path.isfile(nxs_path):
        print(_RED+'Scan %s seems not to exist in recording directory'%(nxs_filename)+_RESET)
        print(('\t\t recording directory : '+recording_dir))
        sys.exit('Nexus not found')

    else:

        if verbose:
            print(_BLUE+" - Open Nexus Data File :"+ _RESET)
        if verbose:
            print('\t'+nxs_path)
        try:
            nexus=PN.PyNexusFile(nxs_path, fast=fast)
        except OSError:
            print(_RED,'\t Nexus file seems not to exist or is not correct',_RESET)
            sys.exit('Nexus not found')

        nbpts=np.int(nexus.get_nbpts())
        if verbose:
            print("\t. Number of data points: ", nbpts)

        # Get stamps
        stamps, data= nexus.extractData()

        nexus.close()

        if show_data_stamps :
            print("\t. Available Counters:")
        for i, stamp in enumerate(stamps):
            if stamp[1] is not None:
                if show_data_stamps:
                    print("\t\t", i, ' -------> ', stamp[1])
            else:
                if show_data_stamps:
                    print("\t\t",i, ' -------> ', stamp[0])


    def extract_and_correct(ind_spectrum):
        '''
        Extract the requested fluospectrum from the nexus file and correct it with ICR/OCR
        '''
        is_icr_found = False
        is_ocr_found = False
        for i, stamp in enumerate(stamps):
            if (stamp[1] is not None and stamp[1].lower() == "fluoicr0"+ind_spectrum):
                fluoicr = data[i]
                is_icr_found = True
            if (stamp[1] is not None and stamp[1].lower() == "fluoocr0"+ind_spectrum):
                fluoocr = data[i]
                is_ocr_found = True
            if (stamp[1] is not None and stamp[1].lower() == "fluospectrum0"+ind_spectrum):
                fluospectrum = data[i]
            if (stamp[1] is None and stamp[0].lower() == "integration_time"):
                integration_time = data[i]

        if is_icr_found:
            ICR = fluoicr
            if is_ocr_found:
                OCR = fluoocr
            else:
                print(_RED+"OCR not found in data. Taking OCR = spectrum_intensity/counting_time."+_RESET)
                OCR = np.array([np.sum(fluospectrum[n])/integration_time[n] for n in range(len(fluospectrum))])

            ratio = np.array([ICR[n]/OCR[n] if (~np.isclose(OCR[n],0.) & ~np.isnan(OCR[n]) & ~np.isnan(ICR[n]))
                              else 0. for n in range(len(ICR))])
            spectrums_corr = np.array([fluospectrum[n]*ratio[n] for n in range(len(ratio))])
            return spectrums_corr

        print(_RED+"ICR not found in data. Check if the box \'Elements\' is right."+_RESET)
        print(_RED+"Try to put 4 in the box \'Elements\' for the single-element detector."+_RESET)
        print(_RED+"Try to put 0, 1, 2, 3 in the box \'Elements\' for the four-elements detector."+_RESET)
        sys.exit('ICR not found.')

    # Correct each chosen element with ICR/OCR and sum them
    allspectrums_corr = np.zeros((nbpts, 2048))

    for i in list_elems:
        allspectrums_corr += extract_and_correct(str(i))

    ind_non_zero_spectrums = np.where(np.sum(allspectrums_corr, axis = 1)>10.)[0]
    first_non_zero_spectrum = ind_non_zero_spectrums[0]
    last_non_zero_spectrum = ind_non_zero_spectrums[-1]

    channels = np.arange(int(first_channel), int(last_channel+1))
    eVs = channels*gain+eV0
    spectrums = allspectrums_corr[0:last_non_zero_spectrum+1,
                                  int(first_channel):int(last_channel+1)]

    return channels, eVs, spectrums, first_non_zero_spectrum, last_non_zero_spectrum


def Extract_mat(filename, recording_dir,
                list_elems, first_channel, last_channel,
                gain, eV0, show_data_stamps, verbose):
    '''
    Extract the spectrums from .mat files and return useful quantities extracted from the .dat file.

    Parameters
    ----------
    filename : str
        Filename without the extension.
    recording_dir : str
        Directory where the mat and dat files are stored.
    list_elems : ndarray
        An array with the elements to extract, for ex. list_elems = [1, 2, 3].
    first_channel : int, optional
        The spectrums will be extracted between first_channel and last_channel.
    last_channel : int
        The spectrums will be extracted between first_channel and last_channel.
    gain : float
        Channels are converted to eVs following eVs = gain*channel+eV0.
    ev0 : float
        Channels are converted to eVs following eVs = gain*channel+eV0.
    show_data_stamps : bool, optional
        Print the list of sensors from the dat file.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    ndarray
        channels, an array containing the channels.
    ndarray
        eVs, an array containing the channels converted to eVs.
    ndarray
        spectrums, an array containing the spectrums.
    int
        first_non_zero_spectrum, index of the first scan extracted.
    int
        last_non_zero_spectrum, index of the last scan extracted.

    Raises
    ------
    SystemExit('.mat not found')
        when mat file is not found.
    SystemExit('.dat not found')
        when dat file is not found.
    SystemExit('ICR not found')
        when no ICR is found (most likely because the wrong elements were given).
    '''
    file_path = recording_dir+filename
    path_to_dat = file_path+'.dat'

    if not os.path.isfile(path_to_dat):
        print(_RED+'Scan %s seems not to exist in recording directory'%(filename+'.dat')+_RESET)
        print(('\t\t recording directory : '+recording_dir))
        sys.exit('.dat not found')

    else:
        
        # Extract list of detector elements available
        SDD_elems_available = []
        nbpts = 0
        for index_element in list_elems:
            path_to_mat = file_path+'_fluospectrum0'+str(index_element)+'.mat'

            if os.path.isfile(path_to_mat):
                SDD_elems_available.append(str(index_element))

                # Extract number of spectrums taken during the scan
                nbpts = np.shape(np.genfromtxt(path_to_mat))[0]
    
        
        if SDD_elems_available == []:
            print(_RED,'\t No .mat file found.',_RESET)
            sys.exit('.mat not found')

        if verbose:
            print("\t. Number of data points: ", nbpts)

        # Get stamps and data
        dat_extracted = np.genfromtxt(path_to_dat, names=True)
        stamps = [(name, name, None) for name in dat_extracted.dtype.names]
        data = np.genfromtxt(path_to_dat, skip_header=1).transpose()
    
        if show_data_stamps :
            print("\t. Available Counters:")
        for i, stamp in enumerate(stamps):
            if stamp[1] is not None:
                if show_data_stamps:
                    print("\t\t", i, ' -------> ', stamp[1])
            else:
                if show_data_stamps:
                    print("\t\t",i, ' -------> ', stamp[0])


    def extract_and_correct(ind_spectrum):
        '''
        Extract the requested fluospectrum from the mat file and correct it with ICR/OCR from the dat file
        '''
        is_icr_found = False
        is_ocr_found = False

        for i, stamp in enumerate(stamps):
            if (stamp[1] is not None and stamp[1].lower() == "fluoicr0"+ind_spectrum):
                fluoicr = data[i]
                is_icr_found = True
            if (stamp[1] is not None and stamp[1].lower() == "fluoocr0"+ind_spectrum):
                fluoocr = data[i]
                is_ocr_found = True
            if (stamp[1] is None and stamp[0].lower() == "integration_time"):
                integration_time = data[i]

            fluospectrum = np.genfromtxt(file_path+'_fluospectrum0'+str(ind_spectrum)+'.mat')
                
        if is_icr_found:
            ICR = fluoicr
            if is_ocr_found:
                OCR = fluoocr
            else:
                print(_RED+"OCR not found in data. Taking OCR = spectrum_intensity/counting_time."+_RESET)
                OCR = np.array([np.sum(fluospectrum[n])/integration_time[n] for n in range(len(fluospectrum))])

            ratio = np.array([ICR[n]/OCR[n] if (~np.isclose(OCR[n],0.) & ~np.isnan(OCR[n]) & ~np.isnan(ICR[n]))
                              else 0. for n in range(len(ICR))])
            spectrums_corr = np.array([fluospectrum[n]*ratio[n] for n in range(len(ratio))])
            return spectrums_corr

        print(_RED+"ICR not found in data. Check if the box \'Elements\' is right."+_RESET)
        print(_RED+"Try to put 4 in the box \'Elements\' for the single-element detector."+_RESET)
        print(_RED+"Try to put 0, 1, 2, 3 in the box \'Elements\' for the four-elements detector."+_RESET)
        sys.exit('ICR not found.')

    # Correct each chosen element with ICR/OCR and sum them
    allspectrums_corr = np.zeros((nbpts, 2048))

    for i in list_elems:
        allspectrums_corr += extract_and_correct(str(i))

    ind_non_zero_spectrums = np.where(np.sum(allspectrums_corr, axis = 1)>10.)[0]
    first_non_zero_spectrum = ind_non_zero_spectrums[0]
    last_non_zero_spectrum = ind_non_zero_spectrums[-1]

    channels = np.arange(int(first_channel), int(last_channel+1))
    eVs = channels*gain+eV0
    spectrums = allspectrums_corr[0:last_non_zero_spectrum+1,
                                  int(first_channel):int(last_channel+1)]

    return channels, eVs, spectrums, first_non_zero_spectrum, last_non_zero_spectrum




def Plot(channels, eVs, spectrums, first_non_zero_spectrum, last_non_zero_spectrum,
         use_eV, arr_peaks, absorbers, logz,
         nxs_filename, plot_spectrogram, plot_sum, plot_first_last):
    '''
    Plot XRF data.

    Parameters
    ----------
    channels : ndarray
        The channels.
    eVs : ndarray
        The channels converted to eVs.
    spectrums : ndarray
        The spectrums.
    first_non_zero_spectrum : int
        Index of the first scan extracted.
    last_non_zero_spectrum : int
        Index of the last scan extracted.
    use_eV : bool
        Convert the channels to eVs.
    arr_peaks : ndarray, optional
        An array with the peaks to display, for ex. arr_peaks = [('Elastic', '12000.'), ('Compton', '11670.')].
    absorbers : str, optional
        Text to display indicating which absorber was used.
    logz : bool
        Log scale on the plots.
    nxs_filename : str
        Nexus filename.
    plot_spectrogram : bool, optional
        Plot the spectrogram.
    plot_sum : bool, optional
        Plot the sum of the spectrums over time.
    plot_first_last : bool, optional
        Plot the first and last spectrum.
    '''
    # Print absorbers
    if absorbers != '':
        print("\t. Absorbers:", str(absorbers))

    if plot_spectrogram:

        fig = plt.figure(figsize=(12,4.6))
        ax1 = fig.add_subplot(111)
        ax1.set_title(nxs_filename.split('\\')[-1], fontsize='x-large')
        ax1.set_xlabel('spectrum index', fontsize='large')
        ax1.set_xlim(left = 0, right = last_non_zero_spectrum)

        if use_eV:
            xx, yy = np.meshgrid(np.arange(0,last_non_zero_spectrum+1), eVs)
            ax1.set_ylabel('eV', fontsize='large')
        else:
            xx, yy = np.meshgrid(np.arange(0,last_non_zero_spectrum+1), channels)
            ax1.set_ylabel('channel', fontsize='large')

        if logz:
            ax1.pcolormesh(xx, yy, spectrums.transpose(), cmap='viridis', shading = 'auto',
                           norm = colors.LogNorm(), rasterized=True)
        else:
            ax1.pcolormesh(xx, yy, spectrums.transpose(), cmap='viridis', shading = 'auto',
                           rasterized=True)

        plt.show()

    if plot_sum:
        fig = plt.figure(figsize=(12,4.5))
        ax1 = fig.add_subplot(111)
        ax1.set_ylabel('counts', fontsize='large')
        if logz:
            ax1.set_yscale('log')
        if use_eV:
            ax1.set_xlabel('eV', fontsize='large')
            line1, = ax1.plot(eVs, np.sum(spectrums, axis = 0), 'b.-', label='Sum of spectrums')
        else:
            ax1.set_xlabel('channel', fontsize='large')
            line1, = ax1.plot(channels, np.sum(spectrums, axis = 0), 'b.-', label='Sum of spectrums')

        if arr_peaks[0][0] is not None :

            # Plot the peak positions

            # Prepare a list of colors and linestyles
            colors_axv = iter(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
                               '#7f7f7f', '#bcbd22', '#17becf']*20)
            linestyles_axv = iter(['--', '-.', '-', ':']*40)

            # Rearrange the peaks to plot them by increasing energy
            arr_peaks = np.array(arr_peaks)
            arg_position_peaks = np.argsort([float(elem[1]) for elem in arr_peaks])
            val_position_peaks = arr_peaks[arg_position_peaks][:,1]
            labels_peaks = arr_peaks[arg_position_peaks][:,0]

            axvlines = []
            for i in range(len(arr_peaks)):
                axvlines.append(ax1.axvline(float(val_position_peaks[i]), label = str(labels_peaks[i]),
                                            color = next(colors_axv), linestyle = next(linestyles_axv)))

            axvlegends = ax1.legend(handles=axvlines, fontsize=10, ncol = len(arr_peaks)//16+1,
                                    bbox_to_anchor=(1.01, 1.), loc='upper left',  borderaxespad=0.)
            plt.gca().add_artist(axvlegends)

        ax1.legend(handles=[line1], fontsize='large', loc='upper left')
        plt.show()

    if plot_first_last:
        #Plot the selected channel range
        fig = plt.figure(figsize=(12,4.5))
        ax1 = fig.add_subplot(111)
        ax1.set_ylabel('counts', fontsize='large')
        if logz:
            ax1.set_yscale('log')
        if use_eV:
            ax1.set_xlabel('eV', fontsize='large')
            line1, = ax1.plot(eVs, spectrums[first_non_zero_spectrum], 'b.-', label='First spectrum')
            line2, = ax1.plot(eVs, spectrums[-1], 'r.-', label='Last spectrum')
        else:
            ax1.set_xlabel('channel', fontsize='large')
            line1, = ax1.plot(channels, spectrums[first_non_zero_spectrum], 'b.-', label='First spectrum')
            line2, = ax1.plot(channels, spectrums[-1], 'r.-', label='Last spectrum')

        if arr_peaks[0][0] is not None :

            # Rearrange the peaks to plot them by increasing energy
            arr_peaks = np.array(arr_peaks)
            arg_position_peaks = np.argsort([float(elem[1]) for elem in arr_peaks])
            val_position_peaks = arr_peaks[arg_position_peaks][:,1]
            labels_peaks = arr_peaks[arg_position_peaks][:,0]

            axvlines = []
            for i in range(len(arr_peaks)):
                axvlines.append(ax1.axvline(float(val_position_peaks[i]), label = str(labels_peaks[i]),
                                            color = next(colors_axv), linestyle = next(linestyles_axv)))

            axvlegends = ax1.legend(handles=axvlines, fontsize=10, ncol = len(arr_peaks)//16+1,
                                    bbox_to_anchor=(1.01, 1.), loc='upper left',  borderaxespad=0.)
            plt.gca().add_artist(axvlegends)

        ax1.legend(handles=[line1, line2], fontsize='large', loc='upper left')
        plt.show()


def Save(nxs_filename, recording_dir, fast, working_dir, verbose):
    '''
    Use the PyNexus library to convert the Nexus file into a .dat file.

    Parameters
    ----------
    nxs_filename : str
        Nexus filename.
    recording_dir : str
        Directory where the nexus file is stored.
    fast : bool
        Trigger fast extract of the nexus file.
    working_dir : str
        Directory where the treated files will be stored.
    verbose : bool
        Verbose mode.
    '''
    savename = working_dir + nxs_filename[:nxs_filename.rfind('.nxs')]

    # We assume extraction was already checked with Extract
    nxs_path = recording_dir+nxs_filename
    nexus=PN.PyNexusFile(nxs_path, fast=fast)

    # Get stamps and Data
    stamps, data = nexus.extractData()

    f = io.StringIO()
    # Avoid printing sensors in the notebook
    with redirect_stdout(f):
        old_nexus_filename = nexus.filename
        # Save in working dir
        nexus.filename = working_dir+nxs_filename
        nexus.savePointExtractedData((stamps, data))
        nexus.saveOneDExtractedData((stamps, data))
        nexus.filename = old_nexus_filename
    out = f.getvalue()

    nexus.close()

    if verbose:
        print('\t. 0D data saved in:')
        print("\t", savename+'.dat')
        print('\t. Spectrum(s) saved in:')
        print("\t", savename+'_fluospectrum*.mat')
