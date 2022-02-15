'''
Class regrouping all methods and attributes for an experiment.
'''
import os
import sys
import time
import subprocess
import csv
from copy import deepcopy
from datetime import datetime
from ast import literal_eval
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from lmfit import fit_report
from IPython.display import clear_output
import ipysheet
from lib.frontend import notebook as nb
from lib.extraction import PyNexus as PN
from lib.fit import funFit

# Define colors for prints
_RED='\x1b[31;01m'
_BOLD="\x1b[01;06m"
_RESET="\x1b[0m"

class Experiment:
    '''
    Class for an experiment.

    Attributes
    ----------
    all_spectrums : ndarray
        2D array of float containing all the spectrums extracted.
    baseline : ndarray
        1D array of float containing the baseline part of a spectrum fit.
    bckg_eVs_inf : float
        Lower bound of the fitting range for the background (if is_bckg_subset).
    bckg_eVs_sup : float
        Upper bound of the fitting range for the background (if is_bckg_subset).        
    beam_energy : float
        Beam energy in eV.
    broad_factor : float
        Broadening factor of the Compton peak's width as compared to an elastic peak.
    compton_part : ndarray
        1D array of float containing the Compton part of a spectrum fit.
    ct : float
        Constant part of the baseline = sl*eVs+ct.
    channels : ndarray
        1D array of int corresponding to the subset of selected channels.
    current_sensorsRelTimestamps : float
        Time stamp of the spectrum fitted.
    data0D : ndarray
        2D array of float containing the values of the sensors recorded during the scan.
    elements : list of Elements
        List of objects Elements.
    elements_fit : list of Elements
        List of objects Elements, used in the fitting process.
    eV0 : float
        Parameter in the conversion of channels to eVs, eVs = channels*gain + eV0.
    eVs : ndarray
        1D array of float containing the channels converted to eVs.
    eVs_fit : ndarray
        1D array of float containing the channels converted to eVs, used in the fitting process.
    files_dir : str
        Directory where the nexus files are.
    first_channel : int
        Index of the first channel after the extraction.
    first_spectrum : int
        Index of the first spectrum after the extraction.
    gain : float
        Parameter in the conversion of channels to eVs, eVs = channels*gain + eV0.
    gaussian_part : ndarray
        1D array of float containing the gaussian part of a spectrum fit.
    highTF0 : float
        Tail fraction of the peak, on the high energy side.
    highTW0 : float
        Tail width of the peak, on the high energy side.
    is_bckg_subset : bool
        True if fitting the background on a subset of the spectrum.
    is_broad_factor : bool
        True if broad_factor is a fit parameter.
    is_ct : bool
        True if ct is a fit parameter.
    is_eV0 : bool
        True if eV0 is a fit parameter.
    is_extract_done : bool
        True if an extraction has been done.
    is_fit_done : book
        True if a fit has been done.
    is_fit_fail : bool
        True if a fit did not converge.
    is_gain : bool
        True if gain is a fit parameter.
    is_highTF0 : bool
        True if highTF0 is a fit parameter.
    is_highTW0 : bool
        True if highTW0 is a fit parameter.
    is_last_fit_a_preview : bool
        True if the last fit done was a preview (single fit).
    is_lowTF0 : bool
        True if lowTF0 is a fit parameter.
    is_lowTW0 : bool
        True if lowTW0 is a fit parameter.
    is_noise : bool
        True if noise is a fit parameter.
    is_SDD_elem_0 : bool
        True if the SDD element 0 is used.
    is_SDD_elem_1 : bool
        True if the SDD element 1 is used.
    is_SDD_elem_2 : bool
        True if the SDD element 2 is used.
    is_SDD_elem_3 : bool
        True if the SDD element 3 is used.
    is_SDD_elem_4 : bool
        True if the SDD element 4 is used.
    is_SF0 : bool
        True if SF0 is a fit parameter.
    is_sheet_loaded : bool
        True if the ipysheet has been loaded.
    is_sl : bool
        True if sl is a fit parameter.
    is_spectrum_empty : bool
        True if the spectrum is considered as empty.
    is_TF0 : bool
        True if TF0 is a fit parameter.
    is_TW0 : bool
        True if TW0 is a fit parameter.
    last_channel : int
        Index of the last channel after the extraction.
    last_non_zero_spectrum : int
        Index of the last non-zeron spectrum before extraction (all spectrums after are empty).
    last_spectrum : int
        Index of the last spectrum after the extraction.
    list_extract_params_files : list of str
        List of csv files with extraction params.
    list_fit_params_files : list of str
        List of csv files with fit params.
    list_isfit : list of str
        List of params names which are active fit parameters.
    list_isfit_str : str
        list_isfit converted to a string, with a comma as a separator.
    list_files_filenames : list of str
        List of files in the 'files' directory.
    list_peaks_params_files : list of str
        List of csv files with peaks params.
    logs : str
        String to be displayed in the logs window.
    logs_lvl : int
        Depth level of the logs.
    lowTF0 : float
        Tail fraction of the peak, on the low energy side.
    lowTW0 : float
        Tail width of the peak, on the low energy side.
    min_strength : float
        Minimum strength of peak to be displayed in the peak selection widget.
    name : str
        Name of the scan.
    nb_allspectrums : int
        Number of spectrums taken during the scan.
    nb_spectrums : int
        Number of spectrums extracted.
    noise : float
        Width of the peak in keV, before contribution from the detector.
    notebook_name : str
        Name of the notebook (necessary for saving in pdf).
    filename : str
        Name of the nexus file.
    path_to_db : str
        Relative path to the xraylib database.
    path_to_extract_params : str
        Relative path to the csv file with loaded extraction parameters.
    path_to_extract_params_default : str
        Relative path to the csv file with default extraction parameters.
    path_to_extract_params_save : str
        Relative path to the csv file with the saved extraction parameters.
    path_to_fit_curve_results : str
        Relative path to the csv file with the fit curve of a given spectrum.
    path_to_fit_folder : str
        Relative path to the folder containing the fit curves.
    path_to_fit_log_results : str
        Relative path to the log file with logs from the fit.
    path_to_fit_params : str
        Relative path to the csv file with loaded fit parameters.
    path_to_fit_params_default : str
        Relative path to the csv file with default fit parameters.
    path_to_fit_params_results : str
        Relative path to the csv file with the parameters resulting from the fit.
    path_to_fit_params_save : str
        Relative path to the csv file with the saved fit parameters.
    path_to_file : str
        Relative path to the scan.
    path_to_peaks_params_default : str
        Relative path to the csv file with default peaks parameters.
    path_to_peaks_params_save : str
        Relative path to the csv file with the saved peaks parameters.
    peaks_params : ndarray
        2D array of str containing the peak parameters.
    peaks_params_filled : ndarray
        2D array of str containing the peak parameters + extra lines filled with ''.
    peaks_params_to_add : ndarray
        2D array of str containing the peak parameters to import from the database.
    report : ndarray
        1D array of str containing the text for the report.
    result : object lmfit.MinimizerResult
        Result of the fit.
    save_dir : str
        Directory where the data will be saved.
    SCF_Si : ndarray
        2D array of float containing the energy, f1, f2 from CXRO for Si.
    SDD_elems_available : list of str
        Indices of SDD elems available (e.g. ['0','1','2','3']).
    SDD_elems_chosen_int : list of int
        Indices of SDD selected (e.g. [0, 1, 2])
    selected_element : str
        Name of the element chosen in the database.
    selected_line : str
        Name of the line chosen in the database.
    sensorsRelTimestamps : ndarray
        1D array of float containing the time stamp of all the spectrum in the nexus file.
    session_id : str
        Session ID based on time.
    SF0 : float
        Shelf fraction of the peak.
    shelf_part : ndarray
        1D array of float containing the shelf part of a spectrum fit.
    sl : float
        Linear part of the baseline = sl*eVs+ct.
    spectrum_model : ndarray
        1D array of float containing the spectrum fit.
    spectrum_to_fit : ndarray
        1D array of float containing the spectrum to fit.
    spectrums : ndarray
        2D array of float corresponding to the subset of selected spectrums.
    spectrums_sum : ndarray
        1D array of float corresponding to the sum of all spectrums over time.
    stamps0D : ndarray
        2D array of list containing the stamps of the sensors recorded during the scan.
    tail_part : ndarray
        1D array of float containing the tail part of a spectrum fit.
    TF0 : float
        Tail fraction of the peak.
    transition_names : list of str
        List with all possible transition names from xraylib database.
    TW0 : float
        Tail width of the peak.


    Methods
    -------
    __init__(notebook_name, save_dir, files_dir):
        Constructor.
    add_str_to_logs(wid, str_to_add):
        Add a string to the logs and update the logs window.
    check_and_init():
        Check if files and folders exist, then create the interactive cell.
    create_fit_params_results(wid, path_to_fit_params_results):
        Create the csv file for the results of the fit and write its header.
    export_nb_to_pdf(wid):
        Export the notebook to pdf using a command line through the OS.
    generate_report(wid):
        Generate the text for the report.
    get_and_save_all_params_in_files(wid):
        Get the parameters from the widgets and save them in files.
    plot_extraction(wid, is_range_spectrums_empty):
        Plot the extraction with the current set of parameters.
    plot_peaks(wid):
        Plot the peaks with the current set of parameters.
    plot_single_fit(wid, title):
        Plot the fit (data, fitting curve, residual).
    plot_all_fit_results(fit_index, path_to_result_folder):
        Function to call each individual plotting subfunctions when printing the report.
    plot_fit_areas_from_file(path_to_fit_params_results):
        Plot the result areas of a fit after loading it from file.
    plot_fit_curve_from_file(path_to_fit_folder, fit_index):
        Plot the result curve of a fit after loading it from file.
    plot_fit_parameters_from_file(path_to_fit_params_save, path_to_fit_params_results):
        Plot the result parameters of a fit after loading it from file.
    plot_fit_positions_from_file(path_to_fit_params_results):
        Plot the result peaks positions of a fit after loading it from file.
    run_single_fit(wid):
        Run the fit.
    save_extract_params_in_file(wid, path_to_extract_params):
        Save extraction parameters in a csv file.
    save_fit_curves_results_in_file(wid, path_to_fit_curve_results):
        Save fit curve in a csv file.
    save_fit_logs_in_file(wid, path_to_fit_log_results, spectrum_index)
         Save fit logs in a txt file.
    save_fit_params_in_file(wid, path_to_fit_params):
        Save fit parameters in a csv file.
    save_fit_params_results_in_file(wid, path_to_fit_params_results, spectrum_index):
        Save fit results in a csv file.
    save_peaks_params_in_file(wid, path_to_peaks_params):
        Save peaks parameters in a csv file.
    set_elements(wid):
        Extract elements/lines/peaks from the sheet.
    set_extract_params_from_file(wid, path_to_extract_params):
        Load extraction parameters from a csv file.
    set_extract_params_from_widgets(wid):
        Set the extraction values of expt from the current values of the widgets.
    set_fit_params_from_file(wid, path_to_fit_params):
        Load fit parameters from a csv file.
    set_fit_params_from_widgets(wid):
        Set the fit values of expt from the current values of the widgets.
    set_list_extract_params_files(wid):
        Set the list of csv files containing extraction parameters.
    set_list_fit_params_files(wid):
        Set the list of csv files containing fit parameters.
    set_list_files(wid):
        Set the list of scans available.
    set_list_peaks_params_files(wid):
        Set the list of csv files containing peaks parameters.
    set_paths(wid):
        Set the path to the different saving files.
    set_peaks_params_from_file(wid, path_to_peaks_params):
        Load peaks parameters from a csv file.
    set_peaks_params_from_sheet(wid, sheet):
        Set the peaks values of expt from the current values of the sheet.
    set_peaks_params_from_widgets(wid):
        Set the peak values of expt from the current values of the widgets.
    set_peaks_relative_strength(wid):
        Set the relative strength of each peak relative to the most intense one within the line.
    set_result_params_to_fit_output(wid):
        Update elements_fit and eVs_fit with the result of the fit, and compute the fitting curve with its sub-parts.
    set_result_params_to_nan(wid):
        Update the results of the fit (params and curves) with NaNs.
    set_scan_info(wid):
        Extract SDD elems available and number of points in the scan.
    set_session_id(wid):
        Set the session ID based on time.
    validate_sheet(wid):
        Validate the sheet in the peak tab.

    Classes
    -------
    Element:
        Class containing lines/peaks parameters relative to an element.
    '''
    def __init__(self, notebook_name, save_dir, files_dir, logs_lvl):
        '''
        Constructor.

        Parameters
        ----------
        notebook_name : str
            The name of the notebook (necessary for saving in pdf).
        save_dir : str
            The directory where the data will be saved.
        files_dir : str
            The directory where the nexus files are.
        logs_lvl : int
            Depth level of the logs.
        '''
        self.notebook_name = notebook_name
        self.save_dir = save_dir
        self.files_dir = files_dir
        self.path_to_extract_params_default = 'lib/params_default/extract_params_default.csv'
        self.path_to_peaks_params_default = 'lib/params_default/peaks_params_default.csv'
        self.path_to_fit_params_default = 'lib/params_default/fit_params_default.csv'
        self.path_to_db = 'lib/frontend/xraylib_lines.pro'

        self.logs = ''
        self.logs_lvl = logs_lvl

        # Construct a list with all possible transition names from xraylib database
        # transition_name = ['KL1', 'KL2', 'KL3', 'KM1', ...]
        self.transition_names = []
        with open(self.path_to_db, "r") as f:
            csvreader = csv.reader(f)
            for row in csvreader:
                if row!=[]:
                    if (row[0][0]!=';' and row[0]!='end'):
                        self.transition_names.append(row[0].split(' = ')[0].split('_')[0])

        self.all_spectrums = np.array([])
        self.baseline = np.array([])
        self.bckg_eVs_inf = 0.
        self.bckg_eVs_sup = 1.
        self.beam_energy = 0.
        self.broad_factor = 0.
        self.compton_part = np.array([])
        self.ct = 0.
        self.channels = np.array([])
        self.current_sensorsRelTimestamps = np.array([])
        self.data0D = np.array([])
        self.elements = None
        self.elements_fit = None
        self.eV0 = 0.
        self.eVs = np.array([])
        self.eVs_fit = np.array([])
        self.first_channel = 0
        self.first_spectrum = 0
        self.gain = 0.
        self.gaussian_part = np.array([])
        self.highTF0 = 0.
        self.highTW0 = 0.
        self.is_bckg_subset = False
        self.is_broad_factor = False
        self.is_ct = False
        self.is_eV0 = False
        self.is_extract_done = False
        self.is_fit_done = False
        self.is_fit_fail = False
        self.is_gain = False
        self.is_highTF0 = False
        self.is_highTW0 = False
        self.is_last_fit_a_preview = False
        self.is_lowTF0 = False
        self.is_lowTW0 = False
        self.is_noise = False
        self.is_SDD_elem_0 = False
        self.is_SDD_elem_1 = False
        self.is_SDD_elem_2 = False
        self.is_SDD_elem_3 = False
        self.is_SDD_elem_4 = False
        self.is_SF0 = False
        self.is_sheet_loaded = False
        self.is_sl = False
        self.is_spectrum_empty = False
        self.is_TF0 = False
        self.is_TW0 = False
        self.last_channel = 0
        self.last_non_zero_spectrum = 0
        self.last_spectrum = 0
        self.list_extract_params_files = []
        self.list_fit_params_files = []
        self.list_isfit = []
        self.list_isfit_str = ''
        self.list_files_filenames = []
        self.list_peaks_params_files = []
        self.lowTF0 = 0.
        self.lowTW0 = 0.
        self.min_strength = 0.
        self.name = ''
        self.nb_allspectrums = 0
        self.nb_spectrums = 0
        self.noise = 0.
        self.filename = ''
        self.path_to_extract_params = ''
        self.path_to_extract_params_save = ''
        self.path_to_fit_curve_results = ''
        self.path_to_fit_folder = ''
        self.path_to_fit_log_results = ''
        self.path_to_fit_params = ''
        self.path_to_fit_params_results = ''
        self.path_to_fit_params_save = ''
        self.path_to_file = ''
        self.path_to_peaks_params_save = ''
        self.peaks_params = np.array([])
        self.peaks_params_filled = np.array([])
        self.peaks_params_to_add = np.array([])
        self.report = np.array([])
        self.result = None
        self.SCF_Si = np.genfromtxt('lib/fit/f-Si') #Requires the file f-si from CXRO.
        self.SDD_elems_available = []
        self.SDD_elems_chosen_int = []
        self.selected_element = ''
        self.selected_line = ''
        self.sensorsRelTimestamps = np.array([])
        self.session_id = ''
        self.SF0 = 0.
        self.shelf_part = np.array([])
        self.sl = 0.
        self.spectrum_model = np.array([])
        self.spectrum_to_fit = np.array([])
        self.spectrums = np.array([])
        self.spectrums_sum = np.array([])
        self.stamps0D = np.array([])
        self.tail_part = np.array([])
        self.TF0 = 0.
        self.TW0 = 0.


    def check_and_init(self):
        '''
        Check if files and folders exist, then create the interactive cell.

        Raises
        ------
        SystemExit('Save directory not found.')
            when save directory not found

        SystemExit('Files directory not found.')
            when files directory not found

        SystemExit('Notebook file not found.')
            when notebook file not found
        '''
        if os.path.exists(self.save_dir):
            print("Results will be saved in the directory:\n%s"%self.save_dir)
        else:
            print(_RED+'Careful, the directory for saving the data was not found.'+_RESET)
            print('Save directory indicated in the first cell: %s'%self.save_dir)
            sys.exit('Save directory not found.')

        if os.path.exists(self.files_dir):
            print("Scans (nexus files) should be in the directory:\n%s"%self.files_dir)
        else:
            print(_RED+"Careful, the directory where the scans are stored was not found."+_RESET)
            print('Files directory indicated in the first cell: %s'%self.files_dir)
            sys.exit('Files directory not found.')

        if not os.path.exists(self.notebook_name):
            print(_RED+"Careful, assign the correct notebook name to self.notebook_name."+_RESET)
            print('Notebook name indicated in the first cell: %s'%self.files_dir)
            sys.exit('Notebook file not found.')

        # Set the tokens
        self.is_extract_done = False
        self.is_sheet_loaded = False

        # Create the interactive cell
        nb.create_cell(code='FE.start_session(expt, wid)', position ='at_bottom', celltype='code', is_print=False)



    def export_nb_to_pdf(self, wid):
        '''
        Export the notebook to pdf using a command line through the OS.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.

        Returns
        -------
        bool
            export_done, True if the export suceeded without error/warning
        '''
        # Save the current state of the notebook (including the widgets)
        nb.save()

        # Export the pdf
        t0 = time.time()
        rc = 1
        while rc>0:
            if (time.time()-t0) > 100:
                # Timeout before PDF export is considered as failed
                export_done = False
                break

            time.sleep(3)
            command = 'jupyter nbconvert '
            command+= self.notebook_name
            command+= ' --to pdf '
            command+= '  --TagRemovePreprocessor.remove_cell_tags \'notPrint\' ' # Remove the widgets from the PDF
            command+= ' --no-input ' # Remove the code cells
            rc = subprocess.call(command,shell=True)
            if rc==0:
                export_done = True

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Export to pdf.'
            self.add_str_to_logs(wid, str_to_add)

        return export_done


    def set_scan_info(self, wid):
        '''
        Extract SDD elems available and number of points in the scan.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''

        if '.nxs' in self.filename:

            nexus = PN.PyNexusFile(self.path_to_file)

            # Extract list of detector elements available
            stamps = nexus.extractStamps()
            SDD_elems_available = []
            for stamp in stamps:
                if (stamp[1] is not None and "fluospectrum0" in stamp[1].lower()):
                    SDD_elems_available.append(stamp[1].lower()[-1])

            # Extract number of spectrums taken during the scan
            nb_allspectrums = int(nexus.get_nbpts())

            self.SDD_elems_available = SDD_elems_available
            self.nb_allspectrums = nb_allspectrums

            # Update logs
            if self.logs_lvl>=1:
                str_to_add = 'Extract information from nexus file.'
                self.add_str_to_logs(wid, str_to_add)

        if '.dat' in self.filename:

            # Extract list of detector elements available
            SDD_elems_available = []
            nb_allspectrums = 0
            for index_element in [0,1,2,3,4]:
                path_to_mat = self.path_to_file[:-4]+'_fluospectrum0'+str(index_element)+'.mat'

                if os.path.isfile(path_to_mat):
                    SDD_elems_available.append(str(index_element))

                    # Extract number of spectrums taken during the scan
                    nb_allspectrums = np.shape(np.genfromtxt(path_to_mat))[0]

            self.SDD_elems_available = SDD_elems_available
            self.nb_allspectrums = nb_allspectrums

            # Update logs
            if self.logs_lvl>=1:
                str_to_add = 'Extract information from mat file.'
                self.add_str_to_logs(wid, str_to_add)

    def plot_extraction(self, wid, is_range_spectrums_empty=False):
        '''
        Plot the extraction with the current set of parameters.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        is_range_spectrums_empty : bool, optional
            True if the selected range of spectrums appears to be empty.
        '''
        # Plot all the spectrums (stopping at the last non-zero one)
        fig = plt.figure(figsize=(12,6))
        ax1 = fig.add_subplot(111)
        ax1.set_title('All the spectrums in the file')
        ax1.set(xlabel = 'spectrum index', ylabel = 'channel')
        ax1.set_xlim(left = -1, right = self.last_non_zero_spectrum+1)
        ax1.axvline(self.first_spectrum, linestyle = '--', color = 'k', linewidth = 3,\
                    label = 'Selected spectrum range')
        ax1.axvline(self.last_spectrum, linestyle = '--', color = 'k', linewidth = 3)
        ax1.imshow(self.all_spectrums.transpose(), cmap = 'viridis', aspect = 'auto', norm=mplcolors.LogNorm())
        plt.legend(fontsize=12)


        # Plot the whole channel range
        fig = plt.figure(figsize=(12,8))
        ax1 = fig.add_subplot(211)
        ax1.set_title('Selected channel range on the sum of all spectrums')
        ax1.set(xlabel = 'channel', ylabel = 'counts')
        ax1.axvline(self.first_channel, linestyle = '--', color = 'r', linewidth = 3, label = 'Selected channel range')
        ax1.axvline(self.last_channel, linestyle = '--', color = 'r', linewidth = 3)
        ax1.plot(np.arange(2048), self.all_spectrums.sum(axis = 0), 'k.-')
        ax1.legend()
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2 = fig.add_subplot(212)
        ax2.set(xlabel = 'channel', ylabel = 'counts')
        ax2.axvline(self.first_channel, linestyle = '--', color = 'r', linewidth = 3)
        ax2.axvline(self.last_channel, linestyle = '--', color = 'r', linewidth = 3)
        ax2.plot(np.arange(2048), self.all_spectrums.sum(axis = 0), 'k.-')
        ax2.set_yscale('log')
        ax2.set_ylim(bottom = 1)
        yticks = ax1.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.subplots_adjust(hspace=.0)

        if not is_range_spectrums_empty:
            #Plot the selected spectrum range
            fig = plt.figure(figsize=(12,6))
            ax1 = fig.add_subplot(111)
            ax1.set_title('Zoom on subset of spectrums [%g:%g]'%(self.first_spectrum,self.last_spectrum))
            ax1.set(xlabel = 'spectrum index', ylabel = 'channel')
            ax1.imshow(self.spectrums.transpose(), cmap = 'viridis', aspect = 'auto', norm=mplcolors.LogNorm(),
                      interpolation='none',
                      extent=[self.first_spectrum-0.5,self.last_spectrum+0.5,
                              self.last_channel+0.5,self.first_channel-0.5])
            ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax1.yaxis.set_major_locator(MaxNLocator(integer=True))

            #Plot the selected channel range
            fig = plt.figure(figsize=(12,8))
            ax1 = fig.add_subplot(211)
            ax1.set_title('Zoom on subset of channels [%g:%g]'%(self.first_channel,self.last_channel))
            ax1.set(xlabel = 'channel', ylabel = 'counts')
            ax1.plot(self.channels, self.spectrums[0], 'r-', label = 'Spectrum %g'%self.first_spectrum)
            ax1.plot(self.channels, self.spectrums[-1], 'b-', label = 'Spectrum %g'%self.last_spectrum)
            ax1.legend(fontsize=12)
            plt.setp(ax1.get_xticklabels(), visible=False)

            ax2 = fig.add_subplot(212)
            ax2.set(xlabel = 'channel', ylabel = 'counts')
            ax2.plot(self.channels, self.spectrums[0], 'r-')
            ax2.plot(self.channels, self.spectrums[-1], 'b-')
            ax2.set_yscale('log')
            ax2.set_ylim(bottom = 1)
            yticks = ax1.yaxis.get_major_ticks()
            yticks[-1].label1.set_visible(False)
            ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
            plt.subplots_adjust(hspace=.0)

        plt.show()

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Plot of the extraction.'
            self.add_str_to_logs(wid, str_to_add)


    def set_extract_params_from_file(self, wid, path_to_extract_params):
        '''
        Load extraction parameters from a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_extract_params : str
            Path to the csv file.
        '''
        if not os.path.exists(path_to_extract_params):
            str_to_add = 'The file %s was not found.'%path_to_extract_params

        else:
            with open(path_to_extract_params, "r") as f:
                reader = csv.DictReader(f, delimiter=';',dialect='excel')
                for row in reader:
                    self.is_SDD_elem_0 = literal_eval(row['#is_SDD_elem_0'])
                    self.is_SDD_elem_1 = literal_eval(row['#is_SDD_elem_1'])
                    self.is_SDD_elem_2 = literal_eval(row['#is_SDD_elem_2'])
                    self.is_SDD_elem_3 = literal_eval(row['#is_SDD_elem_3'])
                    self.is_SDD_elem_4 = literal_eval(row['#is_SDD_elem_4'])
                    self.first_channel = int(row['#first_channel'])
                    self.last_channel = int(row['#last_channel'])
                    self.first_spectrum = int(row['#first_spectrum'])
                    self.last_spectrum = int(row['#last_spectrum'])

            str_to_add = 'Extraction parameters imported from:\n%s'%path_to_extract_params

        # Update logs
        if self.logs_lvl>=0:
            self.add_str_to_logs(wid, str_to_add)

    def set_fit_params_from_file(self, wid, path_to_fit_params):
        '''
        Load fit parameters from a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_fit_params : str
            Path to the csv file.
        '''
        if not os.path.exists(path_to_fit_params):
            str_to_add = 'The file %s was not found.'%path_to_fit_params

        else:
            with open(path_to_fit_params, "r") as f:
                reader = csv.DictReader(f, delimiter=';',dialect='excel')
                for row in reader:
                    self.list_isfit_str = str(row['#list_isfit_str'])
                    self.gain = float(row['#gain'].replace(',', '.'))
                    self.eV0 = float(row['#eV0'].replace(',', '.'))
                    self.sl = float(row['#sl'].replace(',', '.'))
                    self.ct = float(row['#ct'].replace(',', '.'))
                    self.noise = float(row['#noise'].replace(',', '.'))
                    self.SF0 = float(row['#SF0'].replace(',', '.'))
                    self.TF0 = float(row['#TF0'].replace(',', '.'))
                    self.TW0 = float(row['#TW0'].replace(',', '.'))
                    self.broad_factor = float(row['#broad_factor'].replace(',', '.'))
                    self.lowTF0 = float(row['#lowTF0'].replace(',', '.'))
                    self.highTF0 = float(row['#highTF0'].replace(',', '.'))
                    self.lowTW0 = float(row['#lowTW0'].replace(',', '.'))
                    self.highTW0 = float(row['#highTW0'].replace(',', '.'))
                    self.is_bckg_subset = literal_eval(row['#is_bckg_subset'])
                    self.bckg_eVs_inf = float(row['#bckg_eVs_inf'].replace(',', '.'))
                    self.bckg_eVs_sup = float(row['#bckg_eVs_sup'].replace(',', '.'))

            # convert list_isfit_str into a list
            self.list_isfit = self.list_isfit_str.split(',')

            str_to_add = 'Fit parameters imported from:\n%s'%path_to_fit_params

        # Update logs
        if self.logs_lvl>=0:
            self.add_str_to_logs(wid, str_to_add)

    def set_peaks_params_from_file(self, wid, path_to_peaks_params):
        '''
        Load peaks parameters from a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_peaks_params : str
            Path to the csv file.
        '''
        if not os.path.exists(path_to_peaks_params):
            # Update logs
            str_to_add = 'The file %s was not found.'%path_to_peaks_params
            self.add_str_to_logs(wid, str_to_add)

        else:

            peaks_params = np.array([])
            with open(path_to_peaks_params, "r") as f:
                csvreader = csv.reader(f, delimiter=';',dialect='excel')
                # First line is the header
                peaks_header = next(csvreader)
                nb_columns = len(peaks_header)
                for row in csvreader:
                    peaks_params = np.append(peaks_params, row)
            peaks_params = np.reshape(peaks_params, (len(peaks_params)//nb_columns,nb_columns))

            # Extract in peaks_params only the params which will go in the sheet
            self.peaks_params = peaks_params[:,:6]

            # Get the other general parameters
            self.beam_energy = float(peaks_params[0,6].replace(',', '.'))
            self.min_strength = float(peaks_params[0,7].replace(',', '.'))
            self.gain = float(peaks_params[0,8].replace(',', '.'))
            self.eV0 = float(peaks_params[0,9].replace(',', '.'))

            # Update logs
            if self.logs_lvl>=0:
                str_to_add = 'Peaks parameters loaded from the file:\n%s'%path_to_peaks_params
                self.add_str_to_logs(wid, str_to_add)

    def save_peaks_params_in_file(self, wid, path_to_peaks_params):
        '''
        Save peaks parameters in a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_peaks_params : str
            Path to the csv file.
        '''
        # Create a folder for saving all data related to the scan
        if not os.path.exists(self.save_dir+self.name):
            os.mkdir(self.save_dir+self.name)

        # Create a subfolder corresponding to the current session
        if not os.path.exists(self.save_dir+self.name+'/'+self.session_id):
            os.mkdir(self.save_dir+self.name+'/'+self.session_id)

        # Write to the csv file
        header = ['#Peak name', '#Transition name', '#Position (eV)', '#Strength',
                  '#Fit position?', '#Fit peak?', '#Beam energy (eV)', '#Min strength',
                  '#Gain', '#eV0']

        with open(path_to_peaks_params, "w", newline='') as f:
            writer = csv.writer(f,delimiter=';',dialect='excel')
            writer.writerow(header)
            for row in self.peaks_params:
                writer.writerow(np.append(
                                    row,
                                    [self.beam_energy, self.min_strength, self.gain, self.eV0]
                                    )
                                )

        # Update logs
        if self.logs_lvl>=0:
            str_to_add = 'Current peaks parameters saved in:\n%s'%path_to_peaks_params
            self.add_str_to_logs(wid, str_to_add)

    def save_fit_params_in_file(self, wid, path_to_fit_params):
        '''
        Save fit parameters in a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_fit_params : str
            Path to the csv file.
        '''
        # Create a folder for saving all data related to the scan
        if not os.path.exists(self.save_dir+self.name):
            os.mkdir(self.save_dir+self.name)

        # Create a subfolder corresponding to the current session
        if not os.path.exists(self.save_dir+self.name+'/'+self.session_id):
            os.mkdir(self.save_dir+self.name+'/'+self.session_id)


        # Write to the csv file
        with open(path_to_fit_params, "w", newline='') as f:
            writer = csv.writer(f,delimiter=';',dialect='excel')
            header = np.array([
                    '#list_isfit_str',
                    '#gain',
                    '#eV0',
                    '#sl',
                    '#ct',
                    '#noise',
                    '#SF0',
                    '#TF0',
                    '#TW0',
                    '#broad_factor',
                    '#lowTF0',
                    '#highTF0',
                    '#lowTW0',
                    '#highTW0',
                    '#is_bckg_subset',
                    '#bckg_eVs_inf',
                    '#bckg_eVs_sup'
                    ])
            writer.writerow(header)

            writer.writerow([
                    self.list_isfit_str,
                    self.gain,
                    self.eV0,
                    self.sl,
                    self.ct,
                    self.noise,
                    self.SF0,
                    self.TF0,
                    self.TW0,
                    self.broad_factor,
                    self.lowTF0,
                    self.highTF0,
                    self.lowTW0,
                    self.highTW0,
                    self.is_bckg_subset,
                    self.bckg_eVs_inf,
                    self.bckg_eVs_sup
                    ])

        # Update logs
        if self.logs_lvl>=0:
            str_to_add = 'Current fit parameters saved in:\n%s'%path_to_fit_params
            self.add_str_to_logs(wid, str_to_add)

    def set_extract_params_from_widgets(self, wid):
        '''
        Set the extraction values of expt from the current values of the widgets.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.is_SDD_elem_0 = wid.is_SDD_elem_0.value
        self.is_SDD_elem_1 = wid.is_SDD_elem_1.value
        self.is_SDD_elem_2 = wid.is_SDD_elem_2.value
        self.is_SDD_elem_3 = wid.is_SDD_elem_3.value
        self.is_SDD_elem_4 = wid.is_SDD_elem_4.value
        self.first_channel = wid.first_channel.value
        self.last_channel = wid.last_channel.value
        self.first_spectrum = wid.first_spectrum.value
        self.last_spectrum = wid.last_spectrum.value

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'expt updated with extraction params from the widgets.'
            self.add_str_to_logs(wid, str_to_add)

    def set_fit_params_from_widgets(self, wid):
        '''
        Set the fit values of expt from the current values of the widgets.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.gain = wid.gain.value
        self.eV0 = wid.eV0.value
        self.sl = wid.sl.value
        self.ct = wid.ct.value
        self.noise = wid.noise.value
        self.SF0 = wid.SF0.value
        self.TF0 = wid.TF0.value
        self.TW0 = wid.TW0.value
        self.broad_factor = wid.broad_factor.value
        self.lowTF0 = wid.lowTF0.value
        self.highTF0 = wid.highTF0.value
        self.lowTW0 = wid.lowTW0.value
        self.highTW0 = wid.highTW0.value

        self.is_gain = wid.is_gain.value
        self.is_eV0 = wid.is_eV0.value
        self.is_sl = wid.is_sl.value
        self.is_ct = wid.is_ct.value
        self.is_noise = wid.is_noise.value
        self.is_SF0 = wid.is_SF0.value
        self.is_TF0 = wid.is_TF0.value
        self.is_TW0 = wid.is_TW0.value
        self.is_broad_factor = wid.is_broad_factor.value
        self.is_lowTF0 = wid.is_lowTF0.value
        self.is_highTF0 = wid.is_highTF0.value
        self.is_lowTW0 = wid.is_lowTW0.value
        self.is_highTW0 = wid.is_highTW0.value

        self.list_isfit = ['gain'*self.is_gain, 'eV0'*self.is_eV0,
                      'sl'*self.is_sl, 'ct'*self.is_ct, 'noise'*self.is_noise,
                      'SF0'*self.is_SF0, 'TF0'*self.is_TF0,
                      'TW0'*self.is_TW0, 'broad_factor'*self.is_broad_factor,
                      'lowTF0'*self.is_lowTF0, 'highTF0'*self.is_highTF0,
                      'lowTW0'*self.is_lowTW0, 'highTW0'*self.is_highTW0]

        # Remove trailing ""
        while "" in self.list_isfit:
            self.list_isfit.remove("")

        self.list_isfit_str = ','.join(self.list_isfit)

        self.is_bckg_subset = wid.is_bckg_subset.value
        self.bckg_eVs_inf = wid.bckg_eVs_inf.value
        self.bckg_eVs_sup = wid.bckg_eVs_sup.value     
        
        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'expt updated with fit params from the widgets.'
            self.add_str_to_logs(wid, str_to_add)

    def save_extract_params_in_file(self, wid, path_to_extract_params):
        '''
        Save extraction parameters in a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_extract_params : str
            Path to the csv file.
        '''
        # Create a folder for saving all data related to the scan
        if not os.path.exists(self.save_dir+self.name):
            os.mkdir(self.save_dir+self.name)

        # Create a subfolder corresponding to the current session
        if not os.path.exists(self.save_dir+self.name+'/'+self.session_id):
            os.mkdir(self.save_dir+self.name+'/'+self.session_id)

        # Write to the csv file
        with open(path_to_extract_params, "w", newline='') as f:
            writer = csv.writer(f,delimiter=';',dialect='excel')
            header = np.array([
                    '#is_SDD_elem_0',
                    '#is_SDD_elem_1',
                    '#is_SDD_elem_2',
                    '#is_SDD_elem_3',
                    '#is_SDD_elem_4',
                    '#first_channel',
                    '#last_channel',
                    '#first_spectrum',
                    '#last_spectrum',
                    ])
            writer.writerow(header)

            writer.writerow([
                    self.is_SDD_elem_0,
                    self.is_SDD_elem_1,
                    self.is_SDD_elem_2,
                    self.is_SDD_elem_3,
                    self.is_SDD_elem_4,
                    self.first_channel,
                    self.last_channel,
                    self.first_spectrum,
                    self.last_spectrum
                    ])

        # Update logs
        if self.logs_lvl>=0:
            str_to_add = 'Current extraction parameters saved in:\n%s'%path_to_extract_params
            self.add_str_to_logs(wid, str_to_add)


    def set_peaks_params_from_widgets(self, wid):
        '''
        Set the peak values of expt from the current values of the widgets.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.beam_energy = wid.beam_energy.value
        self.min_strength = wid.min_strength.value
        self.gain = wid.gain.value
        self.eV0 = wid.eV0.value

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'expt updated with peaks params from the widgets.'
            self.add_str_to_logs(wid, str_to_add)

    def set_peaks_params_from_sheet(self, wid, sheet):
        '''
        Set the peaks values of expt from the current values of the sheet.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        sheet : object Sheet
            Object from the class Sheet (module ipysheet)
        '''
        # Get the peaks from the sheet
        self.peaks_params_filled = ipysheet.numpy_loader.to_array(ipysheet.easy.current())

        # Remove the empty lines
        self.peaks_params = self.peaks_params_filled[np.where(self.peaks_params_filled[:,0]!='')]

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'expt updated with peaks params from the sheet.'
            self.add_str_to_logs(wid, str_to_add)

    def validate_sheet(self, wid):
        '''
        Validate the sheet in the peak tab.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''

        # Update expt with peaks params from sheet and widgets
        self.set_peaks_params_from_sheet(wid, wid.sheet)
        self.set_peaks_params_from_widgets(wid)

        # Create the elements/lines/peaks
        self.set_elements(wid)

        # Set the strength of each peak relative to the maximum strength within the same line
        self.set_peaks_relative_strength(wid)

        # Convert channels to eVs
        self.eVs = self.gain*self.channels + self.eV0

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Sheet validated.'
            self.add_str_to_logs(wid, str_to_add)

    def set_session_id(self, wid):
        '''
        Set the session ID based on time.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.session_id = datetime.now().strftime('%Y%m%d_%H%M%S')

        # Update logs
        if self.logs_lvl>=0:
            str_to_add = 'Set the session id to %s.'%(self.session_id)
            self.add_str_to_logs(wid, str_to_add)

    def set_paths(self, wid):
        '''
        Set the paths to the different saving files.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.path_to_fit_params_results = self.save_dir+self.name+'/'+self.session_id+'/fit_results.csv'
        self.path_to_extract_params_save = self.save_dir+self.name+'/'+self.session_id+'/'+'extract_params.csv'
        self.path_to_peaks_params_save = self.save_dir+self.name+'/'+self.session_id+'/'+'peaks_params.csv'
        self.path_to_fit_params_save = self.save_dir+self.name+'/'+self.session_id+'/'+'fit_params.csv'
        self.path_to_fit_folder = self.save_dir+self.name+'/'+self.session_id+'/fit_curves/'
        self.path_to_fit_log_results = self.save_dir+self.name+'/'+self.session_id+'/fit_results.log'

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Set the paths to saving files.'
            self.add_str_to_logs(wid, str_to_add)

    def set_list_files(self, wid):
        '''
        Set the list of scans available.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.list_files_filenames = [file for file in sorted(os.listdir(self.files_dir))
                                   if ('.nxs' in file or '.dat' in file)][::-1]

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Set the list of files.'
            self.add_str_to_logs(wid, str_to_add)

    def set_list_extract_params_files(self, wid):
        '''
        Set the list of csv files containing extraction parameters.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.list_extract_params_files = ['Default extraction parameters']
        tmp_files = []
        for root, _, files in os.walk(self.save_dir, topdown=True):
            for name in files:
                if 'extract_params.csv' in name:
                    path_to_csv = os.path.join(root.split('/')[-2],root.split('/')[-1])
                    tmp_files.append(path_to_csv)
        tmp_files.sort(reverse=True)
        self.list_extract_params_files += tmp_files

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Set the list of files with extraction parameters.'
            self.add_str_to_logs(wid, str_to_add)

    def set_list_peaks_params_files(self, wid):
        '''
        Set the list of csv files containing peaks parameters.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.list_peaks_params_files = ['Default peaks parameters']
        tmp_files = []
        for root, _, files in os.walk(self.save_dir, topdown=True):
            for name in files:
                if 'peaks_params.csv' in name:
                    path_to_csv = os.path.join(root.split('/')[-2],root.split('/')[-1])
                    tmp_files.append(path_to_csv)
        tmp_files.sort(reverse=True)
        self.list_peaks_params_files += tmp_files

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Set the list of files with peaks parameters.'
            self.add_str_to_logs(wid, str_to_add)


    def set_list_fit_params_files(self, wid):
        '''
        Set the list of csv files containing fit parameters.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.list_fit_params_files = ['Default fit parameters']
        tmp_files = []
        for root, _, files in os.walk(self.save_dir, topdown=True):
            for name in files:
                if 'fit_params.csv' in name:
                    path_to_csv = os.path.join(root.split('/')[-2],root.split('/')[-1])
                    tmp_files.append(path_to_csv)
        tmp_files.sort(reverse=True)
        self.list_fit_params_files += tmp_files

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Set the list of files with fit parameters.'
            self.add_str_to_logs(wid, str_to_add)

    def set_elements(self, wid):
        '''
        Extract elements/lines/peaks from the sheet.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        # List of objects Element
        self.elements = []

        # Go line by line in the sheet
        for row in self.peaks_params:

            # Treat only the peaks which are fitted
            if row[5]!='no' and row[0]!='':

                # Set the name of the element
                current_element_name = row[0]

                # Set the name of the line by appending the name of the current element to
                # the first character of the transition name (K, L, M)
                # Or the full name if it is not a transition (Elastic peak, Compton peak ...)
                if row[1][0] in ['K', 'L', 'M']:
                    current_line_name = current_element_name+'_'+row[1][0]
                else:
                    current_line_name = current_element_name+'_'+row[1]

                # Check if the element has already been created
                is_new_element = True
                for element in self.elements:
                    if current_element_name == element.name:
                        is_new_element = False

                # Add this element to the list if it did not exist
                if is_new_element:
                    current_element = self.Element(current_element_name)
                    self.elements = np.append(self.elements, current_element)

                # Check if the line has already been created for this element
                is_new_line = True
                for line in current_element.lines:
                    if current_line_name == line.name:
                        is_new_line = False

                # Add this line to the list if it did not exist for this element
                if is_new_line:
                    current_line = current_element.Line(current_line_name)
                    current_element.lines = np.append(current_element.lines, current_line)

                # Create the peak and add it to the line
                current_peak = current_line.Peak(
                            peakName = 'peak'+'_'+row[0]+'_'+row[1],
                            peakPosition = float(row[2]),
                            peakStrength = float(row[3]),
                            peakIs_fitpos = bool(row[4]=='yes'))

                current_line.peaks = np.append(current_line.peaks, current_peak)

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Extract elements from sheet.'
            self.add_str_to_logs(wid, str_to_add)

    def set_peaks_relative_strength(self, wid):
        '''
        Set the relative strength of each peak relative to the most intense one within the line.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        for element in self.elements:
            for line in element.lines:
                max_strength = np.max([peak.strength for peak in line.peaks])
                # Normalize the strengths with the most intense one
                for peak in line.peaks:
                    peak.relative_strength = peak.strength/max_strength

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Set the peaks relative strength.'
            self.add_str_to_logs(wid, str_to_add)

    class Element:
        '''
        Class containing lines/peaks parameters relative to an element.

        Attributes
        ----------
        lines : list of Line
            List of objects Line.
        name : str
            Name of the element (i.e. 'Si', 'P', 'Cl', ...).

        Methods
        -------
        __init__(elementName):
            Constructor.


        Classes
        -------
        Line:
            Class for a fluorescence line (i.e. K, L, M).
        '''
        def __init__(self, elementName):
            '''
            Constructor.

            Parameters
            ----------
            elementName : str
                Name of the element (i.e. Si, P, Cl, ...).
            '''
            self.lines = []
            self.name = elementName

        class Line:
            '''
            Class containing peaks parameters relative to a fluorescence line of an element.

            Attributes
            ----------
            name : str
                Name of the line (i.e. 'Si_K', 'Ar_L', 'K_M', or user defined).
            peaks : list of Peak
                List of objects Peaks.
            area : float
                Area of the line.
            peak_series : ndarray
                1D array containing the part of the spectrum corresponding to the given line, during a fit.

            Methods
            -------
            __init__(lineName):
                Constructor.

            Classes
            -------
            Peak:
                Class for a single peak.
            '''
            def __init__(self, lineName):
                '''
                Constructor.

                Parameters
                ----------
                lineName : str
                    Name of the line (i.e. 'Si_K', 'Ar_L', 'K_M', or user defined).
                '''
                self.peaks = []
                self.name = lineName
                self.area = 1.

            class Peak:
                '''
                Class containing all the information relative to a single peak.

                Attributes
                ----------
                is_fit_pos : bool
                    True if the position of the peak should be a fit param.
                name : str
                    Name of peak ('peak' + element name + line name).
                position : float
                    Position of the peak in eV before fit.
                relative_strength : float
                    Strength of the peak relative to the most intense one of its line.
                strength : float
                    Non-normalized strength of the peak.

                Methods
                -------
                __init__():
                    Constructor.
                '''
                def __init__(self, peakName, peakPosition, peakStrength, peakIs_fitpos):
                    '''
                    Constructor.

                    Parameters
                    ----------
                    peakName : str
                        Name of the peak ('peak' + '_' + element name + '_' + line name).
                    peakPosition : float
                        Position of the peak in eV.
                    peakStrength : float
                        Non-normalized strength of the peak.
                    peakIs_fitpos : bool
                        True if the position of the peak should be a fit param.
                    '''
                    self.name = peakName
                    self.position = peakPosition
                    self.strength = peakStrength
                    self.is_fitpos = peakIs_fitpos


    def plot_peaks(self, wid):
        '''
        Plot the peaks with the current set of parameters.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        # Plot the whole spectrum twice (lin and log y-scale)
        fig = plt.figure(figsize=(15,8))
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True, sharey=False)

        for ax in axs:

            ax.minorticks_on()

            # Iterator to have one color/linestyle per element, in the two plots
            colors = iter(['#006BA4', '#FF800E', '#ABABAB', '#595959', 'k', '#C85200', 'b', '#A2C8EC', '#FFBC79']*200)
            linestyles = iter(['--', '-.', '-', ':']*400)

            # Plot the sum of all spectrums, and the position of the peaks
            ax.plot(self.eVs, self.spectrums_sum, 'k.')
            ax.set(xlabel = 'E (eV)', ylabel = 'counts')
            for element in self.elements:

                color = next(colors)
                linestyle = next(linestyles)

                for line in element.lines:
                    for peak in line.peaks:
                        ax.axvline(x = peak.position,  color = color,
                                    linestyle = linestyle, label = element.name)


        axs[1].set(xlabel = 'E (eV)', ylabel = 'counts')
        axs[1].set_ylim(bottom = 1)
        axs[1].yaxis.set_label_position("right")
        axs[1].yaxis.tick_right()
        axs[1].set_yscale('log')

        # Avoid having multiple times the same label in legend
        handles, labels = axs[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))

        axs[0].legend(by_label.values(), by_label.keys(), bbox_to_anchor=(0., 1.02, 1., .102), loc='lower center',
                   ncol=8, borderaxespad=0.)

        plt.show()

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Plot the peaks.'
            self.add_str_to_logs(wid, str_to_add)


    def plot_single_fit(self, wid=None, title='', is_clear_output=True):
        '''
        Plot the fit (data, fitting curve, residual).

        Parameters
        ----------
        wid : object myWidgets, optional
            Object from the class myWidgets.
        title : str, optional
            Title of the plot.
        is_clear_output : bool, optional
            Clear the previous output or not (used for refreshing during a series of fits).
        '''
        # If gain or eV0 were fitted, the eVs need to be updated
        if ('gain' or 'eV0') in self.list_isfit:
            eVs = self.eVs_fit
        else:
            eVs = self.eVs

        # Plot the whole spectrum twice (lin and log y-scale)
        if is_clear_output:
            clear_output(wait=True)
        fig = plt.figure(figsize=(15,10))
        fig.suptitle(title, fontsize=14)
        fig.subplots_adjust(top=0.95)
        gs = fig.add_gridspec(3, hspace=0, height_ratios=[0.4,0.4,0.2])
        axs = gs.subplots(sharex=True, sharey=False)

        for ax in axs[0:2]:

            ax.minorticks_on()

            # Plot the sum of all spectrums, and the position of the peaks
            ax.plot(eVs, self.spectrum_to_fit, 'k.', label='Data')
            ax.plot(eVs, self.spectrum_model, 'r-', label='Fit')
            if wid is not None:
                if wid.show_gaussian.value:
                    ax.plot(eVs, self.gaussian_part, 'g--', label='Gaussian part')
                if wid.show_shelf.value:
                    ax.plot(eVs, self.shelf_part, 'r-.', label = 'Shelf part')
                if wid.show_tail.value:
                    ax.plot(eVs, self.tail_part, 'm--', label = 'Tail part')
                if wid.show_bckg.value:
                    ax.plot(eVs, self.baseline, 'y:', label = 'Background')
                if (wid.show_compton.value and np.shape(self.compton_part)):
                    ax.plot(eVs, self.compton_part, 'c--', label = 'Compton part')
            ax.set(xlabel = 'E (eV)', ylabel = 'counts')


        axs[0].legend(fontsize=12)
        axs[1].set(xlabel = 'E (eV)', ylabel = 'counts')
        axs[1].set_ylim(bottom = 1)
        axs[1].set_yscale('log')
        axs[1].yaxis.set_label_position("right")
        axs[1].yaxis.tick_right()

        axs[2].set(xlabel = 'E (eV)', ylabel = 'residuals')
        axs[2].plot(eVs, self.spectrum_model-self.spectrum_to_fit, 'k-')

        plt.show()

        # Update logs
        if (self.logs_lvl>=1 and wid is not None):
            str_to_add = 'Plot the fit results.'
            self.add_str_to_logs(wid, str_to_add)

    def set_result_params_to_nan(self, wid):
        '''
        Update the results of the fit (params and curves) with NaNs.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        for element in self.elements_fit:
            for line in element.lines:
                line.area = np.nan
                line.stderr_area = np.nan
                for peak in line.peaks:
                    peak.position = np.nan
                    peak.stderr_position = np.nan

        for param_name in self.list_isfit:
            self.result.params[param_name].value = np.nan
            self.result.params[param_name].stderr = np.nan

        self.spectrum_model = np.nan*self.eVs
        self.gaussian_part = np.nan*self.eVs
        self.shelf_part = np.nan*self.eVs
        self.tail_part = np.nan*self.eVs
        self.baseline = np.nan*self.eVs
        self.compton_part = np.nan*self.eVs
        self.result.residual = np.nan*self.eVs

        # In case of non convergence or empty spectrum, keep the original eVs
        self.eVs_fit = self.eVs

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Results set to nans.'
            self.add_str_to_logs(wid, str_to_add)

    def set_result_params_to_fit_output(self, wid):
        '''
        Update elements_fit and eVs_fit with the result of the fit, and compute the fitting curve with its sub-parts.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        # Extract result from the fit
        for element in self.elements_fit:
            for line in element.lines:
                # It can happen that the fit converges, but do not return stderr
                if line.area.stderr is None:
                    line.stderr_area = np.nan
                else:
                    line.stderr_area = line.area.stderr
                line.area = line.area.value
                for peak in line.peaks:
                    if peak.is_fitpos:
                        if peak.position.stderr is None:
                            peak.stderr_position = np.nan
                        else:
                            peak.stderr_position = peak.position.stderr
                        peak.position = peak.position.value

        # Compute new eVs if gain and/or eV0 were fitted
        self.eVs_fit =  self.result.params['gain'].value*self.channels + self.result.params['eV0'].value

        # Compute the fitting curve and its sub-parts
        self.spectrum_model, self.gaussian_part, self.shelf_part, self.tail_part, self.baseline, self.compton_part =\
            funFit.funSpectrum(self.elements_fit, self.eVs_fit, self.SCF_Si, self.result.params['ct'].value,\
                               self.result.params['sl'].value, self.result.params['noise'].value,\
                               self.result.params['SF0'].value, self.result.params['TF0'].value, \
                               self.result.params['TW0'].value, self.result.params['broad_factor'].value,\
                               self.result.params['lowTF0'].value, self.result.params['highTF0'].value,\
                               self.result.params['lowTW0'].value, self.result.params['highTW0'].value)

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Results updated with fit output.'
            self.add_str_to_logs(wid, str_to_add)

    def run_single_fit(self, wid):
        '''
        Run the fit.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        # Initialize tokens
        self.is_spectrum_empty = False
        self.is_fit_fail = False
        self.is_last_fit_a_preview = True

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Run single fit.'
            self.add_str_to_logs(wid, str_to_add)

        # We do not want to write on self during the fit
        # But funFitSpectrum will write on elements[x].lines[x].area
        # So we use a copy of self.elements for the fit
        self.elements_fit = deepcopy(self.elements)

        # Do the fit of the spectrum
        self.result = funFit.funFitSpectrum(self.spectrum_to_fit, self.list_isfit, self.elements_fit, self.channels,
                                             self.SCF_Si, self.gain, self.eV0, self.ct, self.sl, self.noise,
                                             self.SF0, self.TF0, self.TW0, self.broad_factor,
                                             self.lowTF0, self.highTF0, self.lowTW0, self.highTW0,
                                             self.is_bckg_subset, self.bckg_eVs_inf, self.bckg_eVs_sup)

        # Check if the fit succeeded
        self.is_fit_fail = not self.result.success

        # Check if the spectrum was empty
        if np.sum(self.spectrum_to_fit)<10.:
            self.is_spectrum_empty = True

        if (self.is_spectrum_empty or self.is_fit_fail):
            # Put nans in every fit parameters and in the resulting curves
            self.set_result_params_to_nan(wid)

        else:
            # Set the self.XXX_fit params to the fit output
            # Set the spectrum_model and sub-curves to the fit output
            self.set_result_params_to_fit_output(wid)

    def save_fit_curves_results_in_file(self, wid, path_to_fit_curve_results):
        '''
        Save fit curve in a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_fit_curve_results : str
            Path to the csv file.
        '''

        # Create a subfolder for the fit results
        if not os.path.exists(self.save_dir+self.name+'/'+self.session_id+'/fit_curves/'):
            os.mkdir(self.save_dir+self.name+'/'+self.session_id+'/fit_curves/')

        # Write to the csv file
        with open(path_to_fit_curve_results, "w", newline='') as f:
            writer = csv.writer(f,delimiter=';',dialect='excel')
            header = np.array([
                    '#sensorsRelTimestamps',
                    '#eVs',
                    '#data',
                    '#fit'
                    ])
            writer.writerow(header)

            for i in range(len(self.eVs)):
                writer.writerow([
                        self.current_sensorsRelTimestamps,
                        np.round(self.eVs_fit[i],2),
                        np.round(self.spectrum_to_fit[i],2),
                        np.round(self.spectrum_model[i],2)
                        ])

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Fit curve saved in:\n%s'%path_to_fit_curve_results
            self.add_str_to_logs(wid, str_to_add)


    def create_fit_params_results(self, wid, path_to_fit_params_results):
        '''
        Create the csv file for the results of the fit and write its header.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_fit_params_results : str
            Path to the csv file.
        '''
        # Prepare the header of the csv file
        header = np.array([])

        # Data stamps
        for stamp0D in self.stamps0D:
            if stamp0D[1] is None:
                header =np.append(header, '#'+stamp0D[0])
            else:
                header =np.append(header, '#'+stamp0D[1])

        # Stamps from the fit
        header = np.append(header, '#spectrum_index')

        for element in self.elements:
            for line in element.lines:
                header = np.append(header, '#'+'area_'+line.name)
                header = np.append(header, '#stderr_'+'area_'+line.name)
                for peak in line.peaks:
                    if peak.is_fitpos:
                        header = np.append(header, '#'+'pos_'+peak.name)
                        header = np.append(header, '#stderr_'+'pos_'+peak.name)

        for param_name in self.list_isfit:
            header = np.append(header, '#'+param_name)
            header = np.append(header, '#stderr_'+param_name)

        with open(path_to_fit_params_results, "w", newline='') as f:
            writer = csv.writer(f,delimiter=';')
            writer.writerow(header)

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'File for fit parameters created here:\n%s'%path_to_fit_params_results
            self.add_str_to_logs(wid, str_to_add)


    def save_fit_params_results_in_file(self, wid, path_to_fit_params_results, spectrum_index):
        '''
        Save fit results in a csv file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_fit_params_results : str
            Path to the csv file.
        spectrum_index : int
            Index of the spectrum.
        '''
        # Array to be written
        tbw = np.array([], dtype='float')

        # Put the data0D
        for data in self.data0D:
            tbw = np.append(tbw, data[spectrum_index])

        tbw = np.append(tbw, spectrum_index)

        for element in self.elements_fit:
            for line in element.lines:
                tbw = np.append(tbw, np.round(line.area, 4))
                tbw = np.append(tbw, np.round(line.stderr_area, 4))

                for peak in line.peaks:
                    if peak.is_fitpos:
                        tbw = np.append(tbw, np.round(peak.position, 4))
                        tbw = np.append(tbw, np.round(peak.stderr_position, 4))

        for param_name in self.list_isfit:
            tbw = np.append(tbw, np.round(self.result.params[param_name].value, 4))
            # It can happen that the fit converges, but do not return stderr
            if self.result.params[param_name].stderr is None:
                tbw = np.append(tbw, np.nan)
            else:
                tbw = np.append(tbw, np.round(self.result.params[param_name].stderr, 4))

        with open(path_to_fit_params_results, "a+", newline='') as f:
            writer = csv.writer(f,delimiter=';')
            writer.writerow(tbw)

        # Update logs
        if self.logs_lvl>=0:
            str_to_add = 'Fit results saved in:\n%s'%path_to_fit_params_results
            self.add_str_to_logs(wid, str_to_add)


    def save_fit_logs_in_file(self, wid, path_to_fit_log_results, spectrum_index):
        '''
        Save fit logs in a txt file.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        path_to_fit_log_results : str
            Path to the log file.
        spectrum_index : int
            Index of the spectrum.
        '''
        with open(path_to_fit_log_results, "a+", newline='') as f:
            f.write('Fit of spectrum '+str(spectrum_index)+'\n')
            if self.is_fit_fail:
                f.write('#FIT DID NOT CONVERGE\n')
            f.write(fit_report(self.result))
            f.write('\n')
            f.write('\n')

        # Update logs
        if self.logs_lvl>=0:
            str_to_add = 'Fit logs saved in:\n%s'%path_to_fit_log_results
            self.add_str_to_logs(wid, str_to_add)

    def get_and_save_all_params_in_files(self, wid):
        '''
        Get the parameters from the widgets and save them in files.
        Do not update the parameters from extract, since this should be done by extracting the scan.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''

        # Set the paths of the different saving files
        self.set_paths(wid)

        ############################
        # Extraction parameters
        # Do not update the extract params from the widgets,
        # in case the user changed the parameters in the widget AFTER extracting the scan
        # (this would cause a mismatch between the info from the widget and the actual size of the spectrums)

        # Save
        self.save_extract_params_in_file(wid, self.path_to_extract_params_save)

        ############################
        # Peaks parameters

        if self.is_sheet_loaded:
            # Validate the sheet (update expt with widget params)
            self.validate_sheet(wid)

            # Save
            self.save_peaks_params_in_file(wid, self.path_to_peaks_params_save)

            # Update list of available params files and the widget
            self.set_list_peaks_params_files(wid)
            wid.select_peaks_params_file.options = self.list_peaks_params_files

            # Change the default option of the file selection window to the most recent one
            wid.select_peaks_params_file.value = \
                      [x for x in self.list_peaks_params_files if self.session_id in x][0]


        ############################
        # Fit parameters

        # Update expt from the widget params
        self.set_fit_params_from_widgets(wid)

        # Save
        self.save_fit_params_in_file(wid, self.path_to_fit_params_save)

        # Update list of available params files and the widget
        self.set_list_fit_params_files(wid)
        wid.select_fit_params_file.options = self.list_fit_params_files

        # Change the default option of the file selection window to the most recent one
        wid.select_fit_params_file.value = \
                  [x for x in self.list_fit_params_files if self.session_id in x][0]

        # Update logs
        if self.logs_lvl>=1:
            str_to_add = 'Plot the peaks.'
            self.add_str_to_logs(wid, str_to_add)


    def plot_fit_curve_from_file(self, path_to_fit_folder, fit_index):
        '''
        Plot the result curve of a fit after loading it from file.
        We do not use the attributes of the current expt because this function should be callable even
        after the kernel has been restarted.

        Parameters
        ----------
        path_to_fit_folder : str
            Relative path to the folder containing the fit curves.
        fit_index : int
            Index of the fitted spectrum to display.
        '''
        # Import the fitting curves from file
        fit_curve_results = np.genfromtxt(path_to_fit_folder+\
                            'spectrum_'+str(fit_index)+'.csv', delimiter=';', names=True)
        self.eVs = fit_curve_results['eVs']
        self.eVs_fit = fit_curve_results['eVs']
        self.spectrum_to_fit =  fit_curve_results['data']
        self.spectrum_model =  fit_curve_results['fit']

        # Plot the result
        self.plot_single_fit(title='Fit of spectrum '+str(fit_index), is_clear_output=False)


    def plot_fit_areas_from_file(self, path_to_fit_params_results):
        '''
        Plot the result areas of a fit after loading it from file.
        We do not use the attributes of the current expt because this function should be callable even
        after the kernel has been restarted.

        Parameters
        ----------
        path_to_fit_params_results : str
            Relative path to the csv file with the parameters resulting from the fit.
        '''
        # Import results and params
        fit_results = np.genfromtxt(path_to_fit_params_results, delimiter=';', names=True)

        # ['area_Cl_K', 'area_X_xx', 'area_Ar_K', 'area_Compton_Co']
        names_area_line = [x for x in fit_results.dtype.names if x.startswith('area_')]

        # ['Cl_K', 'X_xx', 'Ar_K', 'Compton_Co']
        names_line = [x.split('area_')[1] for x in names_area_line]

        # Total number of plots to print
        nb_plots = len(names_line)

        # Divide the plots in figures of nb_rows plots
        # To avoid overflowing the pdf with all plots on one page
        nb_rows = 5

        if nb_plots > 0:

            # Token to print the title only once
            is_title = True

            for k in range(0, nb_plots, nb_rows):

                # Print all the packs of nb_rows plots
                if (nb_plots-k)//nb_rows>0:

                    # Create a new figure for each pack of nb_rows plots
                    if k%nb_rows == 0:
                        fig, ax = plt.subplots(figsize=(15,4*nb_rows), nrows=nb_rows)

                        # ax is a float when there is only one row
                        if nb_plots == 1:
                            ax = [ax]

                    # Plot all the plots of the corresponding pack
                    for i in range(nb_rows):
                        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%g'))
                        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

                        if is_title:
                            ax[i].set_title('AREAS OF FLUORESCENCE LINES', pad = 15, y=1.)
                            ax[i].title.set_fontsize(16)
                            ax[i].title.set_fontweight('bold')
                            is_title = False

                        name_area_line = names_area_line[k+i]
                        name_line = names_line[k+i]
                        ax[i].plot(fit_results['spectrum_index'],
                                 fit_results[name_area_line],
                                 'r.-')
                        ax[i].set_ylabel('Area %s [%s]'%(name_line.split('_')[0],name_line.split('_')[1]))

                    plt.xlabel('Spectrum index')
                    plt.show()

                # Plot the last pack with a number of elements < nb_rows
                else:

                    # Number of plots remaining
                    nb_plots_left = nb_plots%nb_rows

                    fig, ax = plt.subplots(figsize=(15,4*nb_plots_left), nrows=nb_plots_left)

                    if nb_plots_left == 1:
                        ax = [ax]

                    for i in range(nb_plots_left):
                        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%g'))
                        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

                        if is_title:
                            ax[i].set_title('AREAS OF FLUORESCENCE LINES', pad = 15, y=1.)
                            ax[i].title.set_fontsize(16)
                            ax[i].title.set_fontweight('bold')
                            is_title = False

                        name_area_line = names_area_line[k+i]
                        name_line = names_line[k+i]
                        ax[i].plot(fit_results['spectrum_index'],
                                 fit_results[name_area_line],
                                 'r.-')
                        ax[i].set_ylabel('Area %s [%s]'%(name_line.split('_')[0],name_line.split('_')[1]))

                    plt.xlabel('Spectrum index')
                    plt.show()


    def plot_fit_parameters_from_file(self, path_to_fit_params_save, path_to_fit_params_results):
        '''
        Plot the result parameters of a fit after loading it from file.
        We do not use the attributes of the current expt because this function should be callable even
        after the kernel has been restarted.

        Parameters
        ----------
        path_to_fit_params_save : str
            Relative path to the csv file with the saved fit parameters.
        path_to_fit_params_results : str
            Relative path to the csv file with the parameters resulting from the fit.
        '''
        # Import results and params
        fit_params = np.genfromtxt(path_to_fit_params_save, delimiter=';', names=True, dtype=None, encoding=None)
        fit_results = np.genfromtxt(path_to_fit_params_results, delimiter=';', names=True)

        # Construct the list of fitted parameters ('False' if empty)
        list_isfit_str = str(fit_params['list_isfit_str']).split(',')

        if list_isfit_str[0]!='False':
            # Total number of plots to print
            nb_plots = len(list_isfit_str)
        else:
            nb_plots = 0

        # Divide the plots in figures of nb_rows plots
        # To avoid overflowing the pdf with all plots on one page
        nb_rows = 5

        if nb_plots > 0:

            # Token to print the title only once
            is_title = True

            for k in range(0, nb_plots, nb_rows):

                # Print all the packs of nb_rows plots
                if (nb_plots-k)//nb_rows>0:

                    # Create a new figure for each pack of nb_rows plots
                    if k%nb_rows == 0:
                        fig, ax = plt.subplots(figsize=(15,4*nb_rows), nrows=nb_rows)

                        # ax is a float when there is only one row
                        if nb_plots == 1:
                            ax = [ax]

                    # Plot all the plots of the corresponding pack
                    for i in range(nb_rows):
                        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%g'))
                        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

                        if is_title:
                            ax[i].set_title('OTHER FIT PARAMETERS', pad = 15, y=1.)
                            ax[i].title.set_fontsize(16)
                            ax[i].title.set_fontweight('bold')
                            is_title = False

                        name_param = list_isfit_str[k+i]
                        ax[i].plot(fit_results['spectrum_index'],
                                 fit_results[name_param],
                                 'r.-')
                        ax[i].set_ylabel(name_param)

                    # Add the xlabel only on the last plot of the the pack
                    plt.xlabel('Spectrum index')
                    plt.show()

                # Plot the last pack with a number of elements < nb_rows
                else:

                    # Number of plots remaining
                    nb_plots_left = nb_plots%nb_rows

                    fig, ax = plt.subplots(figsize=(15,4*nb_plots_left), nrows=nb_plots_left)

                    if nb_plots_left == 1:
                        ax = [ax]

                    for i in range(nb_plots_left):
                        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%g'))
                        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

                        if is_title:
                            ax[i].set_title('OTHER FIT PARAMETERS', pad = 15, y=1.)
                            ax[i].title.set_fontsize(16)
                            ax[i].title.set_fontweight('bold')
                            is_title = False

                        name_param = list_isfit_str[k+i]
                        ax[i].plot(fit_results['spectrum_index'],
                                 fit_results[name_param],
                                 'r.-')
                        ax[i].set_ylabel(name_param)

                    plt.xlabel('Spectrum index')
                    plt.show()


    def plot_fit_positions_from_file(self, path_to_fit_params_results):
        '''
        Plot the result peaks positions of a fit after loading it from file.
        We do not use the attributes of the current expt because this function should be callable even
        after the kernel has been restarted.

        Parameters
        ----------
        path_to_fit_params_results : str
            Relative path to the csv file with the parameters resulting from the fit.
        '''
        # Import results
        fit_results = np.genfromtxt(path_to_fit_params_results, delimiter=';', names=True)

        # ['pos_peak_X_xx', 'pos_peak_Ar_KM3']
        names_pos_peak = [x for x in fit_results.dtype.names if x.startswith('pos_')]

        # ['X_xx', 'Ar_KM3']
        names_peak = [x.split('pos_peak_')[1] for x in names_pos_peak]

        # Total number of plots to print
        nb_plots = len(names_peak)

        # Divide the plots in figures of nb_rows plots
        # To avoid overflowing the pdf with all plots on one page
        nb_rows = 5

        if nb_plots > 0:

            # Token to print the title only once
            is_title = True

            for k in range(0, nb_plots, nb_rows):

                # Print all the packs of nb_rows plots
                if (nb_plots-k)//nb_rows>0:

                    # Create a new figure for each pack of nb_rows plots
                    if k%nb_rows == 0:
                        fig, ax = plt.subplots(figsize=(15,4*nb_rows), nrows=nb_rows)

                        # ax is a float when there is only one row
                        if nb_plots == 1:
                            ax = [ax]

                    # Plot all the plots of the corresponding pack
                    for i in range(nb_rows):
                        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%g'))
                        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

                        if is_title:
                            ax[i].set_title('POSITIONS OF FLUORESCENCE PEAKS', pad = 15, y=1.)
                            ax[i].title.set_fontsize(16)
                            ax[i].title.set_fontweight('bold')
                            is_title = False

                        name_pos_peak = names_pos_peak[k+i]
                        name_peak = names_peak[k+i]
                        ax[i].plot(fit_results['spectrum_index'],
                                 fit_results[name_pos_peak],
                                 'r.-')
                        ax[i].set_ylabel('Position %s [%s] (eV)'%(name_peak.split('_')[0],name_peak.split('_')[1]))

                    plt.xlabel('Spectrum index')
                    plt.show()

                # Plot the last pack with a number of elements < nb_rows
                else:

                    # Number of plots remaining
                    nb_plots_left = nb_plots%nb_rows

                    fig, ax = plt.subplots(figsize=(15,4*nb_plots_left), nrows=nb_plots_left)

                    if nb_plots_left == 1:
                        ax = [ax]

                    for i in range(nb_plots_left):
                        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%g'))
                        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

                        if is_title:
                            ax[i].set_title('POSITIONS OF FLUORESCENCE LINES', pad = 15, y=1.)
                            ax[i].title.set_fontsize(16)
                            ax[i].title.set_fontweight('bold')
                            is_title = False

                        name_pos_peak = names_pos_peak[k+i]
                        name_peak = names_peak[k+i]
                        ax[i].plot(fit_results['spectrum_index'],
                                 fit_results[name_pos_peak],
                                 'r.-')
                        ax[i].set_ylabel('Position %s [%s] (eV)'%(name_peak.split('_')[0],name_peak.split('_')[1]))

                    plt.xlabel('Spectrum index')
                    plt.show()

    def plot_all_fit_results(self, fit_index, path_to_result_folder):
        '''
        Function to call each individual plotting subfunctions when printing the report.
        We do not use the attributes of the current expt because this function should be callable even
        after the kernel has been restarted.

        Parameters
        ----------
        fit_index : int
            Index of the fitted spectrum to display.
        path_to_result_folder : str
            Relative path to the folder where the results are stored.
        '''

        # Reconstruct the different paths
        # e.g. path_to_result_folder = 'save/SIRIUS_Fluo_2020_02_16_02289/20211102_135304/'
        path_to_fit_folder = path_to_result_folder+'fit_curves/'
        path_to_fit_params_results = path_to_result_folder+'fit_results.csv'
        path_to_fit_params_save = path_to_result_folder+'fit_params.csv'

        # Plot the fit resulting curve in the tab
        self.plot_fit_curve_from_file(path_to_fit_folder, fit_index)

        # Plot the time series of areas in the tab
        self.plot_fit_areas_from_file(path_to_fit_params_results)

        # Plot the time series of peaks positions in the tab
        self.plot_fit_positions_from_file(path_to_fit_params_results)

        # Plot the time series of fit parameters in the tab
        self.plot_fit_parameters_from_file(path_to_fit_params_save,
                                           path_to_fit_params_results)

    def extract_avg_from_fit(self, wid):
        '''
        Get averages on fitted parameters and set the corresponding values in the widgets.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        if self.is_last_fit_a_preview:
            # Get the results directly from the fit output
            fit_params_results = self.result.params

            str_to_add = 'Widgets updated with averages from last fit preview.'

        else:
            # Get the results from file
            fit_params_results = np.genfromtxt(self.path_to_fit_params_results, delimiter=';', names=True)

            str_to_add = 'Widgets updated with averages from last series of fits.'


        if 'gain' in self.list_isfit:
            wid.gain.value = np.round(np.nanmean(fit_params_results['gain']),4)

        if 'eV0' in self.list_isfit:
            wid.eV0.value = np.round(np.nanmean(fit_params_results['eV0']),4)

        if 'ct' in self.list_isfit:
            wid.ct.value = np.round(np.nanmean(fit_params_results['ct']),4)

        if 'sl' in self.list_isfit:
            wid.sl.value = np.round(np.nanmean(fit_params_results['sl']),4)

        if 'noise' in self.list_isfit:
            wid.noise.value = np.round(np.nanmean(fit_params_results['noise']),4)

        if 'SF0' in self.list_isfit:
            wid.SF0.value = np.round(np.nanmean(fit_params_results['SF0']),4)

        if 'TF0' in self.list_isfit:
            wid.TF0.value = np.round(np.nanmean(fit_params_results['TF0']),4)

        if 'TW0' in self.list_isfit:
            wid.TW0.value = np.round(np.nanmean(fit_params_results['TW0']),4)

        if 'broad_factor' in self.list_isfit:
            wid.broad_factor.value = np.round(np.nanmean(fit_params_results['broad_factor']),4)

        if 'lowTF0' in self.list_isfit:
            wid.lowTF0.value = np.round(np.nanmean(fit_params_results['lowTF0']),4)

        if 'highTF0' in self.list_isfit:
            wid.highTF0.value = np.round(np.nanmean(fit_params_results['highTF0']),4)

        if 'lowTW0' in self.list_isfit:
            wid.lowTW0.value = np.round(np.nanmean(fit_params_results['lowTW0']),4)

        if 'highTW0' in self.list_isfit:
            wid.highTW0.value = np.round(np.nanmean(fit_params_results['highTW0']),4)

        # Update logs
        if self.logs_lvl>=0:
            self.add_str_to_logs(wid, str_to_add)

    def add_str_to_logs(self, wid, str_to_add):
        '''
        Add a string to the logs and update the logs window.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        str_to_add : str
            String to add to the logs.
        '''

        # Add the current date
        date_to_add = _BOLD+datetime.now().strftime('%d/%m/%Y, %H:%M:%S')+_RESET

        # We update in reversed order, there is no way to scroll automatically
        # to the bottom of the output widget
        self.logs = date_to_add+ '\n'+ str_to_add + '\n' + self.logs
        with wid.out_logs:
            wid.out_logs.clear_output()
            print(self.logs)



    def generate_report(self, wid):
        '''
        Generate the text for the report.

        Parameters
        ----------
        wid : object myWidgets
            Object from the class myWidgets.
        '''
        self.report = ['# '+self.name]
        self.report.append('break')
        self.report.append('## Session ID: '+self.session_id)

        self.report.append('Fit results for file ```%s```'%self.filename)
        self.report.append('Parameters and results saved in:  \n```%s```'%\
                           (self.save_dir+self.name+'/'+self.session_id+'/'))
        self.report.append('break')

        self.report.append('**Extraction parameters**')
        self.report.append('Spectrum interval = [%g, %g]'%(self.first_spectrum,self.last_spectrum))
        self.report.append('Channel interval = [%g, %g]'%(self.first_channel,self.last_channel))
        self.report.append('SDD elements used = %s'%(str(self.SDD_elems_chosen_int)))
        self.report.append('break')

        self.report.append('**Fit parameters**')
        self.report.append('List of fitted parameters: %s'%str(self.list_isfit))
        self.report.append('beam energy = %g'%self.beam_energy +'; min. strength = %g'%self.min_strength)
        self.report.append('**Params for conversion to eVs**')
        self.report.append('gain = %g'%self.gain +'; eV0 = %g'%self.eV0)
        self.report.append('**Params for linear background**')
        self.report.append('slope = %g'%self.sl+'; constant = %g'%self.ct)
        if self.is_bckg_subset:
            self.report.append('Fit the background on the subset [%g eV, %g eV]'%(self.bckg_eVs_inf, self.bckg_eVs_sup))
        self.report.append('**Params for elastic peaks**')
        self.report.append('noise = %g'%self.noise
                           +'; tail fraction (low energy side) = %g'%self.TF0)
        self.report.append('tail width (low energy side) = %g'%self.TW0+'; shelf fraction = %g'%self.SF0)
        self.report.append('**Params for Compton peaks**')
        self.report.append('broadening factor = %g'%self.broad_factor+'; tail fraction (low energy side) = %g'%self.lowTF0\
                          +'; tail fraction (high energy side) = %g'%self.highTF0)
        self.report.append('tail width (low energy side) = %g'%self.lowTW0+'; tail width (high energy side) = %g'%self.highTW0)

        # Update logs
        if self.logs_lvl>=0:
            str_to_add = 'Report generated.'
            self.add_str_to_logs(wid, str_to_add)
