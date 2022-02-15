'''
Module for extraction of the spectrum from the scan.
'''
import os
import numpy as np
import ipywidgets as widgets
from IPython.display import display

from lib.extraction import PyNexus as PN
from lib.extraction import XRF

# Define colors for prints
_RED='\x1b[31;01m'
_RESET="\x1b[0m"

def create_tab(expt, wid):
    '''
    Create the widgets for the corresponding tab.

    Parameters
    ----------
    expt : object
        object from the class Experiment
    wid : object myWidgets
        object from the class myWidgets
    '''
    def finteractive_select_file(filename):
        '''
        Called by wid.interactive_select_file, the interactive widget for scan selection.

        Parameters
        ----------
        filename : str
            name of the file with extension (i.e. 'SIRIUS_Fluo_2020_07_07_0070.nxs' or 'SIRIUS_Fluo_2020_07_07_0070.nxs')
        '''
        if expt.list_files_filenames == []:

            str_to_add = _RED+'There is no nexus or dat file in the files directory.'+'\n'\
                             +'Files directory: %s'%expt.files_dir+_RESET

            with out_level_3:
                out_level_3.clear_output()
                print(str_to_add)

            # Update logs
            if expt.logs_lvl>=0:
                expt.add_str_to_logs(wid, str_to_add)

        else:

            expt.filename = filename
            expt.path_to_file = expt.files_dir+expt.filename
            expt.name = expt.filename[:-4]

            # Quick extraction of scan info
            expt.set_scan_info(wid)

            with out_level_3:
                out_level_3.clear_output()

            str_to_add = 'There are %g data points in the scan.'%(expt.nb_allspectrums)

            if expt.SDD_elems_available == []:
                str_to_add += '\n'+'No SDD element available in the file. Not a XRF scan?'

            else:
                str_to_add += '\n'+'List of available SDD elements: %s'%str(['%s'%s for s in expt.SDD_elems_available])

            with out_level_3:
                print(str_to_add)

            # Update logs
            if expt.logs_lvl>=0:
                expt.add_str_to_logs(wid, str_to_add)

            #####################################################################
            # Import the latest parameters if they exist, if not the default ones

            expt.path_to_extract_params = expt.path_to_extract_params_default

            # Import previous params if they are found
            if os.path.exists(expt.save_dir+expt.name):
                list_save_dir = sorted(os.listdir(expt.save_dir+expt.name))

                if list_save_dir != []:

                    path_to_previous_extract_params = \
                    expt.save_dir+expt.name+'/'+list_save_dir[-1]+'/'+'extract_params.csv'

                    if os.path.exists(path_to_previous_extract_params):
                        expt.path_to_extract_params = path_to_previous_extract_params

            # Import the csv file and update expt and widgets
            expt.set_extract_params_from_file(wid, expt.path_to_extract_params)
            wid.set_extract_params_from_expt(expt)


    def on_button_select_extract_params_to_load_clicked(b):
        '''
        Display a widget for choosing the file and load it on click.
        '''
        def on_button_import_extract_params_file(b):
            '''
            Import csv file and update values of expt and widgets.
            '''
            # Choose between default params or a previous set of params
            if wid.select_extract_params_file.value == 'Default extraction parameters':
                expt.path_to_extract_params = expt.path_to_extract_params_default
            else:
                expt.path_to_extract_params = expt.save_dir+wid.select_extract_params_file.value+'/extract_params.csv'

            # Import the csv file and update expt and widgets
            expt.set_extract_params_from_file(wid, expt.path_to_extract_params)
            wid.set_extract_params_from_expt(expt)


        #Set the list of csv files containing extraction parameters
        expt.set_list_extract_params_files(wid)

        # Update default value on the widget for selecting the csv file
        wid.select_extract_params_file.options = expt.list_extract_params_files

        # Define the widget
        wid.button_import_extract_params_file.on_click(on_button_import_extract_params_file)

        with out_level_3:
            out_level_3.clear_output()
            display(widgets.HBox([wid.select_extract_params_file, wid.button_import_extract_params_file]))

        with out_level_4:
            out_level_4.clear_output()


    def on_button_save_extract_params_as_default_clicked(b):
        '''
        Save the current extraction parameters in the default csv file.
        '''
        # Get params from widgets
        expt.set_extract_params_from_widgets(wid)

        # Save
        expt.save_extract_params_in_file(wid, expt.path_to_extract_params_default)

    def on_button_extract_clicked(b):
        '''
        Extract the selected scan.
        '''
        # Update the default value in the peak selection widget
        wid.set_default_peak_value(expt)

        # Update the default value in the fit selection widget
        wid.set_default_fit_value(expt)

        # Get params from widgets
        expt.set_extract_params_from_widgets(wid)

        #######################
        # Create list of SDD elements from the booleans
        SDD_elems_chosen_bool = np.array([expt.is_SDD_elem_0,expt.is_SDD_elem_1,
                                          expt.is_SDD_elem_2,expt.is_SDD_elem_3,
                                          expt.is_SDD_elem_4])

        # Convert the list of booleans to a list of integers
        tmp_elems = SDD_elems_chosen_bool*[10,11,12,13,14]
        expt.SDD_elems_chosen_int = [i-10 for i in tmp_elems if i>0]

        # Test if all the selected SDD elems are present in the nexus file
        test_SDD_elems_chosen = all(elem in expt.SDD_elems_available \
                                    for elem in [str(x) for x in expt.SDD_elems_chosen_int])

        if test_SDD_elems_chosen:

            if '.nxs' in expt.filename:

                # Extract the XRF over the whole range of channels and non-zero spectrums
                _ , _ , expt.all_spectrums, _ , expt.last_non_zero_spectrum = \
                XRF.Extract(nxs_filename = expt.filename, recording_dir = expt.files_dir,
                            list_elems = expt.SDD_elems_chosen_int,
                            first_channel = 0, last_channel = 2048,
                            gain = 1., eV0 = 0.,
                            fast = True, show_data_stamps = False, verbose = False)


                # Extract stamps and data from the sensors recorded in the scan
                nexus = PN.PyNexusFile(expt.path_to_file, fast=True)
                expt.stamps0D, expt.data0D= nexus.extractData(which='0D')
                for i, stamp0D in enumerate(expt.stamps0D):
                    if (stamp0D[1] is None \
                        and stamp0D[0] in ['sensorsRelTimestamps', 'sensors_rel_timestamps']):
                        # Time stamps used for later stitching of the fit results to the sensors data
                        expt.sensorsRelTimestamps = expt.data0D[i]

            if '.dat' in expt.filename:

                # Extract the XRF over the whole range of channels and non-zero spectrums
                _ , _ , expt.all_spectrums, _ , expt.last_non_zero_spectrum = \
                XRF.Extract_mat(filename = expt.filename[:-4], recording_dir = expt.files_dir,
                                list_elems = expt.SDD_elems_chosen_int,
                                first_channel = 0, last_channel = 2048,
                                gain = 1., eV0 = 0.,
                                show_data_stamps = False, verbose = False)


                path_to_dat = expt.files_dir+expt.filename[:-4]+'.dat'
                dat_extracted = np.genfromtxt(path_to_dat, names=True)
                expt.stamps0D = [(name, name, None) for name in dat_extracted.dtype.names]
                expt.data0D = np.genfromtxt(path_to_dat, skip_header=1).transpose()
                for name in dat_extracted.dtype.names:
                    if 'sensorsRelTimestamps' in name:
                        expt.sensorsRelTimestamps = dat_extracted['sensorsRelTimestamps']
                    if 'sensors_rel_timestamps' in name:
                        expt.sensorsRelTimestamps = dat_extracted['sensors_rel_timestamps']

            # Subset of channels and spectrums defined by user
            expt.channels = np.arange(expt.first_channel, expt.last_channel+1)
            expt.spectrums = expt.all_spectrums[expt.first_spectrum:expt.last_spectrum+1,
                                                expt.first_channel:expt.last_channel+1]
            expt.nb_spectrums = len(expt.spectrums)



            # Set the maximum value of the widget for selecting a spectrum in Fit tab
            # Need an intermediate step to avoid max<min

            wid.fit_index.min = 0
            wid.fit_index.value = 0
            wid.fit_index.max = 0

            wid.fit_index.max = expt.last_spectrum
            wid.fit_index.min = expt.first_spectrum
            wid.fit_index.value = expt.first_spectrum


            # Sum on all the spectrums, used for plotting the peaks
            expt.spectrums_sum = np.sum(expt.spectrums, axis=0)

            # Give the info that the extraction was done
            expt.is_extract_done = True

            # Check if the range of selected spectrums is empty
            if expt.spectrums.sum()<10.:
                with out_level_3:
                    out_level_3.clear_output()
                    print('Extraction of:\n%s'%expt.path_to_file)
                    print(_RED+'The selected range of spectrums appears to be empty!'+_RESET)
                    print('Select another range of spectrums.')

                # Update logs
                if expt.logs_lvl>=0:
                    str_to_add = 'Extraction of:\n%s'%expt.path_to_file+'\n'\
                                +_RED+'The selected range of spectrums appears to be empty!'+_RESET+'\n'\
                                +'Select another range of spectrums.'
                    expt.add_str_to_logs(wid, str_to_add)

                with out_level_4:
                    out_level_4.clear_output()
                    # Plot the result of the extraction (but only the first part)
                    expt.plot_extraction(wid, is_range_spectrums_empty=True)

            else:
                str_to_add = 'Extraction of:\n%s'%expt.path_to_file+'\n'\
                             +'Spectrums empty after data point %g.'%expt.last_non_zero_spectrum+'\n'\
                             +_RED+'Extraction in progress...'+_RESET

                with out_level_3:
                    out_level_3.clear_output()
                    print(str_to_add)

                # Update logs
                if expt.logs_lvl>=0:
                    expt.add_str_to_logs(wid, str_to_add)


                with out_level_4:
                    out_level_4.clear_output()

                    # Plot the result of the extraction
                    expt.plot_extraction(wid)

                    # Update logs
                    if expt.logs_lvl>=0:
                        str_to_add = 'Extraction done.'
                        expt.add_str_to_logs(wid, str_to_add)

        else:

            str_to_add = _RED+'Selected SDD elements not present in the nexus file.'+_RESET
            with out_level_3:
                out_level_3.clear_output()
                print(str_to_add)

            with out_level_4:
                out_level_4.clear_output()

            # Update logs
            if expt.logs_lvl>=0:
                expt.add_str_to_logs(wid, str_to_add)



    def on_button_refresh_clicked(b):
        '''
        Update the list of files in the corresponding folder.
        '''

        # Update logs
        if expt.logs_lvl>=0:
            str_to_add = 'Refresh the list of nexus files.'
            expt.add_str_to_logs(wid, str_to_add)

        expt.set_list_files(wid)
        wid.set_list_files(expt)


    ####################################
    # Define widgets
    ####################################


    # It is more convenient to define interactive widgets here than in myWidgets
    wid.interactive_select_file = widgets.interactive(finteractive_select_file,
                                                     filename = wid.select_file)

    # Link widgets to functions
    wid.button_save_extract_params_as_default.on_click(on_button_save_extract_params_as_default_clicked)
    wid.button_select_extract_params_to_load.on_click(on_button_select_extract_params_to_load_clicked)
    wid.button_extract.on_click(on_button_extract_clicked)
    wid.button_refresh.on_click(on_button_refresh_clicked)


    ##########################################################################
    # Level 3 DYNAMIC
    # Console
    # Needs to be defined before level 1 (because of the interactive widget)
    ##########################################################################
    out_level_3 = widgets.Output()


    ####################################
    # Level 1 STATIC
    # Scan selection and general options
    ####################################

    out_level_1 = widgets.Output()
    wid_level_1 = widgets.HBox([
                          wid.interactive_select_file,
                          widgets.VBox([wid.button_extract,
                                        wid.button_select_extract_params_to_load,
                                        wid.button_save_extract_params_as_default,
                                        wid.button_refresh
                                        ])
                              ])

    with out_level_1:
        display(wid_level_1)

    ##########################################
    # Level 2 STATIC
    # Extraction parameters
    ##########################################

    out_level_2 = widgets.Output()

    wid_level_2 = widgets.VBox([
                          widgets.HBox([
                                  widgets.Label('Include SDD elements:'),
                                  wid.is_SDD_elem_0, wid.is_SDD_elem_1, wid.is_SDD_elem_2,
                                  wid.is_SDD_elem_3, wid.is_SDD_elem_4
                                      ]),
                          widgets.HBox([
                                  wid.first_channel, wid.last_channel, wid.first_spectrum, wid.last_spectrum
                                      ])
                              ],layout={'border': '1px solid black'})

    with out_level_2:
        display(wid_level_2)

    ######################################################
    # Level 4 DYNAMIC
    # Plots
    ######################################################
    out_level_4 = widgets.Output()


    display(widgets.VBox([out_level_1,out_level_2,
                          out_level_3,out_level_4]))
