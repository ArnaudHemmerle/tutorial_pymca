'''
Module for fitting a series of spectrum.
'''
import os
import shutil
import ipywidgets as widgets
from IPython.display import display

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
    def on_button_import_fit_params_file_clicked(b):
        '''
        Import csv file and update values of expt and widgets.
        '''
        if (expt.is_extract_done and expt.is_sheet_loaded):
            # Choose between default params or a previous set of params
            if wid.select_fit_params_file.value == 'Default fit parameters':
                expt.path_to_fit_params = expt.path_to_fit_params_default
            else:
                expt.path_to_fit_params = expt.save_dir+wid.select_fit_params_file.value+'/fit_params.csv'

            # Import the csv file and update the sheet
            expt.set_fit_params_from_file(wid, expt.path_to_fit_params)
            wid.set_fit_params_from_expt(expt)

        else:

            # Update logs
            if expt.logs_lvl>=0:
                str_to_add = _RED+'Extract a scan and define the peaks first.'+_RESET
                expt.add_str_to_logs(wid, str_to_add)

    def on_button_save_fit_params_as_default_clicked(b):
        '''
        Save the current fit parameters in the default csv file.
        '''
        # Get params from widgets
        expt.set_fit_params_from_widgets(wid)

        # Save
        expt.save_fit_params_in_file(wid, expt.path_to_fit_params_default)

    def on_button_run_single_fit_clicked(b):
        '''
        Run a fit on the current spectrum.
        '''
        if (expt.is_extract_done and expt.is_sheet_loaded):

            # Update logs
            if expt.logs_lvl>=0:
                str_to_add = _RED+'Fit in progress...'+_RESET
                expt.add_str_to_logs(wid, str_to_add)

            # Validate the sheet
            expt.validate_sheet(wid)

            # Set the fit params from the values in the widgets
            expt.set_fit_params_from_widgets(wid)

            # Get the spectrum to fit from the widget
            expt.spectrum_to_fit = expt.spectrums[wid.fit_index.value-expt.first_spectrum]

            # Do the fit
            expt.run_single_fit(wid)

            with out_level_5:
                out_level_5.clear_output()

                # Plot the result of the fit
                expt.plot_single_fit(wid, title='Fit of spectrum %g'%(wid.fit_index.value))

                if expt.is_spectrum_empty:
                    print('Spectrum is empty.')
                elif expt.is_fit_fail:
                    print('Fit failed to converge.')
                else:
                    display(expt.result)

            # Update logs
            if expt.logs_lvl>=0:
                if expt.is_spectrum_empty:
                    str_to_add = 'Spectrum %s is empty.'%str(wid.fit_index.value)
                elif expt.is_fit_fail:
                    str_to_add = 'Fit of spectrum %s failed to converge.'%str(wid.fit_index.value)
                else:
                    str_to_add = 'Fit done.'
                expt.add_str_to_logs(wid, str_to_add)

        else:
            # Update logs
            if expt.logs_lvl>=0:
                if not expt.is_extract_done:
                    str_to_add = _RED+'Extract a scan first.'+_RESET
                if not expt.is_sheet_loaded:
                    str_to_add = _RED+'Load peaks first.'+_RESET
                expt.add_str_to_logs(wid, str_to_add)

    def run_series_fit():
        '''
        Procedure for preparing the series of fits, doing it, and saving the results.
        '''

        # Update logs
        if expt.logs_lvl>=0:
            str_to_add = _RED+'Fits in progress...'+_RESET
            expt.add_str_to_logs(wid, str_to_add)

        # Create the csv file with the results
        expt.create_fit_params_results(wid, expt.path_to_fit_params_results)

        # Empty the folder fit_curves
        if os.path.exists(expt.path_to_fit_folder):
            shutil.rmtree(expt.path_to_fit_folder)
        os.mkdir(expt.path_to_fit_folder)

        # Get the current logs lvl
        old_logs_lvl = expt.logs_lvl

        for spectrum_index in range(expt.first_spectrum, expt.last_spectrum+1):

            # Get the spectrum to fit
            expt.spectrum_to_fit = expt.spectrums[spectrum_index-expt.first_spectrum]

            # Do the fit
            expt.run_single_fit(wid)

            # Save the fit curves in csv files, using the timestamps to help user matching with the fit_results.csv file
            expt.current_sensorsRelTimestamps = expt.sensorsRelTimestamps[spectrum_index]
            expt.path_to_fit_curve_results = expt.path_to_fit_folder+'spectrum_'+str(spectrum_index)+'.csv'
            expt.save_fit_curves_results_in_file(wid, expt.path_to_fit_curve_results)

            # Save the fit results
            expt.save_fit_params_results_in_file(wid, expt.path_to_fit_params_results, spectrum_index)

            # Save the log of the fit results
            expt.save_fit_logs_in_file(wid, expt.path_to_fit_log_results, spectrum_index)

            # Set the token
            expt.is_fit_done = True
            expt.is_last_fit_a_preview = False

            with out_level_4:
                out_level_4.clear_output()

            with out_level_5:
                # Plot the result of the fit

                if ((spectrum_index-expt.first_spectrum)%wid.plot_frequency.value==0
                    or spectrum_index==expt.last_spectrum):
                    expt.plot_single_fit(wid, title='Fit of spectrum %g/%g'%(spectrum_index,expt.last_spectrum))
                    if expt.is_spectrum_empty:
                        print('Spectrum is empty.')
                    elif expt.is_fit_fail:
                        print('Fit failed to converge.')

            # Update logs
            if expt.logs_lvl>=0:
                if expt.is_spectrum_empty:
                    str_to_add = 'Spectrum %s is empty.'%str(spectrum_index)
                    expt.add_str_to_logs(wid, str_to_add)
                if expt.is_fit_fail:
                    str_to_add = 'Fit of spectrum %s failed to converge.'%str(spectrum_index)
                    expt.add_str_to_logs(wid, str_to_add)


            # Deactivate the logs of subfunction for the fit routine, if they are not intensive
            if expt.logs_lvl < 1:
                expt.logs_lvl = -1

        expt.logs_lvl = old_logs_lvl

        # Update logs
        if expt.logs_lvl>=0:
            str_to_add = 'Fits done.'
            expt.add_str_to_logs(wid, str_to_add)

    def on_button_replace_all_files_and_fit_clicked(b):
        '''
        Replace the files with the current params and clean the previous fit results.
        Then start the fit.
        '''
        # Get params from widgets and save in files
        expt.get_and_save_all_params_in_files(wid)

        # Empty and recreate the folder fit_curves
        if os.path.exists(expt.path_to_fit_folder):
            shutil.rmtree(expt.path_to_fit_folder)
        os.mkdir(expt.path_to_fit_folder)

        # Update logs
        if expt.logs_lvl>=0:
            str_to_add = 'Fit curves saved in the folder:\n%s'%expt.path_to_fit_folder
            expt.add_str_to_logs(wid, str_to_add)

        # Remove previous fit results
        if os.path.exists(expt.path_to_fit_params_results):
            os.remove(expt.path_to_fit_params_results)

        # Remove the fit logs
        if os.path.exists(expt.path_to_fit_log_results):
            os.remove(expt.path_to_fit_log_results)

        # Start the fit
        run_series_fit()

    def on_button_create_all_files_and_fit_clicked(b):
        '''
        Create a new session and save all params in csv files.
        Then start the fit.
        '''
        # Reset the session id and update the widget
        expt.set_session_id(wid)
        wid.set_session_id(expt)

        # Set the path to the different saving files
        expt.set_paths(wid)

        # Get params from widgets and save in files
        expt.get_and_save_all_params_in_files(wid)

        # Update list of files
        wid.select_fit_params_file.options = expt.list_fit_params_files

        # Start the fit
        run_series_fit()


    def on_button_run_all_fit_clicked(b):
        '''
        Run a fit for each spectrum extracted, and save the results in files.
        '''
        # Set the path to the different saving files
        expt.set_paths(wid)

        if (expt.is_extract_done and expt.is_sheet_loaded):

            # Set the fit params from the values in the widgets
            expt.set_fit_params_from_widgets(wid)

            ###############################################
            # Check if any params file exist
            # If yes ask the user if they should be all erased or if new ones should be created
            # If no create the folder and the csv file

            if (os.path.exists(expt.path_to_extract_params_save) or os.path.exists(expt.path_to_peaks_params_save)\
                or os.path.exists(expt.path_to_fit_params_save) or os.path.exists(expt.path_to_fit_folder)\
                or os.path.exists(expt.path_to_fit_params_results)):

                wid.button_replace_all_files_and_fit.on_click(on_button_replace_all_files_and_fit_clicked)
                wid.button_create_all_files_and_fit.on_click(on_button_create_all_files_and_fit_clicked)

                with out_level_4:
                    out_level_4.clear_output()
                    print(_RED+"Files already exist."+_RESET)
                    display(widgets.HBox([wid.button_replace_all_files_and_fit,
                                          wid.button_create_all_files_and_fit]))


            else:
                # Get params from widgets and save in files
                expt.get_and_save_all_params_in_files(wid)

                with out_level_4:
                    out_level_4.clear_output()

                # Start the fit
                run_series_fit()


        else:
            # Update logs
            if expt.logs_lvl>=0:
                if not expt.is_extract_done:
                    str_to_add = _RED+'Extract a scan first.'+_RESET
                    expt.add_str_to_logs(wid, str_to_add)
                if not expt.is_sheet_loaded:
                    str_to_add = _RED+'Load peaks first.'+_RESET
                    expt.add_str_to_logs(wid, str_to_add)

    ####################################
    # Define the widgets
    ####################################

    # Widget for selecting the file
    # Create the list of files available and display them
    expt.set_list_fit_params_files(wid)
    wid.select_fit_params_file.options = expt.list_fit_params_files

    wid.button_import_fit_params_file.on_click(on_button_import_fit_params_file_clicked)
    wid.button_save_fit_params_as_default.on_click(on_button_save_fit_params_as_default_clicked)
    wid.button_run_single_fit.on_click(on_button_run_single_fit_clicked)
    wid.button_run_all_fit.on_click(on_button_run_all_fit_clicked)


    #####################################
    # Level 1 STATIC
    # Load file and general options
    #####################################

    out_level_1 = widgets.Output()

    wid_level_1 = widgets.VBox([
                    widgets.HBox([
                        wid.select_fit_params_file,
                        widgets.VBox([wid.button_import_fit_params_file,
                                      wid.button_save_fit_params_as_default])
                                ]),
                             ])

    with out_level_1:
        display(wid_level_1)


    ######################################################
    # Level 2 DYNAMIC
    # Console of file selection
    ######################################################

    out_level_2 = widgets.Output()


    ##########################################
    # Level 3 STATIC
    # Fit parameters
    ##########################################

    out_level_3 = widgets.Output()
    wid_level_3 = widgets.VBox([
                    widgets.VBox([
                          widgets.VBox([widgets.Label(r'\(\textbf {Params for conversion to eVs.}\)'),
                                        widgets.Label('eVs = gain*channels + eV0')]),
                          widgets.HBox([wid.gain, wid.is_gain, wid.eV0, wid.is_eV0]),
                          widgets.VBox([widgets.Label(' '),
                                        widgets.Label(r'\(\textbf {Params for linear background.}\)')]),
                          widgets.VBox([widgets.HBox([wid.sl, wid.is_sl, wid.ct, wid.is_ct]),
                                        widgets.HBox([wid.is_bckg_subset, wid.bckg_eVs_inf, wid.bckg_eVs_sup])
                                        ]),
                          widgets.VBox([widgets.Label(' '),
                                        widgets.Label(r'\(\textbf {Params for elastic peaks.}\)')]),
                          widgets.VBox([widgets.HBox([wid.noise, wid.is_noise, wid.SF0, wid.is_SF0]),
                                        widgets.HBox([ wid.TW0, wid.is_TW0, wid.TF0, wid.is_TF0])
                                        ]),
                          widgets.VBox([widgets.Label(' '),
                                        widgets.Label(r'\(\textbf {Params for Compton peaks.}\)')]),
                          widgets.VBox([widgets.HBox([wid.lowTF0, wid.is_lowTF0, wid.highTF0, wid.is_highTF0]),
                                        widgets.HBox([wid.lowTW0, wid.is_lowTW0, wid.highTW0, wid.is_highTW0]),
                                        widgets.HBox([wid.broad_factor, wid.is_broad_factor])
                                        ])
                              ],layout={'border': '1px solid black', 'margin': '0px 0px 10px 0px'}),
                    widgets.HBox([wid.show_gaussian, wid.show_shelf, wid.show_tail, wid.show_compton, wid.show_bckg]),
                    widgets.HBox([wid.button_run_single_fit, wid.fit_index]),
                    widgets.HBox([wid.button_run_all_fit, wid.plot_frequency])
                                 ])

    with out_level_3:
        display(wid_level_3)


    ##########################################
    # Level 4 DYNAMIC
    # Console save
    ##########################################

    out_level_4 = widgets.Output()


    ##########################################
    # Level 5 DYNAMIC
    # Output fit
    ##########################################

    out_level_5 = widgets.Output()

    display(widgets.VBox([out_level_1, out_level_2, out_level_3, out_level_4, out_level_5]))
