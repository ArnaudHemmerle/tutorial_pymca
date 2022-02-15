'''
Module for defining the widgets in the top bar.
'''
import ipywidgets as widgets
from IPython.display import display
from lib.frontend import notebook as nb

# Define colors for prints
_RED='\x1b[31;01m'
_RESET="\x1b[0m"

def create_bar(expt, wid):
    '''
    Create the widgets for the top bar.

    Parameters
    ----------
    expt : object
        object from the class Experiment.
    wid : object myWidgets
        object from the class myWidgets.
    '''
    def on_button_export_clicked(b):
        '''
        Export the notebook to PDF.
        '''
        # Update logs
        if expt.logs_lvl>=0:
            str_to_add = _RED+'Export to pdf in progress...'+_RESET
            expt.add_str_to_logs(wid, str_to_add)

        export_done = expt.export_nb_to_pdf(wid)

        if export_done:
            str_to_add = 'Notebook exported to %s.pdf'%expt.notebook_name.split('.')[0]
        else:
            str_to_add = 'There was something wrong with the export to pdf.'+'\n'\
                        +'Did you rename the Notebook? If yes:'+'\n'\
                        +'1) Change the value of expt.notebook_name in the first cell (top of the Notebook).'+'\n'\
                        +'2) Re-execute the first cell.'+'\n'\
                        +'3) Try to export the pdf again in the last cell (bottom of the Notebook).'

        if expt.logs_lvl>=0:
            expt.add_str_to_logs(wid, str_to_add)

    def on_button_save_all_params_as_default_clicked(b):
        '''
        Save all parameters in the corresponding default csv file.
        '''
        ##########################
        # Extraction parameters

        # Save
        expt.save_extract_params_in_file(wid, expt.path_to_extract_params_default)

        ##########################
        # Peaks parameters
        if expt.is_sheet_loaded:
            # Validate the sheet
            expt.validate_sheet(wid)

            # Save
            expt.save_peaks_params_in_file(wid, expt.path_to_peaks_params_default)

        else:
            # Update logs
            if expt.logs_lvl>=0:
                str_to_add = _RED+'Peaks are not loaded yet.'+_RESET
                expt.add_str_to_logs(wid, str_to_add)

        ##########################
        # Fit parameters

        # Set the fit params from the values in the widgets
        expt.set_fit_params_from_widgets(wid)

        # Save
        expt.save_fit_params_in_file(wid, expt.path_to_fit_params_default)

    def on_button_analyze_next_scan_clicked(b):
        '''
        Close the current main widget and open a new one at the bottom of the notebook.
        '''

        cell = \
         'expt = Experiment(notebook_name, save_dir, files_dir, logs_lvl)\n'\
        +'wid = myWidgets()\n'\
        +'FE.start_session(expt, wid)'
        nb.create_cell(code=cell, position ='at_bottom', celltype='code', is_print=False)

        nb.delete_current_cell()

    def on_button_save_and_close_clicked(b):
        '''
        Save and delete the current widget (so that the notebook is clean next time it is opened).
        '''

        nb.delete_current_cell()

        nb.save()

    def on_button_extract_avg_clicked(b):
        '''
        Get averages on fitted parameters and set the corresponding values in the widgets.
        '''
        if (expt.is_fit_done or expt.is_last_fit_a_preview):
            # Extract avg and update widgets
            expt.extract_avg_from_fit(wid)

        else:
            # Update logs
            if expt.logs_lvl>=0:
                str_to_add = _RED+'Do a fit first.'+_RESET
                expt.add_str_to_logs(wid, str_to_add)

    ####################################
    # Define widgets
    ####################################

    ####################################
    # Level 1 STATIC
    # General options
    ####################################

    wid.button_save_all_params_as_default.on_click(on_button_save_all_params_as_default_clicked)
    wid.button_export.on_click(on_button_export_clicked)
    wid.button_analyze_next_scan.on_click(on_button_analyze_next_scan_clicked)
    wid.button_save_and_close.on_click(on_button_save_and_close_clicked)
    wid.button_extract_avg.on_click(on_button_extract_avg_clicked)

    out_level_1 = widgets.Output()
    wid_level_1 = widgets.VBox([
                        widgets.HBox([
                          wid.button_export,
                          wid.button_save_all_params_as_default,
                          wid.button_extract_avg
                                     ]),
                        widgets.HBox([
                          wid.button_analyze_next_scan,
                          wid.button_save_and_close
                                     ])
                              ])

    with out_level_1:
        display(wid_level_1)

    display(widgets.VBox([out_level_1]))
