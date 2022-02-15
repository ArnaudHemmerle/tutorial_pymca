'''
Module for printing the report in the notebook.
'''
import ipywidgets as widgets

from IPython.display import display, Markdown
from lib.frontend import notebook as nb

# Define colors for prints
_RED='\x1b[31;01m'
_RESET="\x1b[0m"

def create_tab(expt, wid):
    '''
    Create the widgets for the corresponding tab.

    Parameters
    ----------
    expt : object
        object from the class Experiment.
    wid : object myWidgets
        object from the class myWidgets.
    '''
    def on_button_preview_report_clicked(b):
        '''
        Preview the report.
        '''

        if expt.is_fit_done:

            with out_level_2:
                out_level_2.clear_output()

                # Generate the report, split it, and display it in the report tab
                expt.generate_report(wid)
                split_report='  \n'.join(expt.report).split('break')
                display(Markdown('  \n'.join(split_report)))

                # Plot the fit resulting curve in the tab
                expt.plot_fit_curve_from_file(expt.path_to_fit_folder, wid.fit_index.value)

                # Plot the time series of areas in the tab
                expt.plot_fit_areas_from_file(expt.path_to_fit_params_results)

                # Plot the time series of peaks positions in the tab
                expt.plot_fit_positions_from_file(expt.path_to_fit_params_results)

                # Plot the time series of fit parameters in the tab
                expt.plot_fit_parameters_from_file(expt.path_to_fit_params_save,
                                                   expt.path_to_fit_params_results)


                str_to_add = 'Preview of the notebook cells generated.'

        else:
            with out_level_2:
                out_level_2.clear_output()

                str_to_add = _RED+'Do a fit first.'+_RESET

                print(str_to_add)

        # Update logs
        if expt.logs_lvl>=0:
            expt.add_str_to_logs(wid, str_to_add)

    def on_button_print_report_clicked(b):
        '''
        Print the report in the notebook.
        The cells have to be generated in reversed order.
        '''
        # Generate the report and update the preview
        on_button_preview_report_clicked(b)

        if expt.is_fit_done:

            path_to_result_folder = '/'.join(expt.path_to_fit_params_save.split('/')[:-1])+'/'

            # Add and execute a code cell for the plotting all the fit results
            command = \
            'fit_index='+str(wid.fit_index.value)+'\n'\
            +'path_to_result_folder=\''+str(path_to_result_folder)+'\'\n'\
            +'expt.plot_all_fit_results(fit_index, path_to_result_folder)'

            nb.create_cell(command, position='below', celltype='code', is_print = True, is_execute = True)

            # Add the markdown cells (cells are separated by the key word 'break')
            split_report='  \n'.join(expt.report).split('break')
            for cell in split_report[::-1]:
                nb.create_cell(cell, position='below', celltype='markdown', is_print = True, is_execute = True)

            # Update logs
            if expt.logs_lvl>=0:
                str_to_add = 'Notebook cells with report of latest fit generated.'
                expt.add_str_to_logs(wid, str_to_add)

    ####################################
    # Define widgets
    ####################################

    wid.button_preview_report.on_click(on_button_preview_report_clicked)
    wid.button_print_report.on_click(on_button_print_report_clicked)


    ####################################
    # Level 1 STATIC
    # Buttons
    ####################################

    out_level_1 = widgets.Output()
    wid_level_1 = widgets.HBox([wid.button_preview_report,
                                wid.button_print_report,
                                wid.fit_index])

    with out_level_1:
        display(wid_level_1)


    ######################################################
    # Level 2 DYNAMIC
    # Report
    ######################################################
    out_level_2 = widgets.Output()


    display(widgets.VBox([out_level_1,out_level_2]))
