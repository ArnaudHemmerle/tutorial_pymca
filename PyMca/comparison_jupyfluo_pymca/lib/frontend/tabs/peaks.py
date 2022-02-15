'''
Module for setting the peaks on the spectrum.
'''
import numpy as np
import ipywidgets as widgets
import ipysheet
import xraylib
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
        object from the class Experiment.
    wid : object myWidgets
        object from the class myWidgets.
    '''
    def set_sheet():
        '''
        Set the sheet and display it.
        '''
        column_headers = ['Peak name', 'Transition name', 'Position (eV)', 'Strength',
                          'Fit position?', 'Fit peak?']
        nb_columns = len(column_headers)

        # Determine the number of rows to have a fixed number of empty rows
        nb_filled_rows = len([elem for elem in expt.peaks_params if elem[0]!=''])
        nb_empty_rows = len([elem for elem in expt.peaks_params if elem[0]==''])
        if nb_empty_rows<15:
            nb_rows = nb_filled_rows+15
        else:
            nb_rows = np.shape(expt.peaks_params)[0]

        # ipysheet does not work correctly with empy cells
        # it is necessary to fill first the cells with something
        nb_rows_to_fill = nb_rows - np.shape(expt.peaks_params)[0]
        fake_cells = np.reshape(np.array(nb_rows_to_fill*nb_columns*['']),
                               (nb_rows_to_fill, nb_columns))
        expt.peaks_params_filled = np.vstack((expt.peaks_params, fake_cells))

        # Create, set, and display the sheet
        wid.sheet = ipysheet.easy.sheet(columns=nb_columns, rows=nb_rows, column_headers = column_headers)
        for i in range(nb_columns):
            ipysheet.easy.column(i,  expt.peaks_params_filled[:,i])

        # Token
        expt.is_sheet_loaded = True

        with out_level_3:
            out_level_3.clear_output()
            display(wid.sheet)


    def on_button_import_peaks_params_file_clicked(b):
        '''
        Import csv file and update values of expt and widgets.
        '''
        if expt.is_extract_done:

            # Choose between default params or a previous set of params
            if wid.select_peaks_params_file.value == 'Default peaks parameters':
                expt.path_to_peaks_params = expt.path_to_peaks_params_default
            else:
                expt.path_to_peaks_params = expt.save_dir+wid.select_peaks_params_file.value+'/peaks_params.csv'

            # Import the csv file and update the sheet
            expt.set_peaks_params_from_file(wid, expt.path_to_peaks_params)
            wid.set_peaks_params_from_expt(expt)

            with out_level_3:
                out_level_3.clear_output()
                set_sheet()

            with out_level_4:
                out_level_4.clear_output()
                display(wid_level_5)

            with out_level_5:
                out_level_5.clear_output()

            with out_level_6:
                out_level_6.clear_output()

            with out_level_7:
                out_level_7.clear_output()

            # Validate the sheet and plot the spectrum with peaks
            on_button_validate_sheet_clicked(b)

        else:

            # Update logs
            if expt.logs_lvl>=0:
                str_to_add = _RED+'Extract a scan first.'+_RESET
                expt.add_str_to_logs(wid, str_to_add)



    def on_button_save_peaks_params_as_default_clicked(b):
        '''
        Save the peaks extraction parameters in the default csv file.
        '''
        if expt.is_sheet_loaded:
            # Update expt with peaks params from sheet and widgets
            expt.set_peaks_params_from_sheet(wid, wid.sheet)
            expt.set_peaks_params_from_widgets(wid)

            # Save peaks params in file
            expt.save_peaks_params_in_file(wid, expt.path_to_peaks_params_default)

        else:
            # Update logs
            if expt.logs_lvl>=0:
                str_to_add = _RED+'Load peaks first.'+_RESET
                expt.add_str_to_logs(wid, str_to_add)

    def on_button_validate_sheet_clicked(b):
        '''
        Validate the sheet and plot the corresponding peaks on the sum of spectrums.
        '''
        expt.validate_sheet(wid)

        with out_level_7:
            out_level_7.clear_output()
            expt.plot_peaks(wid)


    def on_button_select_peaks_from_db_to_load_clicked(b):
        '''
        Display the widgets for adding peaks from the database.
        '''
        def on_button_import_peaks_in_sheet_clicked(b):
            '''
            Update sheet with selected peaks from the database.
            '''
            # Prepare the array to add
            if expt.peaks_params_to_add != []:
                expt.peaks_params = np.vstack((expt.peaks_params, expt.peaks_params_to_add))
                expt.peaks_params_to_add = []

            # Update sheet
            set_sheet()

            if expt.logs_lvl>=0:
                str_to_add = 'Import peaks from the database.'
                expt.add_str_to_logs(wid, str_to_add)

        def on_button_import_peaks_from_db_clicked(b):
            '''
            Extract selected peaks parameters from the database.
            '''
            with out_level_6:
                out_level_6.clear_output()

            if expt.logs_lvl>=1:
                str_to_add = 'Open the database.'
                expt.add_str_to_logs(wid, str_to_add)

            expt.peaks_params_to_add = []

            # Extract selection from widget
            expt.selected_element = wid.select_element.value
            expt.selected_line = wid.select_line.value

            # Based on the numbering in xrayd_lines.pro
            if expt.selected_line == 'K':
                ind_min = -29
                ind_max = 0
            if expt.selected_line == 'L':
                ind_min = -113
                ind_max = -29
            if expt.selected_line == 'M':
                ind_min = -219
                ind_max = -113

            # Get the atom number from the database
            Z = xraylib.SymbolToAtomicNumber(expt.selected_element)

            # Look in the database and extract the info
            for i in range(ind_min,ind_max):

                # try/except to skip the absent lines for a given atom
                try:
                    # energy in xraylib is in keV
                    # energy in this notebook is in eV
                    strength = xraylib.CS_FluorLine_Kissel(Z, i, expt.beam_energy/1000.)
                    energy = xraylib.LineEnergy(Z, i)*1000.

                    # Put an absolute cut-off on the strength
                    if strength>expt.min_strength:
                        expt.peaks_params_to_add.append([
                                    expt.selected_element,
                                    expt.transition_names[-i-1],
                                    '{:.1f}'.format(energy),
                                    '{:f}'.format(strength),'no','yes'])

                        with out_level_6:
                            print('Transition name: %s, energy (eV): %g, strength: %g'\
                                  %(expt.transition_names[-i-1], energy, strength))

                except ValueError:
                    pass

            with out_level_6:
                display(wid.button_import_peaks_in_sheet)


        # Define widgets
        wid.button_import_peaks_from_db.on_click(on_button_import_peaks_from_db_clicked)
        wid.button_import_peaks_in_sheet.on_click(on_button_import_peaks_in_sheet_clicked)


        with out_level_5:
            out_level_5.clear_output()
            display(widgets.HBox([wid.select_element, wid.select_line,
                                    wid.button_import_peaks_from_db]))

        with out_level_6:
            out_level_6.clear_output()

    ####################################
    # Define the widgets
    ####################################

    # Widget for selecting the file
    # Create the list of files available and display them
    expt.set_list_peaks_params_files(wid)
    wid.select_peaks_params_file.options = expt.list_peaks_params_files

    wid.button_import_peaks_params_file.on_click(on_button_import_peaks_params_file_clicked)
    wid.button_save_peaks_params_as_default.on_click(on_button_save_peaks_params_as_default_clicked)
    wid.button_validate_sheet.on_click(on_button_validate_sheet_clicked)
    wid.button_select_peaks_from_db_to_load.on_click(on_button_select_peaks_from_db_to_load_clicked)


    #####################################
    # Level 1 STATIC
    # Peaks params selection
    #####################################

    out_level_1 = widgets.Output()
    wid_level_1 = widgets.VBox([
                    widgets.HBox([
                        wid.select_peaks_params_file,
                        widgets.VBox([wid.button_import_peaks_params_file,
                                      wid.button_save_peaks_params_as_default
                                        ])
                                ])
                              ])


    with out_level_1:
        display(wid_level_1)




    ######################################################
    # Level 2 STATIC
    # Widgets for peaks params
    ######################################################

    out_level_2 = widgets.Output()
    wid_level_3 = widgets.VBox([
                        widgets.Label('Params for conversion to eVs (eVs = gain*channels + eV0):'),
                        widgets.HBox([wid.gain, wid.eV0]),
                        widgets.Label('Params for importing peaks from the database:'),
                        widgets.HBox([wid.beam_energy, wid.min_strength])
                                ],layout={'border': '1px solid black', 'margin': '0px 0px 10px 0px'})



    with out_level_2:
        display(wid_level_3)


    ######################################################
    # Level 3 DYNAMIC
    # Sheet
    ######################################################

    out_level_3 = widgets.Output()


    ######################################################
    # Level 4 DYNAMIC
    # Widgets below sheet
    ######################################################

    out_level_4 = widgets.Output()
    wid_level_5 = widgets.HBox([wid.button_validate_sheet, wid.button_select_peaks_from_db_to_load])


    #########################################################
    # Level 5 DYNAMIC
    # Selection of peaks in database
    #########################################################

    out_level_5 = widgets.Output()

    #########################################################
    # Level 6 DYNAMIC
    # Console of database
    #########################################################

    out_level_6 = widgets.Output()

    #########################################################
    # Level 7 DYNAMIC
    # Plot
    #########################################################

    out_level_7 = widgets.Output()

    display(widgets.VBox([out_level_1, out_level_2, out_level_3,
                          out_level_4, out_level_5, out_level_6, out_level_7]))
