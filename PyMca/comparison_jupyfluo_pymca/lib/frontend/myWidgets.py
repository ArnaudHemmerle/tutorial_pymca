'''
Class regrouping all methods and attributes concerning the widgets.
'''
import os
import ipywidgets as widgets
import xraylib

class myWidgets:
    """
    Class for the widgets.

    Attributes
    ----------
    bckg_eVs_inf : widget FloatText
        Lower bound of the fitting range for the background (if is_bckg_subset).
    bckg_eVs_sup : widget FloatText
        Upper bound of the fitting range for the background (if is_bckg_subset).
    beam_energy : widget FloatText
        Beam energy in eV.
    broad_factor : widget FloatText
        Broadening factor of the Compton peak's width as compared to an elastic peak.
    button_analyze_next_scan : widget Button
        Close the current main widget and open a new one at the bottom of the notebook.
    button_create_all_params_file : widget Button
        Reset session ID and save all params in csv files.
    button_create_extract_params_file : widget Button
        Reset session ID and save extraction parameters in a csv file.
    button_create_fit_params_file  : widget Button
        Reset session ID and save fit parameters in a csv file.
    button_create_peaks_params_file : widget Button
        Reset session ID and save peaks parameters in a csv file.
    button_export : widget Button
        Export to pdf.
    button_extract : widget Button
        Extract the selected scan.
    button_extract_avg : widget Button
        Get averages on fitted parameters and set the corresponding values in the widgets.
    button_import_extract_params_file : widget Button
        Load extraction parameters from file.
    button_import_fit_params_file : widget Button
        Import csv file and update values of expt and widgets.
    button_import_peaks_from_db : widget Button
        Extract selected peaks parameters from the database.
    button_import_peaks_in_sheet :  widget Button
        Update sheet with selected peaks from the database.
    button_import_peaks_params_file : widget Button
        Import csv file and update values of expt and widgets.
    button_preview_report : widget Button
        Preview the report in the corresponding tab.
    button_print_report : widget Button
        Print the report in the Notebook.
    button_refresh : widget Button
        Refresh the list of nexus files.
    button_replace_all_params_file : widget Button
        Save all params in csv files, replacing existing ones.
    button_replace_extract_params_file : widget Button
        Replace the csv file with the current extraction params.
    button_replace_fit_params_file : widget Button
        Replace the csv file with the current fit params.
    button_replace_peaks_params_file : widget Button
        Replace the csv file with the current peaks params.
    button_run_all_fit : widget Button
        Run a fit for each spectrum extracted, and save the results in files.
    button_run_single_fit : widget Button
        Run a fit on the current spectrum.
    button_save_all_params_as_default : widget Button
        Save all params in the default csv files.
    button_save_and_close : widget Button
        Save and delete the current widget (so that the notebook is clean next time it is opened).
    button_save_extract_params_as_default : widget Button
        Save extraction parameters in the default csv file.
    button_save_fit_params_as_default : widget Button
        Save fit parameters in the default csv file.
    button_save_peaks_params_as_default : widget Button
        Save peaks parameters in the default csv file.
    button_select_extract_params_to_load : widget Button
        Display a widget for choosing the extraction params file and load it on click.
    button_select_peaks_from_db_to_load : widget Button
        Display the widgets for adding peaks from the database.
    button_validate_sheet : widget Button
        Validate the peaks displayed in the sheet.
    ct : widget FloatText
        Constant part of the baseline = sl*eVs+ct.
    eV0 : widget FloatText
        Parameter in the conversion of channels to eVs, eVs = channels*gain + eV0.
    first_channel : widget IntText
        Index of the first channel after the extraction.
    first_spectrum : widget IntText
        Index of the first spectrum after the extraction.
    fit_index : widget BoundedIntText
        Index of the spectrum to be fitted.
    gain : widget FloatText
        Parameter in the conversion of channels to eVs, eVs = channels*gain + eV0.
    highTF0 : widget FloatText
        Tail fraction of the peak, on the high energy side.
    highTW0 : widget FloatText
        Tail width of the peak, on the high energy side.
    interactive_select_file : widget interactive
        Interactive selection of the file.
    is_bckg_subset : widget Checkbox
        True if fitting the background on a subset of the spectrum.
    is_broad_factor : widget Checkbox
        True if broad_factor is a fit parameter.
    is_ct : widget Checkbox
        True if ct is a fit parameter.
    is_eV0 : widget Checkbox
        True if eV0 is a fit parameter.
    is_gain : widget Checkbox
        True if gain is a fit parameter.
    is_highTF0 : widget Checkbox
        True if highTF0 is a fit parameter.
    is_highTW0 : widget Checkbox
        True if highTW0 is a fit parameter.
    is_lowTF0 : widget Checkbox
        True if lowTF0 is a fit parameter.
    is_lowTW0 : widget Checkbox
        True if lowTW0 is a fit parameter.
    is_noise : widget Checkbox
        True if noise is a fit parameter.
    is_SDD_elem_0 : widget Checkbox
        True if the SDD element 0 is used.
    is_SDD_elem_1 : widget Checkbox
        True if the SDD element 1 is used.
    is_SDD_elem_2 : widget Checkbox
        True if the SDD element 2 is used.
    is_SDD_elem_3 : widget Checkbox
        True if the SDD element 3 is used.
    is_SDD_elem_4 : widget Checkbox
        True if the SDD element 4 is used.
    is_SF0 : widget Checkbox
        True if SF0 is a fit parameter.
    is_sl : widget Checkbox
        True if sl is a fit parameter.
    is_TF0 : widget Checkbox
        True if TF0 is a fit parameter.
    is_TW0 : widget Checkbox
        True if TW0 is a fit parameter.
    last_channel : widget BoundedIntText
        Index of the last channel after the extraction.
    last_spectrum : widget BoundedIntText
        Index of the last spectrum after the extraction.
    lowTF0 : widget FloatText
        Tail fraction of the peak, on the low energy side.
    lowTW0 : widget FloatText
        Tail width of the peak, on the low energy side.
    min_strength : widget FloatText
        Minimum strength of peak to be displayed in the peak selection widget.
    noise : widget FloatText
        Width of the peak in keV, before contribution from the detector.
    out_logs : widget Output
        Window for the logs.
    plot_frequency : widget BoundedIntText
        During a series of fit, plot the result every X spectrums.
    select_element : widget Dropdown
        Selection list of elements in the database.
    select_extract_params_file : widget Select
        Selection list of csv files with extraction params.
    select_fit_params_file : widget Select
        Selection list of csv files with fit params.
    select_line : widget Dropdown
        Select a fluorescence line K, L or M.
    select_file : widget Select
        Selection list of nexus files.
    select_peaks_params_file : widget Select
        Selection list of csv files with peak params.
    session_id : widget HTML
        Display the session ID.
    SF0 : widget FloatText
        Shelf fraction of the peak.
    sheet : widget Sheet
        Interactive sheet with the peaks parameters, for the module ipysheet.
    show_bckg : widget Checkbox
        True if the background part should be plot with the fit.
    show_compton : widget Checkbox
        True if the Compton part should be plot with the fit.
    show_gaussian : widget Checkbox
        True if the gaussian part should be plot with the fit.
    show_shelf : widget Checkbox
        True if the shelf part should be plot with the fit.
    show_tail : widget Checkbox
        True if the tail part should be plot with the fit.
    sl : widget FloatText
        Linear part of the baseline = sl*eVs+ct.
    tab : widget Tab
        Main tab.
    TF0 : widget FloatText
        Tail fraction of the peak.
    TW0 : widget FloatText
        Tail width of the peak.

    Methods
    -------
    __init__():
        Constructor.
    define_widgets_bar():
        Define the widgets used in the top bar.
    define_widgets_extract():
        Define the widgets used in the extraction.
    define_widgets_fit():
        Define the widgets used in the fit.
    define_widgets_main():
        Define the widgets used in the main file.
    define_widgets_peaks():
        Define the widgets used in peaks definition.
    define_widgets_report():
        Define the widgets used in the report tab.
    set_default_peak_value(expt):
        Set the option of the widget select_peaks_params with the most recent peak file from the current scan.
    set_extract_params_from_expt(expt):
        Update the values of the widgets from the extraction values of expt.
    set_fit_params_from_expt(expt):
        Update the values of the widgets from the fit values of expt.
    set_list_files(expt):
        Update the widget select_file.
    set_peaks_params_from_expt(expt):
        Update the values of the widgets from the peaks values of expt.
    set_session_id(expt):
        Update the widget session_id.
    """
    def __init__(self):
        '''
        Constructor.
        '''
        self.define_widgets_main()
        self.define_widgets_bar()
        self.define_widgets_extract()
        self.define_widgets_fit()
        self.define_widgets_peaks()
        self.define_widgets_report()


    def define_widgets_main(self):
        '''
        Define the widgets used in the main file.
        '''
        self.session_id = widgets.HTML()


    def define_widgets_bar(self):
        '''
        Define the widgets used in the top bar.
        '''
        self.button_export = widgets.Button(description="Export to PDF",
                                            layout=widgets.Layout(width='150px'))

        self.button_analyze_next_scan = widgets.Button(description="Analyze next scan",
                                            layout=widgets.Layout(width='150px'))

        self.button_extract_avg = widgets.Button(description="Extract averages from fit",
                                            layout=widgets.Layout(width='200px'))

        self.button_save_and_close = widgets.Button(description="Save & close the notebook",
                                            layout=widgets.Layout(width='200px'))

        self.button_save_all_params_as_default = widgets.Button(description="Save all params as default",
                                            layout=widgets.Layout(width='200px'))

        self.button_replace_all_params_file = widgets.Button(description="Replace files")

        self.button_create_all_params_file = widgets.Button(description="Create new files")

    def define_widgets_extract(self):
        '''
        Define the widgets used in the extraction.
        '''
        self.is_SDD_elem_0 = widgets.Checkbox(
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='50px'),
                    description='0')

        self.is_SDD_elem_1 = widgets.Checkbox(
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='50px'),
                    description='1')

        self.is_SDD_elem_2 = widgets.Checkbox(
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='50px'),
                    description='2')

        self.is_SDD_elem_3 = widgets.Checkbox(
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='50px'),
                    description='3')

        self.is_SDD_elem_4 = widgets.Checkbox(
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='50px'),
                    description='4')

        self.first_channel = widgets.BoundedIntText(
                    min=0,
                    max=2047,
                    step=1,
                    description='First channel',
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='200px'))

        self.last_channel = widgets.BoundedIntText(
                    min=0,
                    max=2047,
                    step=1,
                    description='Last channel',
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='200px'))

        self.first_spectrum = widgets.IntText(
                    step=1,
                    description='First spectrum',
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='200px'))

        self.last_spectrum = widgets.IntText(
                    step=1,
                    description='Last spectrum',
                    style = {'description_width': 'initial'},
                    layout=widgets.Layout(width='200px'))

        self.button_refresh = widgets.Button(description="Refresh list",
                                        layout=widgets.Layout(width='250px'))

        self.button_extract = widgets.Button(description="Extract scan",
                                        layout=widgets.Layout(width='250px'))

        self.button_select_extract_params_to_load = widgets.Button(description="Import extract params",
                                            layout=widgets.Layout(width='250px'))

        self.button_save_extract_params_as_default = widgets.Button(description="Save extract params as default",
                                                    layout=widgets.Layout(width='250px'))

        self.button_import_extract_params_file = widgets.Button(description="Import selected params",
                                                               layout=widgets.Layout(width='190px'))

        self.button_replace_extract_params_file = widgets.Button(description="Replace file")

        self.button_create_extract_params_file = widgets.Button(description="Create new file")

        self.select_file = widgets.Select(
                                  options=['Empty'],
                                  rows=10,
                                  description=' ',
                                  style = {'description_width': 'initial'},
                                  layout=widgets.Layout(width='660px'))

        self.select_extract_params_file = widgets.Select(
                            options=['Empty'],
                            rows=10,
                            description=' ',
                            style = {'description_width': 'initial'},
                            layout=widgets.Layout(width='700px'))


    def define_widgets_peaks(self):
        '''
        Define the widgets used in the peaks definition.
        '''
        self.button_import_peaks_params_file = widgets.Button(description="Import peaks params",
                                                layout=widgets.Layout(width='250px'))

        self.button_save_peaks_params_as_default = widgets.Button(description="Save peaks params as default",
                                                layout=widgets.Layout(width='250px'))

        self.select_peaks_params_file = widgets.Select(
                                  options=['Default peaks parameters'],
                                  rows=10,
                                  description=' ',
                                  style = {'description_width': 'initial'},
                                  layout=widgets.Layout(width='660px'))


        self.button_replace_peaks_params_file = widgets.Button(description="Replace file")

        self.button_create_peaks_params_file = widgets.Button(description="Create new file")

        self.button_validate_sheet = widgets.Button(description="Validate sheet")

        self.button_select_peaks_from_db_to_load = widgets.Button(description="Add peaks from database",
                                                       layout=widgets.Layout(width='200px'))

        self.select_element = widgets.Dropdown(
                                options=[str(xraylib.AtomicNumberToSymbol(i)) for i in range(1,99)],
                                value='Ar',
                                description='Select element:',
                                layout=widgets.Layout(width='200px'),
                                style={'description_width': 'initial'})

        self.button_import_peaks_in_sheet = widgets.Button(description="Import peaks in sheet",
                                                           layout=widgets.Layout(width='200px'))


        self.button_import_peaks_from_db = widgets.Button(description="Get peaks info",
                                                       layout=widgets.Layout(width='200px'))


        self.select_line = widgets.Dropdown(
                                options=['K', 'L', 'M'],
                                description='Select line:',
                                value='K',
                                layout=widgets.Layout(width='200px'),
                                style={'description_width': 'initial'})

        self.beam_energy = widgets.FloatText(
                                description='Energy (eV)',
                                layout=widgets.Layout(width='150px'),
                                style={'description_width': 'initial'})

        self.min_strength = widgets.FloatText(
                                    description='Minimum peak strength',
                                    layout=widgets.Layout(width='200px'),
                                    style={'description_width': 'initial'})

        self.gain = widgets.FloatText(
            description='gain',
            layout=widgets.Layout(width='120px'),
            style={'description_width': 'initial'})

        self.eV0 = widgets.FloatText(
                description='eV0',
                layout=widgets.Layout(width='100px'),
                style={'description_width': 'initial'})



    def define_widgets_fit(self):
        '''
        Define the widgets used in the fit.
        '''
        self.button_import_fit_params_file = widgets.Button(description="Import fit params",
                                                layout=widgets.Layout(width='250px'))

        self.button_save_fit_params_as_default = widgets.Button(description="Save fit params as default",
                                                layout=widgets.Layout(width='250px'))

        self.select_fit_params_file = widgets.Select(
                                  options=['Default peaks parameters'],
                                  rows=10,
                                  description=' ',
                                  style = {'description_width': 'initial'},
                                  layout=widgets.Layout(width='660px'))


        self.button_replace_fit_params_file = widgets.Button(description="Replace file")

        self.button_create_fit_params_file = widgets.Button(description="Create new file")

        self.button_run_single_fit = widgets.Button(description="Preview fit on single spectrum",
                                              layout=widgets.Layout(width='300px'))

        self.button_run_all_fit = widgets.Button(description="Run fit on all spectrums",
                                              layout=widgets.Layout(width='300px'))

        self.button_replace_all_files_and_fit = widgets.Button(description="Replace files")

        self.button_create_all_files_and_fit = widgets.Button(description="Create new files")


        self.fit_index = widgets.BoundedIntText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='Preview spectrum X:',
            min=0)

        self.plot_frequency = widgets.BoundedIntText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='Plot fit every X spectrum:',
            value=10,
            min=1)


        self.is_gain = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_eV0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_lowTW0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_highTW0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_lowTF0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_highTF0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_broad_factor = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_TW0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_TF0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_SF0 = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_sl = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_ct = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.is_noise = widgets.Checkbox(
            style={'description_width': '0px'},
            description=' ')

        self.lowTW0 = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='300px'),
            description='tail width (low energy side)')

        self.highTW0 = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='300px'),
            description='tail width (high energy side)')

        self.lowTF0 = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='300px'),
            description='tail fraction (low energy side)')

        self.highTF0 = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='300px'),
            description='tail fraction (high energy side)')

        self.broad_factor = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='broadening factor')

        self.TW0 = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='300px'),
            description='tail width (low energy side)')

        self.TF0 = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='300px'),
            description='tail fraction (low energy side)')

        self.SF0 = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='shelf fraction')

        self.sl = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='slope')

        self.ct = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='constant')

        self.noise = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='150px'),
            description='noise')

        self.show_gaussian = widgets.Checkbox(
            style={'description_width': 'initial'},
            description='Show gaussian part',
            value=True)

        self.show_shelf = widgets.Checkbox(
            style={'description_width': 'initial'},
            description='Show shelf part',
            value=True)

        self.show_tail = widgets.Checkbox(
            style={'description_width': 'initial'},
            description='Show tail part',
            value=True)

        self.show_compton = widgets.Checkbox(
            style={'description_width': 'initial'},
            description='Show Compton part',
            value=True)

        self.show_bckg = widgets.Checkbox(
            style={'description_width': 'initial'},
            description='Show background',
            value=True)

        self.bckg_eVs_inf = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='lower bound (eV)')    
        
        self.bckg_eVs_sup = widgets.FloatText(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='upper bound (eV)')          

        self.is_bckg_subset = widgets.Checkbox(
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='200px'),
            description='Fit background on a subset:')        
        
    def set_extract_params_from_expt(self, expt):
        '''
        Update the values of the widgets from the extraction values of expt.

        Parameters
        ----------
        expt : object Experiment
            object from the class Experiment
        '''
        self.is_SDD_elem_0.value = expt.is_SDD_elem_0
        self.is_SDD_elem_1.value = expt.is_SDD_elem_1
        self.is_SDD_elem_2.value = expt.is_SDD_elem_2
        self.is_SDD_elem_3.value = expt.is_SDD_elem_3
        self.is_SDD_elem_4.value = expt.is_SDD_elem_4
        self.first_channel.value = expt.first_channel
        self.last_channel.value = expt.last_channel
        self.first_spectrum.value = expt.first_spectrum
        self.last_spectrum.value = expt.last_spectrum

        if expt.logs_lvl>=1:
            str_to_add = 'Update extraction parameters widget.'
            expt.add_str_to_logs(self, str_to_add)

    def set_fit_params_from_expt(self, expt):
        '''
        Update the values of the widgets from the fit values of expt.

        Parameters
        ----------
        expt : object Experiment
            object from the class Experiment
        '''
        self.gain.value = expt.gain
        self.eV0.value = expt.eV0
        self.sl.value = expt.sl
        self.ct.value = expt.ct
        self.noise.value = expt.noise
        self.SF0.value = expt.SF0
        self.TF0.value = expt.TF0
        self.TW0.value = expt.TW0
        self.broad_factor.value = expt.broad_factor
        self.lowTF0.value = expt.lowTF0
        self.highTF0.value = expt.highTF0
        self.lowTW0.value = expt.lowTW0
        self.highTW0.value = expt.highTW0
        self.bckg_eVs_inf.value = expt.bckg_eVs_inf
        self.bckg_eVs_sup.value = expt.bckg_eVs_sup

        self.is_gain.value = 'gain' in expt.list_isfit
        self.is_eV0.value = 'eV0' in expt.list_isfit
        self.is_sl.value = 'sl' in expt.list_isfit
        self.is_ct.value = 'ct' in expt.list_isfit
        self.is_noise.value = 'noise' in expt.list_isfit
        self.is_SF0.value = 'SF0' in expt.list_isfit
        self.is_TF0.value = 'TF0' in expt.list_isfit
        self.is_TW0.value = 'TW0' in expt.list_isfit
        self.is_broad_factor.value = 'broad_factor' in expt.list_isfit
        self.is_lowTF0.value = 'lowTF0' in expt.list_isfit
        self.is_highTF0.value = 'highTF0' in expt.list_isfit
        self.is_lowTW0.value = 'lowTW0' in expt.list_isfit
        self.is_highTW0.value = 'highTW0' in expt.list_isfit
        self.is_bckg_subset.value = expt.is_bckg_subset

        if expt.logs_lvl>=1:
            str_to_add = 'Update fit parameters widget.'
            expt.add_str_to_logs(self, str_to_add)

    def set_peaks_params_from_expt(self, expt):
        '''
        Update the values of the widgets from the peaks values of expt.

        Parameters
        ----------
        expt : object Experiment
            object from the class Experiment
        '''
        self.beam_energy.value = expt.beam_energy
        self.min_strength.value = expt.min_strength
        self.gain.value = expt.gain
        self.eV0.value = expt.eV0

        if expt.logs_lvl>=1:
            str_to_add = 'Update peaks parameters widget.'
            expt.add_str_to_logs(self, str_to_add)

    def set_session_id(self, expt):
        '''
        Update the widget session_id.

        Parameters
        ----------
        expt : object Experiment
            object from the class Experiment
        '''
        self.session_id.value = f"<h style=\"font-size:20px;\"><b>{'Session ID: '+expt.session_id}</b></h> "

        if expt.logs_lvl>=1:
            str_to_add = 'Update session id widget.'
            expt.add_str_to_logs(self, str_to_add)

    def set_list_files(self, expt):
        '''
        Update the widget select_file.

        Parameters
        ----------
        expt : object Experiment
            object from the class Experiment
        '''
        self.select_file.options = expt.list_files_filenames

        # Update logs
        if expt.logs_lvl>=1:
            str_to_add = 'Update nexus file selection widget.'
            expt.add_str_to_logs(self, str_to_add)


    def set_default_peak_value(self, expt):
        '''
        Set the option of the widget select_peaks_params_file with the most recent peak file from the current scan.

        Parameters
        ----------
        expt : object Experiment
            object from the class Experiment
        '''
        # Create the list of existing peaks files from the current scan
        list_current_id_files = []
        for root, _, files in os.walk(expt.save_dir, topdown=True):
            for name in files:
                if 'peaks_params' in name:
                    path_to_csv = os.path.join(root.split('/')[-2],root.split('/')[-1])
                    if expt.name in root:
                        list_current_id_files.append(path_to_csv)
        list_current_id_files.sort(reverse=False)

        # Get the most recent peak file from the current scan
        if list_current_id_files == []:
            self.select_peaks_params_file.value = 'Default peaks parameters'
        else:
            self.select_peaks_params_file.value = list_current_id_files[-1]

        # Update logs
        if expt.logs_lvl>=1:
            str_to_add = 'Set default value for peaks file selection widget.'
            expt.add_str_to_logs(self, str_to_add)

    def set_default_fit_value(self, expt):
        '''
        Set the option of the widget select_fit_params_file with the most recent fit file from the current scan.

        Parameters
        ----------
        expt : object Experiment
            object from the class Experiment
        '''
        # Create the list of existing fit files from the current scan
        list_current_id_files = []
        for root, _, files in os.walk(expt.save_dir, topdown=True):
            for name in files:
                if 'fit_params' in name:
                    path_to_csv = os.path.join(root.split('/')[-2],root.split('/')[-1])
                    if expt.name in root:
                        list_current_id_files.append(path_to_csv)
        list_current_id_files.sort(reverse=False)

        # Get the most recent fit file from the current scan
        if list_current_id_files == []:
            self.select_fit_params_file.value = 'Default fit parameters'
        else:
            self.select_fit_params_file.value = list_current_id_files[-1]

        # Update logs
        if expt.logs_lvl>=1:
            str_to_add = 'Set default value for fit file selection widget.'
            expt.add_str_to_logs(self, str_to_add)

    def define_widgets_report(self):
        '''
        Define the widgets used in the report tab.
        '''
        self.button_preview_report = widgets.Button(description="Preview report",
                                                    layout=widgets.Layout(width='200px'))


        self.button_print_report = widgets.Button(description="Print report",
                                                  layout=widgets.Layout(width='200px'))
