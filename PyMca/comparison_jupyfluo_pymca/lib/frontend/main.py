'''
Main module of the frontend. Set the session and widgets.
'''
from IPython.display import display
import matplotlib.pyplot as plt
import ipywidgets as widgets

from lib.frontend import bar as wbar
from lib.frontend.tabs import extract as wextract
from lib.frontend.tabs import fit as wfit
from lib.frontend.tabs import peaks as wpeaks
from lib.frontend.tabs import report as wreport

# General parameters for the plots
plot_params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
plt.rcParams.update(plot_params)

def start_session(expt, wid):
    '''
    Start a new session.

    Parameters
    ----------
    expt : object Experiment
        object from the class Experiment.
    wid : object myWidgets
        object from the class myWidgets.
    '''
    ##################################
    # Settings of the logs window
    ##################################
    wid.out_logs = widgets.Output(layout={'border': '1px solid black', \
                                          'max_height': '150px', 'overflow': 'scroll'})

    if expt.logs_lvl>=0:
        expt.add_str_to_logs(wid, 'Session starts.')

    ##################################
    # Initialize the session
    ##################################

    # Generate a session id based on time
    expt.set_session_id(wid)
    wid.set_session_id(expt)

    # Set the list of files in the files directory
    expt.set_list_files(wid)
    wid.set_list_files(expt)

    ##################################
    # Settings of the top bar
    ##################################
    out_bar = widgets.Output()
    with out_bar:
        wbar.create_bar(expt, wid)

    ##################################
    # Settings of the tabs
    ##################################
    # Tab "Extract"
    out_extract = widgets.Output()
    with out_extract:
        wextract.create_tab(expt, wid)


    # Tab "Peaks"
    out_peaks = widgets.Output()
    with out_peaks:
        wpeaks.create_tab(expt, wid)


    # Tab "Fit"
    out_fit = widgets.Output()
    with out_fit:
        wfit.create_tab(expt, wid)


    # Tab "Report"
    out_report = widgets.Output()
    with out_report:
        wreport.create_tab(expt, wid)

    # Create the main widget tab
    wid.tab = widgets.Tab(children=[out_extract, out_peaks, out_fit, out_report])

    wid.tab.set_title(0, 'Extraction')
    wid.tab.set_title(1, 'Peaks')
    wid.tab.set_title(2, 'Fit')
    wid.tab.set_title(3, 'Report')


    ##################################
    # Display of the widgets
    ##################################
    display(widgets.VBox([wid.session_id, wid.out_logs, out_bar, wid.tab]))
