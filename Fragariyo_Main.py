"""
IMAnnotatorv3: PeakMatcher for Terminal Fragments
Author: Carolina Rojas Ramirez
Date: 05/22/2020
In silico fragmentation of proteins using mass.fast_mass2 from pyteomics
Pieces of code from Daniel A. Polasky's IMAnnotatorv2
Pieces of code from Carolina Rojas Ramirez' Fragmentor (internal fragment analysis + Disulfinator)
"""


import terminalFragmentor_Main
import InternalFragmentor_Main
import Findpeaks
import os
import sys
import pygubu
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import messagebox
import subprocess
import multiprocessing
import logging
from logging.handlers import RotatingFileHandler
import SimpleToolTip
import matplotlib
import RenameIMTBXoutputs
import Modifications
import OutputAnalysis_v2
import PeakPicking_and_Deisotoping



matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.use('Agg')

# Load resource file paths, supporting both live code and code bundled by PyInstaller
if getattr(sys, 'frozen', False):
    root_dir = sys._MEIPASS
    program_data_dir = r"C:\Program Files (x86)\fragmentor_build"
else:
    root_dir = os.path.dirname(__file__)
    program_data_dir = root_dir

guardpath = os.path.join(program_data_dir, 'guardfile.txt')
log_file = os.path.join(program_data_dir, 'Fragariyo.log')
PROJECT_PATH = os.path.join(root_dir, 'UI')
hard_file_path_ui = os.path.join(root_dir, 'UI', 'Fragariyo_GUI.ui')
hard_tooltips_file = os.path.join(root_dir, 'tooltips.txt')
help_file = os.path.join(root_dir, 'Fragariyo_SOP.pdf')
about_file = os.path.join(root_dir, 'README.txt')
hard_datanalysis_ui = os.path.join(root_dir, 'UI', 'DataAnalysis_GUI.ui')

# print(hard_file_path_ui)

class Fragariyo(object):
    """
    Primary graphical class for running CIU2 user interface. Uses PyGuBu builder to create
    interface and handles associated events.
    """
    def __init__(self, tk_root_window):
        """
        Create a new GUI window, connect feedback to buttons, and load parameters
        """
        self.tk_root = tk_root_window

        # create a Pygubu builder
        self.builder = builder = pygubu.Builder()
        builder.add_resource_path(PROJECT_PATH)

        self.builder.add_from_file(hard_file_path_ui)


        # create widget using provided root (Tk) window
        self.mainwindow = builder.get_object('CIU_app_top')
        self.mainwindow.protocol('WM_DELETE_WINDOW', self.on_close_window)


        callbacks = {
            'on_button_help_clicked': self.on_button_help_clicked,
            'on_button_about_clicked': self.on_button_about_clicked,

            #IMTBX output renaming
            'on_button_IMTBXrenamed_clicked': self.on_button_IMTBXrenamed_clicked,


            #Now it is Modifications Library upload button
            'on_button_modsrepo_clicked': self.on_button_modsrepo_clicked,



            'on_button_changedir_clicked': self.on_button_changedir_clicked,


            'on_terminal_run_clicked': self.on_button_terminal_run_clicked,
            'on_internal_run_clicked': self.on_button_internal_run_clicked,
            'on_button_importxy_clicked': self.on_button_importxy_clicked,

            'on_button_dataanalysis_clicked':self.on_button_dataanalysis_clicked,
            'on_singleseqcov_clicked': self.singleseqcoverage,
            'on_combinedeqcov_clicked': self.comboseqcoverage,
            'on_merginginternal_clicked':self.combining_internalouttsv,
            'on_averaginginternal_clicked':self.averaging_internalouttsv,
            'on_termfragplotterinput_clicked':self.on_termfragplotterinput_clicked,
            'on_internfragplotter_clicked':self.on_internfragplotter_clicked,
            "on_button_peakpickdeiso_fromraw_clicked": self.on_button_peakpickdeiso_fromraw_clicked


        }
        builder.connect_callbacks(callbacks)

        # Peak picking extra params
        self.noisesections_val = self.builder.get_object('numberofsectionsvalue_entry')
        self.signalabovenoise_val = self.builder.get_object('timesabovenoise_entry')
        # Internal fragment extra params
        self.massresolution = self.builder.get_object('massresolution')
        self.toleranceerror = self.builder.get_object('toleranceerror')

        self.initialize_tooltips()


        # holder for feature information in between user assessment - plan to replace with better solution eventually
        self.temp_feature_holder = None

        self.analysis_file_list = []
        self.output_dir = root_dir
        self.output_dir_override = False
        self.modificationsrepo = {}





    def run(self):
        """
        Run the main GUI loop
        :return: void
        """
        self.mainwindow.mainloop()

    def on_close_window(self):
        """
        Close (destroy) the app window and the Tkinter root window to stop the process.
        :return: void
        """
        self.mainwindow.destroy()
        self.tk_root.destroy()

    def initialize_tooltips(self):
        """
        Register tooltips for all buttons/widgets that need tooltips
        :return: void
        """
        tooltip_dict = parse_tooltips_file(hard_tooltips_file)
        for tip_key, tip_value in tooltip_dict.items():
            SimpleToolTip.create(self.builder.get_object(tip_key), tip_value)



    def on_button_help_clicked(self):
        """
        Open the manual
        :return: void
        """
        subprocess.Popen(help_file, shell=True)

    def on_button_about_clicked(self):
        """
        Open the license/about information file
        :return: void
        """
        subprocess.Popen(about_file, shell=True)


    def on_button_IMTBXrenamed_clicked(self):
        """
        Open a filechooser for the user to select raw files, then process them
        :return:
        """
        self.progress_started()

        RenameIMTBXoutputs.main()

        self.progress_done()

    def on_button_importxy_clicked(self):
        """

        :return:
        """

        self.progress_started()

        Findpeaks.isotope_xtractor(main_outdir=self.output_dir)

        self.progress_done()

    def on_button_dataanalysis_clicked(self):
        """

        :return:
        """

        self.progress_started()

        watersimport_app = OutputAnalysis_v2.DataAnalysisUI(hard_datanalysis_ui)
        returncode = watersimport_app.run()

        if returncode:
            self.progress_done()

        self.progress_done()



    def on_button_changedir_clicked(self):
        """
        Open a file chooser to change the output directory and update the display
        :return: void
        """
        newdir = filedialog.askdirectory()
        if newdir == '':
            # user hit cancel - don't set directory
            return
        self.output_dir = newdir
        self.output_dir_override = True     # stop changing directory on file load if the user has specified a directory
        self.update_dir_entry()

    def update_dir_entry(self):
        """
        Update the graphical display of the output directory
        :return: void
        """
        self.builder.get_object('Text_outputdir').config(state=tk.NORMAL)
        self.builder.get_object('Text_outputdir').delete(1.0, tk.INSERT)
        self.builder.get_object('Text_outputdir').insert(tk.INSERT, self.output_dir)
        self.builder.get_object('Text_outputdir').config(state=tk.DISABLED)


    def on_button_terminal_run_clicked(self):
        """
        Run terminal fragmentor
        :return: saves new .ciu files
        """
        self.progress_started()

        terminalFragmentor_Main.main_batch_multipass(main_outdir=self.output_dir,
                                                     modificationsrepo = self.modificationsrepo)

        self.progress_done()

    def on_button_internal_run_clicked(self):
        """
        Run terminal fragmentor
        :return: saves new .ciu files
        """

        self.progress_started()

        InternalFragmentor_Main.internal_main_batch_multipass(main_outdir=self.output_dir,repoofmods=self.modificationsrepo, massresolution=self.massresolution.get(),error_ppm=self.toleranceerror.get())


        self.progress_done()



    def progress_print_text(self, text_to_print, prog_percent):
        """
        Set the progress bar to print text (e.g. 'analysis in progress...')
        :param text_to_print: text to display in the progress bar
        :param prog_percent: percent to fill the progress bar
        :return: void
        """
        # update progress
        progress_bar = self.builder.get_object('Progressbar_main')
        progress_bar['value'] = prog_percent
        # display text
        self.builder.get_object('Entry_progress').config(state=tk.NORMAL)
        self.builder.get_object('Entry_progress').delete(0, tk.END)
        self.builder.get_object('Entry_progress').insert(0, text_to_print)
        self.builder.get_object('Entry_progress').config(state=tk.DISABLED)
        # refresh display
        self.mainwindow.update()

    def progress_started(self):
        """
        Display a message showing that the operation has begun
        :return: void
        """
        self.builder.get_object('Entry_progress').config(state=tk.NORMAL)
        self.builder.get_object('Entry_progress').delete(0, tk.END)
        self.builder.get_object('Entry_progress').insert(0, 'Processing... (This window will not respond until processing completes!)')
        self.builder.get_object('Entry_progress').config(state=tk.DISABLED)

        self.builder.get_object('Progressbar_main')['value'] = 10
        self.mainwindow.update()

    def update_progress(self, current_analysis, num_analyses):
        """
        Update the progress bar to display the current progress through the analysis list
        :param current_analysis: the file NUMBER currently being worked on by the program
        :param num_analyses: the total number of files in the current analysis
        :return: void
        """
        current_prog = (current_analysis + 1) / float(num_analyses) * 100
        prog_string = 'Processed {} of {}'.format(current_analysis + 1, num_analyses)

        progress_bar = self.builder.get_object('Progressbar_main')
        progress_bar['value'] = current_prog

        self.builder.get_object('Entry_progress').config(state=tk.NORMAL)
        self.builder.get_object('Entry_progress').delete(0, tk.END)
        self.builder.get_object('Entry_progress').insert(0, prog_string)
        self.builder.get_object('Entry_progress').config(state=tk.DISABLED)
        self.mainwindow.update()

    def progress_done(self):
        """
        Called after methods finished to ensure GUI mainloop continues
        :return: void
        """
        self.builder.get_object('Entry_progress').config(state=tk.NORMAL)
        self.builder.get_object('Entry_progress').delete(0, tk.END)
        self.builder.get_object('Entry_progress').insert(0, 'Done!')
        self.builder.get_object('Entry_progress').config(state=tk.DISABLED)

        self.builder.get_object('Progressbar_main')['value'] = 100
        # added to keep program from exiting when run from command line - not sure if there's a better way to do this, but seems to work
        self.run()
        # return 0

    def open_files(self, filetype):
        """
        Open a tkinter filedialog to choose files of the specified type. NOTE: withdraws the main window of the
        UI (self) to prevent users from clicking on any other GUI elements while the filedialog is active.
        :param filetype: filetype filter in form [(name, extension)]
        :return: list of selected files
        """
        self.mainwindow.withdraw()
        files = filedialog.askopenfilenames(filetypes=filetype)

        # check for excessively long paths and warn the user to shorten them or face possible crashes
        warn_flag = False
        for file in files:
            if len(file) > 200:
                warn_flag = True
        if warn_flag:
            messagebox.showwarning('Warning! Long File Path', 'Warning! At least one loaded file has a path length greater than 200 characters. Windows will not allow you to save files with paths longer than 260 characters. It is strongly recommended to shorten the name of your file and/or the path to the folder containing it to prevent crashes if the length of analysis files/outputs exceed 260 characters.')

        self.mainwindow.deiconify()
        return files

    def on_button_modsrepo_clicked(self):
        """
        Open the Agilent CIU extractor app with a specified output directory. The extractor app will generate _raw.csv
        files into the specified directory. Then runs Agilent data converter script to handle duplicate DTs and
        edit the CV header to be the correct values. Finally, loads data into CIUSuite 2 using the standard load
        raw files method.
        :return: void
        """

        self.progress_started()

        modsinputfile = filedialog.askopenfilename(title='Modification Info', filetypes=[('Tab-separated Mods Repository File', '.txt')])

        self.modificationsrepo = Modifications.modsrepo_creator(modsinputfile)

        self.progress_done()
        return

    def singleseqcoverage(self):
        """
        Terminal Fragments Output Analysis Method
        Plots sequence coverage base on fragment intensity
        :return:
        """
        self.progress_started()
        hitsfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('Hits', '.hits')])

        OutputAnalysis_v2.frips_main_seq_cov(hitsfiles, seqcov=True, seqcov_iontype=False, combination=False)
        self.progress_done()
        return

    def comboseqcoverage(self):
        """
        Terminal Fragments Putput Analysis Method
        :return:
        """
        hitsfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('Hits', '.hits')])

        OutputAnalysis_v2.frips_main_seq_cov(hitsfiles, seqcov=False, seqcov_iontype=False, combination=True)

    def on_termfragplotterinput_clicked(self):
        """
        Terminal Fragments Sequence Coverage Plotter, ClipsMS style
        :return:
        """

        hitsfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('Hits', '.hits')])

        OutputAnalysis_v2.ClipsMS_FragmentorPipe(hitsfiles,self.output_dir)

        inputfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('CSV', '_clipsfig3.csv')])

        OutputAnalysis_v2.fragment_plotter(inputfiles, outputdir=self.output_dir,internal=False)

    def combining_internalouttsv(self):
        """
        Internal Fragments Putput Analysis Method
        :return:
        """

        tsvfiles = filedialog.askopenfilenames(title='InternalFrags', filetypes=[('Internal Matches', '.tsv')])

        OutputAnalysis_v2.merging_interalfrag_tsv(tsvfiles, avg=None)

    def averaging_internalouttsv(self):
        """
        Internal Fragments Putput Analysis Method
        :return:
        """

        tsvfiles = filedialog.askopenfilenames(title='InternalFrags', filetypes=[('Internal Matches', '.tsv')])

        OutputAnalysis_v2.merging_interalfrag_tsv(tsvfiles, avg=True)

    def on_internfragplotter_clicked(self):
        """
        Internal Fragments Plotter (taken from CLipsMS: "Don't fix what is not broken")
        :return:
        """

        inputfiles = filedialog.askopenfilenames(title='Load Hits Files',
                                                 filetypes=[('Internal Fragments', '.tsv')])
        OutputAnalysis_v2.fragment_plotter(inputfiles,outputdir=self.output_dir, internal=True)

    def on_button_peakpickdeiso_fromraw_clicked(self):
        """
        Function to convert raw data to a peaklist
        """

        PeakPicking_and_Deisotoping.peaklist_fromrawcsv(int(self.noisesections_val.get()), int(self.signalabovenoise_val.get()), 190)





class CIU2ConsoleFormatter(logging.Formatter):
    """
    Custom format handler for printing info/warnings/errors with various formatting to console
    """
    info_format = '%(message)s'     # for INFO level, only pass the message contents
    err_format = '[%(levelname)s]: %(message)s'     # for warnings and errors, pass additional information

    def __init__(self):
        super().__init__()

    def format(self, record):
        """
        override formatter to allow for custom formats for various levels
        :param record: log record to format
        :return: formatted string
        """
        # save original format to restore later
        original_fmt = self._style._fmt

        if record.levelno == logging.INFO or record.levelno == logging.DEBUG:
            self._style._fmt = CIU2ConsoleFormatter.info_format

        elif record.levelno == logging.WARNING or record.levelno == logging.ERROR or record.levelno == logging.CRITICAL:
            self._style._fmt = CIU2ConsoleFormatter.err_format

        output = logging.Formatter.format(self, record)
        self._style._fmt = original_fmt

        return output


# ****** CIU Main I/O methods ******
def save_existing_output_string(full_output_path, string_to_save):
    """
    Write an existing (e.g. combined from several objects) string to file
    :param full_output_path: full path filename in which to save output
    :param string_to_save: string to save
    :return: void
    """
    try:
        with open(full_output_path, 'w') as outfile:
            outfile.write(string_to_save)
    except PermissionError:
        messagebox.showerror('Please Close the File Before Saving',
                             'The file {} is open or being used by another program! Please close it, THEN press the OK button to retry saving'.format(full_output_path))
        with open(full_output_path, 'w') as outfile:
            outfile.write(string_to_save)


def parse_subclass_inputs():
    """
    Helper method to handle user input for subclass label strings
    :return: list of strings
    """
    subclass_string = simpledialog.askstring('Enter State Labels, Separated by Commas',
                                             'Enter the Labels for each State, separated by commas. NOTE: the labels entered must EXACTLY match a part of the file name for each state or the data will not be loaded properly!')
    if subclass_string is not None:
        subclass_splits = subclass_string.split(',')
        subclass_labels = []
        for split in subclass_splits:
            subclass_labels.append(split.strip())
    else:
        subclass_labels = []
    if len(subclass_labels) < 2:
        logger.warning('No (or only 1) state labels read from input! Classification will NOT use states')
    return subclass_labels


def parse_tooltips_file(tooltips_file):
    """
    Parse the tooltips.txt file for all tooltips to display in the GUI. Returns a dictionary
    with all object names and descriptions to display
    :param tooltips_file: File to parse (tab-delimited text, headers = #)
    :return: Dictionary, key=component name to pass to Pygubu, value=description
    """
    tooltip_dict = {}
    try:
        with open(tooltips_file, 'r') as tfile:
            lines = list(tfile)
            for line in lines:
                # skip headers and blank lines
                if line.startswith('#') or line.startswith('\n'):
                    continue
                line = line.replace('"', '')
                splits = line.rstrip('\n').split('\t')
                try:
                    key = splits[0].strip()
                    value = splits[1].strip()
                    tooltip_dict[key] = value
                except ValueError:
                    logger.warning('Tooltip not parsed for line: {}'.format(line))
                    continue
        return tooltip_dict

    except FileNotFoundError:
        logger.error('params file not found!')


def init_logs():
    """
    Initialize logging code for CIUSuite 2. Logs debug information to file, warning and above to file and
    console output. NOTE: using the WARNING level really as INFO to avoid GUI builder logs from
    image creation being displayed on program start.
    :return: logger
    """
    mylogger = logging.getLogger('main')

    if mylogger.hasHandlers():
        mylogger.handlers.clear()
    # create file handler which logs even debug messages
    file_handler = RotatingFileHandler(log_file, maxBytes=1 * 1024 * 1024, backupCount=1)

    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(file_formatter)
    mylogger.addHandler(file_handler)

    # create console handler with a higher log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = CIU2ConsoleFormatter()
    console_handler.setFormatter(console_formatter)
    mylogger.addHandler(console_handler)

    mylogger.setLevel(logging.DEBUG)
    mylogger.propagate = False
    return mylogger


if __name__ == '__main__':


    try:
        with open(guardpath, 'w') as guardfile:
            guardfile.write('guarding')
    except IOError:
        pass

    try:
        with open(guardpath, 'w') as guardfile:
            guardfile.write('guarding')
    except IOError:
        pass
    logger = init_logs()
    multiprocessing.freeze_support()

    root = tk.Tk()
    root.withdraw()
    ciu_app = Fragariyo(root)
    logger.info('Starting FRagariyo')
    ciu_app.run()



    # closer handlers once finished
    for handler in logger.handlers:
        handler.close()
        logger.removeFilter(handler)

    # remove guard file
    try:
        os.remove(guardpath)
    except IOError:
        pass









