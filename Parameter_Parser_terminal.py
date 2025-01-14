"""
Author: Most code from DP
Date: Feb 19, 2020
"""

import combination
import os
import pickle


# parameter position in template file dictionary for quick reference. Key = parameter attribute name, value = position in splits
ppos = {'Analysis Num':0,
        'Analysis Name': 1,
        'seq': 2,
        'frag_chem':3,
        'iontypes': 4,
        'maxcharge': 5,
        'neutral_loss_bool':6,
        'disulfides': 7,
        'mods_array': 8,
        'r': 9,
        'uniprot_offset': 10,
        'ss_allowbroken': 11,
        'disulfides_ls': 12,
        'naturally_redcys':13,
        'mod_bool': 14,
        'noncys_mods': 15,
        'init_tol':16,
        'final_tol':17,
        'cal_bool':18
        }


TERMIONSECD = ['c', 'c-1', 'z+1', 'zdot', 'adot', 'b', 'y']
TERMIONSCID = ['a', 'b', 'y']


def parse_disulf_ls(disulf_str, uniprot_offset):
    """
    :param disulf_str: a list of strings of cysteine location involve in disulfide bonds (e.g ['15-18', ''])
    :return: A list of cystine locations that are involve in disulfide bonds (set), e.g. [{18, 15}]
    """

    spl = disulf_str.split(";")

    print(spl)
    ssls = []
    for ssbond in spl:
        if ssbond == '':
            continue
        else:
            # print(ssbond)
            enlace = ssbond.split("-")
            bondset = set()
            for num in enlace:
                #If statement helps to get rid off and empty set
                if num:
                    # print(f"num = {num}")
                    # print(f"num type = {type(num)}")
                    intnum = int(num) - uniprot_offset
                    bondset.add(intnum)
                    # print(f"bondset = {bondset}")
                # else:
                #     print("Not a number!")

        ssls.append(bondset)

    # print(f"ssls = {ssls}")
    return ssls


def str_to_ls(modstr):
    """
    :param modstr:
    :return:
    """
    # Parameteres for modification permutations for disulfide breakage

    modssplit = modstr.split(';')
    mods_ls = []
    for mod in modssplit:
        mods_ls.append((mod))
    return mods_ls

def parse_param_template_batch_multipass(param_file):
    """
    Read template csv file for all parameter information.
    :param param_file: (string) full system path to csv template file
    """

    params_dict = {}

    with open(param_file, 'r') as pfile:
        processed_analysis = 0
        for line in list(pfile):
            # print(line)

            # print(f"Processed: {processed_analysis}")
            if line.startswith('#'):
                continue

            splits = line.rstrip('\n').split(',')
            # print(splits)


            current_analysis = splits[ppos['Analysis Num']]
            # print(f"Current: {current_analysis}")
            if current_analysis != processed_analysis:
                params_dict[current_analysis] = []
            #Initilize params object
            params = Parameters()

            params.analysisName = splits[ppos['Analysis Name']]
            params.analysisNum = splits[ppos['Analysis Num']]
            params.seq = splits[ppos['seq']].strip()
            params.maxcharge = int(splits[ppos['maxcharge']])


            #Ion types
            iontypes_str = splits[ppos['iontypes']]
            iontypes_ls = []
            # print(iontypes_str)
            iontypes_strsplit = iontypes_str.split(';')


            fragmentation_type = splits[ppos["frag_chem"]]

            if fragmentation_type == "CID":
                for ion in TERMIONSCID:
                    for type in iontypes_strsplit:
                        if ion[0] == type:
                            iontypes_ls.append(ion)

            elif fragmentation_type == "ECD" or fragmentation_type == "ETD":
                for ion in TERMIONSECD:
                    for type in iontypes_strsplit:
                        if ion[0] == type:
                            iontypes_ls.append(ion)
            else:
                print("Only CID and ECD/ETD are currently supported!")

            print(iontypes_ls)
            params.iontypes = iontypes_ls

            #Neutral losses
            params.neutraloss_bool = parse_bool(splits[ppos['neutral_loss_bool']])


            #Parameteres for modification permutations
            # Parameteres for modification permutations for disulfide breakage
            params.mod_bool = parse_bool(splits[ppos['mod_bool']])
            modstr = splits[ppos['noncys_mods']]
            params.noncysmods = str_to_ls(modstr)


            #Disulfide_analysis
            params.disulfide_bool = parse_bool(splits[ppos['disulfides']])
            modstr = splits[ppos['mods_array']]
            params.arr = str_to_ls(modstr)

            natredcysstr = splits[ppos['naturally_redcys']]
            params.naturally_redcys = str_to_ls(natredcysstr)
            params.r = splits[ppos['r']]
            try:
                params.uniprot_offset = int(splits[ppos['uniprot_offset']])
            except ValueError:
                params.uniprot_offset = 0
                print(f"No uniprot offset given.")


            try:
                params.ss_allowbroken = int(splits[ppos['ss_allowbroken']])
            except ValueError:
                params.ss_allowbroken = 0
                print("Number of disulfides (within a fragment) allow to be broken not given.")



            #Parsing disulfides
            disulfides_str = splits[ppos['disulfides_ls']]
            try:
                params.disulfide_ls = parse_disulf_ls(disulfides_str, params.uniprot_offset)
            except ValueError:
                params.disulfide_ls = ''
                print("Number of disulfides ls not given.")


            params_dict[params.analysisNum].append(params)

            processed_analysis = current_analysis
            # print(f"after processed: {processed_analysis}")

            #Matching
            params.init_tol = int(splits[ppos['init_tol']])
            params.final_tol = int(splits[ppos['final_tol']])
            params.cal_bool = parse_bool(splits[ppos['cal_bool']])


    return params.seq, params_dict

def parse_bool(param_string):
    """
    Parse input strings to boolean
    :param param_string: string
    :return: bool
    """
    if param_string.lower() in ['t', 'true', 'yes', 'y']:
        return True
    elif param_string.lower() in ['f', 'false', 'no', 'n']:
        return False
    else:
        raise ValueError('Invalid boolean: {}'.format(param_string))

class Parameters(object):
    """
    Container to hold all parameters for searches to simplify method calls/etc
    """

    def __init__(self):
        """
        No parameters initialized immediately
        """
        self.params_dict = {}

        # ion prediction parameters
        self.analysisName = None
        self.analysisNum = None
        self.seq = None
        self.iontypes = None
        self.neutraloss_bool = None
        self.maxcharge = None
        self.mod_bool = None
        self.noncysmods = None

        # Disulfide Analysis
        self.arr = None
        self.r = None
        self.disulfide_bool = None
        self.uniprot_offset = None
        self.ss_allowbroken = None
        self.naturally_redcys = None

        #Matching
        self.init_tol = None
        self.final_tol = None
        self.cal_bool = None



    def set_params(self, params_dict):
        """
        Set a series of parameters given a dictionary of (parameter name, value) pairs
        :param params_dict: Dictionary, key=param name, value=param value
        :return: void
        """
        for name, value in params_dict.items():
            try:
                # only set the attribute if it is present in the object - otherwise, raise attribute error
                self.__getattribute__(name)
                self.__setattr__(name, value)
            except AttributeError:
                # no such parameter
                print('No parameter name for param: ' + name)
                continue
        self.update_dict()

    def combodict_calc(self):

        try:
            combo_dict = combination.batch_combos(self.arr, int(self.r))
        except ValueError:
            combo_dict = {}

        return combo_dict

    def update_dict(self):
        """
        Build (or rebuild) a dictionary of all attributes contained in this object
        :return: void
        """
        for field in vars(self):
            value = self.__getattribute__(field)
            self.params_dict[field] = value

    def __str__(self):
        """
        string
        :return: string
        """
        return '<Params> protein {}'.format(self.analysisName)
    __repr__ = __str__

class ExpIon:
    """
    Container for experimental data, corresponding to a peak cluster with monoisotopic peak information.
    Designed to handle several different input types with various levels of information, but most easily
    used with outputs from IMTBX/Grppr
    terminalFragmentor
    """
    def __init__(self, init_data, mmass_bool=None):
        """
        Initialize a new ExpIon container
        :param init_data = list of input data to generate the cluster. Must contain the following entries
            0: mz_mono: monoisotopic peak m/z
            1: dt_mono: monoisotopic peak drift time (bin)
            2: pkht_mono: monoisotopic peak height
            3: pkar_mono: monoisotopic peak area
            4: mz_toppk: most abundant peak m/z
            5: dt_toppk: most abundant peak drift time (bin)
            6: pkht_toppk: most abundant peak m/z height
            7: pkar_toppk: most abundant peak m/z area
            8: num_pks: number of peaks in cluster
            9: idx_top: index of most abundant peak in cluster
            10: charge: cluster charge
            11: mz_avg: cluster average m/z
            12: pkht_cluster: cluster height (summed)
            13: pkar_cluster: cluster area (summed)
            14: correlation: averagine correlation score (if doing averagine modeling, added in v2.4.0.0)
            15: noise: noise level used in determining SNR (added in v2.5.0.0)
        """
        if mmass_bool is None or mmass_bool is False:
            self.mz_mono = round(init_data[0], 3)
            self.dt_mono = init_data[1]
            self.pkht_mono = init_data[2]
            self.pkar_mono = init_data[3]
            self.mz_toppk = init_data[4]
            self.dt_toppk = init_data[5]
            self.pkht_toppk = init_data[6]
            self.pkar_toppk = init_data[7]
            self.num_pks = init_data[8]
            self.idx_top = init_data[9]
            self.charge = init_data[10]
            self.mz_avg = init_data[11]
            self.pkht_cluster = init_data[12]
            self.pkar_cluster = init_data[13]
            self.correlation = init_data[14]
            self.noise = init_data[15]
            self.data_list = init_data
        else:
            # mMass data, which does not have many of these fields
            self.mz_mono = round(init_data[0], 4)
            self.pkht_cluster = init_data[1]
            self.charge = init_data[4]

            # mMass specific params
            self.sig_noise = init_data[3]
            self.fwhm = init_data[5]
            self.resolution = init_data[6]
            self.data_list = init_data
            # allow easy use of Dmitry-based proc methods
            self.pkar_cluster = self.pkht_cluster
            self.isotopeenvelope = init_data[7]
        self.cal_mz_mono = None

    def __lt__(self, other):
        return self.mz_mono < other.mz_mono

    def __eq__(self, other):
        return self.mz_mono == other.mz_mono

    def __hash__(self):
        # print(hash(str(self)))
        return hash(self.mz_mono)

    def __round__(self, n=2):

        return round(self.mz_mono,n)



    def __str__(self):
        """
        string representation
        :return: string
        """
        return f'<Exp> mz: {self.mz_mono}, z: {self.charge}, int: {self.pkht_cluster}'
    __repr__ = __str__

def load_hits_file(filename):
    """
    Load a saved (pickled) .hits file and return the stored list of FragmtSites
    :param filename: full path to file to load
    :return: list of FragmtSite containers with hits information
    :rtype: list[FragmentSite]
    """
    with open(filename, 'rb') as loadfile:
        sitelist = pickle.load(loadfile)
    return sitelist

def unified_exp_parser(input_file):
    """
    Single entry point for all experimental input types. Determines file type and parses it if possible,
    raises TypeError if the file is not of a known type.
    :param input_file: full system path to input file to parse
    :return: list of ExpIons with peaklist information, base filename without extension
    """
    filetype = input_file.split('.')[-1]
    print(f"Filetype {filetype}")
    short_filename = os.path.basename(input_file).rstrip('.' + filetype)
    if filetype == 'isotopes':
        # IMTBX/Grppr file
        exp_ions = parse_peaks_grppr(input_file)
    elif filetype == 'txt':
        # mMass file
        exp_ions = parse_mmass_peaklist(input_file)
    elif filetype == 'csv':
        exp_ions = csvfile_parser(input_file)
    elif filetype == 'unmatched':
        exp_ions = load_hits_file(input_file)
    else:
        # unknown filetype
        raise TypeError(filetype)
    return exp_ions, short_filename, filetype

def parse_xtract_peaklist(input_file):
    """
    Parse csv file generated by copying the monoisotopic peak list from Thermo's Xtract output
    into excel and saving as csv.
    :param input_file: csv file with Xtract monoisotopic output from xtract
    :return: list of ExpIons
    """
    peak_list = []
    with open(input_file, 'r') as peaksfile:
        for line in list(peaksfile):
            splits = line.split(',')
            try:
                mz = float(splits[0])
                intensity = float(splits[1])
                # current file doesn't have charge info, which is needed. Will try to find Xtract output with charge info...
            except ValueError:
                # header or other non peak line, skip
                continue




def parse_peaks_grppr(input_isotopes_file):
    """
    Parse .isotopes files from Dmitry's tool. Returns a list of ExpCluster objects containing peak info
    :param input_isotopes_file: (string) full system path to .isotopes file to parse for peaks
    :return: list of ExpIon objects containing parsed information
    :rtype: list[ExpIon]
    """
    peak_list = []
    # line_count=0

    # read the file
    with open(input_isotopes_file, 'r') as input_file:
        lines = list(input_file)
        for line in lines:
            # line_count+=1
            # file is tab delimited
            splits = line.split('\t')

            # ignore headers by ensuring values passed are floats
            arg_list = []
            try:
                for value in splits:
                    myval = float(value)
                    arg_list.append(myval)
            except ValueError:
                # this is a header line, continue to the next line
                continue

            # Don't make ExpIon objects from peaks which have a cluster area of less than 2500
            pkar_cluster = arg_list[13]
            # print(pkar_cluster)

            if pkar_cluster > 2500:
                # create a new exp cluster object using the input data from this line and append it to the peak_list
                peak_list.append(ExpIon(arg_list))
            else:
                continue
    # print(line_count)
    return peak_list


def get_filename(isotopes_file_input):
    """
    Returns the sample name for a Dmitry tool-style .isotopes file for naming outputs
    :param isotopes_file_input: .isotopes file full system path
    :return: original sample name
    """
    filesplits = isotopes_file_input.split('/')
    filename = filesplits[len(filesplits) - 1]
    filename = filename.rstrip('.isotopes')
    return filename


def parse_mmass_peaklist(input_file):
    """
    Parse an mMass peaklist and return a list of cluster objects for peak matching
    :param input_file: the .txt file in mMass format to read (comma separated)
    :return: list of cluster objects with all fields mMass can provide
    """
    peak_list = []

    with open(input_file, 'r') as peaks_file:
        lines = list(peaks_file)
        for line in lines:
            line = line.rstrip('\n')
            splits = line.split(',')
            # splits = line.split('\t')
            arg_list = []
            try:
                for value in splits:
                    if value != '':
                        myval = float(value)
                        arg_list.append(myval)
                    else:
                        arg_list.append(value)
            except ValueError:
                # this is a header line, continue to the next line
                continue

            # create a new exp cluster object using the input data from this line and append it to the peak_list
            peak = ExpIon(arg_list, mmass_bool=True)
            peak_list.append(peak)
    return peak_list


def parse_single_exp_mmass_zcheck(input_file):
    """
    Open an mmass file with filechooser and return its peaklist, filename. Reorders input data
    for charge checking.
    :param input_file: full system path to input file to read
    :return: list of ExpIon objects with mmass inits
    """
    peak_list = []
    with open(input_file, 'r') as peaks_file:
        lines = list(peaks_file)
        for line in lines:
            line = line.rstrip('\n')
            splits = line.split(',')
            arg_list = []
            try:
                for value in splits:
                    if value != '':
                        myval = float(value)
                        arg_list.append(myval)
                    else:
                        arg_list.append(value)
            except ValueError:
                # this is a header line, continue to the next line
                continue
            # re-organize arg list to match expected input format

            ordered_list = [arg_list[2], arg_list[5], arg_list[3]]

            # create a new exp cluster object using the input data from this line and append it to the peak_list
            peak = ExpIon(ordered_list, True)
            peak_list.append(peak)

    return peak_list

def csvfile_parser(csv_file):
    """
    :param csv_file: A csv file with int and mono_mz columns
    :return: A list of ExpIon objects
    """

    peak_list = []
    with open(csv_file, 'r') as peaks_file:
        lines = list(peaks_file)
        # print(lines)
        for line in lines:
            if line.startswith("#"):
                continue
            line = line.rstrip('\n')
            splits = line.split(',')
            # print(splits)

            arg_list = []
            for value in splits:
                if value:
                    arg_list.append(value)

            # print(arg_list)
            # re-organize arg list to match expected input format

            mz = arg_list[0]
            charge = arg_list[1]
            intensity = arg_list[2]
            # If using the in-house peak picking
            try:
                isotopicenvelope_centroid = arg_list[4]
            except IndexError:
                isotopicenvelope_centroid = ""


            if str(charge).isnumeric():

                try:
                    ordered_list = [float(mz), float(intensity), "", "", int(charge), "", "", isotopicenvelope_centroid]
                except ValueError:
                    print(" One of the peak attributes is the wrong type. Check that no value for each peak is missing.")

                # create a new exp cluster object using the input data from this line and append it to the peak_list
                peak = ExpIon(ordered_list, True)

                peak_list.append(peak)

            else:
                print(f"The charge values are no numeric! Please check the column format of the experimental ion CSV")
                break
    print(f"peak_list = {peak_list}")
    return peak_list



if __name__ == "__main__":
    pass