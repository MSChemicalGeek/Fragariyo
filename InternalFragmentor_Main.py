"""
Author: Carolina Rojas Ramirez
Date: 07/11/2018
In silico fragmentation of proteins using mass.fast_mass2 from pyteomics
"""
import time
import combination
from pyteomics import mass
from pyteomics import parser
from tkinter import filedialog
from tkinter import simpledialog
import tkinter as tk
import os
import pickle
import Parameter_Parser
import re
from tkinter import messagebox
from Modifications import mods_repo
from PeakMatch import matchmaker2_multipass
from Parameter_Parser import expion_parser
from Parameter_Parser import isotope_xtractor
import terminalFragmentor_Main as tfm


# Setting up tkinter
root = tk.Tk()
root.withdraw()


class interfrag:
    """
    Container object for internal fragments and their properties
    :param mz_mono: int, monoisotope m/z
    :param charge: int, ion charge
    :param ion type: str
    :param sequence: str, internal fragment
    :param cys_num: int, cysteines available for modification
    :param mono_neutral: int, ion neutral mass
    :param mods: ls, mods if any
    :param cysloc: set, locations of cysteines
    :param ss_count: int, amount of disulfide bonds in the internal fragment
    :param ifragstart: Index where the fragment starts in the protein sequence
    :param ifragend: Index where the fragment ends in the protein sequence
    :param reverse: bool, is the sequence a reverse version.
    """

    def __init__(self, mz_mono, charge, ion_type, sequence, cys_num, mono_neutral, mods, cysloc, ss_count, cysmods, ifragstart, ifragend, reverse):
        self.mz_mono = round(mz_mono, 8)
        self.charge = charge
        self.ion_type = ion_type
        self.sequence = sequence
        self.cys_num = cys_num
        self.mono_neutral = round(mono_neutral,8)
        self.mods = mods
        self.cysmods = cysmods
        self.cysloc = cysloc
        self.ss_count = ss_count

        self.ifragstart = ifragstart
        self.ifragend = ifragend
        self.reverse = reverse


    #Attribute to be able to compare to other objects
    def __eq__(self, other):
        return self.mono_neutral == other.mono_neutral

    #Attribute to represent objects in dictionaries
    def __hash__(self):
        return hash(self.mono_neutral)

    # Attribute to represent objects in stdout
    def __str__(self):
        return f"{self.mono_neutral}\t{self.sequence}\t{self.charge}\t{self.mz_mono}\t{self.mods}\t{self.ion_type}\t{self.cysloc}\t{self.ss_count}\t{self.cys_num}\t{self.cysmods}\t{self.ifragstart}\t{self.ifragend}\t{self.reverse}"

    # Attribute to represent objects in other objects
    __repr__ = __str__

def residue_content(frag, seqstart):
    """
    Function to create fragment site objects
    :param frag: str, sequence
    :param protein_sequence: str, full protein sequence
    :param terminal: str, 'N' or 'C'
    :param ss_ls: ls, a list of the disulfide bonds (set) in the protein
    :param intrabroken_ss: int, how many disulfide bonds (kept intact inside the fragment site) are meant to be broken
    :param natredcys: ls, a list of the location (int) of cysteines  that do not participate in disulfide bonds
    :return: fragment site object
    """

    residue_dict = {"A":[],"R":[],"N":[],"D":[],"C":[],"E":[],"Q":[],"G":[],"H":[],"I":[],"L":[],"K":[],"M":[],"F":[],"P":[],"S":[],"T":[],"W":[],"Y":[],"V":[]}


    # # Use regex to locate initial residue location in the whole protein sequence
    # interindex = re.compile(frag, re.I)
    # mo = interindex.search(protein_sequence)
    #
    # # Add one to "translate from python indexing to residue number
    # seqstart = mo.start() + 1
    #
    # seqend = mo.end()

    # Analysis per Residue
    for residue in residue_dict:
        amino_acid = re.compile(fr'{residue}', re.I)
        # residue locations
        resloc = set()
        iterr = amino_acid.finditer(frag)
        for res in iterr:
            # print(cys.group(), cys.start())
            resloc.add(res.start() + seqstart)
        residue_dict[residue] = resloc
        # print(cysloc)

    # print(f"This is a {terminal}-term site = {frag}_{iseqstart}-{iseqend}. Amino Acid composition = {residue_dict}")

    return residue_dict

def varmod_combos(dict_of_varmods):
    """
    :param dict_of_varmods: A dictionary of the variable modifications to analyze
    :return: A dictionary of the variable modifications to analyze with the number of modifications
    """
    varmodls = []
    for varmod in dict_of_varmods:
        for varmod_num in dict_of_varmods[varmod]:
            varmodls.append(f"{varmod}_{varmod_num}")

    # print(f"varmodls = {varmodls}")
    totvarmods = len(dict_of_varmods)
    varmodcombos = combination.rSubset(varmodls, totvarmods)
    # print(varmodcombos)
    return varmodcombos

def varmods_processing(residue_survey, modification, amino, var_mods_dict):
    """
    Adding variable modifications
    :param FragmentSiteObj: Fragment Site container
    :param modification: Modification object
    :param amino: Residue to be modified
    :param var_mods_dict: Dictionary with the modification and their amount
    :return: A dictionary of variable modification tailored to the current fragment site
    """
    # Variable mod part
    # print(f"before adding varmods = {mods_repo[modification].current_num}")
    varmods_num = mods_repo[modification].current_num

    if len(residue_survey[amino]) >= mods_repo[modification].max_num:

        for num_mod in range(1, mods_repo[modification].max_num + 1):

            # print(f"Starting at each num_mod = {mods_repo[modification].current_num}")

            # print(f"put {num_mod} times = m/z {(mods_repo[modification].mass) * num_mod}")
            varmods_num += num_mod

            # print(f"mod current num = {varmods_num}")
            if varmods_num > mods_repo[modification].max_num:
                continue
            else:
                # var_mods_dict[f"{modification}[{num_mod}]"] = (mods_repo[modification].mass) * num_mod

                var_mods_dict[f"{modification}"].append(num_mod)
            varmods_num -= num_mod

        # mods.remove(modification)

    else:
        # print(f"The places for var_mods {FragmentSiteObj.resi_dict[amino]} is less than the allowed mod max")

        # for num_mod in range(1, len(FragmentSiteObj.resi_dict[amino]) + 1):
        for num_mod in range(1, len(residue_survey[amino]) + 1):
            # print(f"Available aminos = {len(FragmentSiteObj.resi_dict[amino])}")

            # print(f"put {num_mod} times = m/z {(mods_repo[modification].mass) * num_mod}")
            varmods_num += num_mod
            # print(f"mod current num = {varmods_num}")
            if varmods_num > len(residue_survey[amino]):
                continue
            else:
                # var_mods_dict[f"{modification}[{num_mod}]"] = (mods_repo[
                #                                                   modification].mass) * num_mod
                var_mods_dict[f"{modification}"].append(num_mod)
            varmods_num  -= num_mod
        # mods.remove(modification)

    mods_repo[modification].current_num = 0
    return var_mods_dict

def modificator(residue_dict, mod_ls, charge, var_mods_dict, mz_mono, neutral_mono, mods_repo):
    """
    Places modifications on a Fragment Site object using all the information necessary
    :param FragmentSiteObj: container for the Fragment Site info
    :param mod_ls: A list of modifications to be considered
    :param charge: Charge state of the theoretical ion to be calculated
    :param var_mods_dict: Dictionary of the modification to be palced with their respective amounts
    :param mz_mono: float, m/z value
    :param neutral_mono: float, deconvoluted mass
    :return: mz_mono, neutral_mono, mods, var_mods_dict for a theoretical ion
    """
    # If modifications need to be considered

    # mods = mod_ls.copy()

    #Intitialized list and masses to be calculated for the theoretical ion
    mods = []
    accum_mz_fx = 0
    accum_neutral_fx = 0

    for modification in mod_ls:
        mods_repo[modification].current_num = 0
        # print(f"the modObj is {mods_repo[modification]} with type {type(mods_repo[modification])}")
        # print(f"the FragSite has the an amino acid composition of {FragmentSiteObj.resi_dict}")
        # print(f"the modification {mods_repo[modification]} can be placed on {mods_repo[modification].target_aas}")


        if mods_repo[modification].target_aas:
            for amino in mods_repo[modification].target_aas:
                # print(f"In the  fragment, {amino} is at positions {FragmentSiteObj.resi_dict[amino]}")

                if residue_dict[amino]:
                    if mods_repo[modification].fixed:
                        for res in mods_repo[modification].fixed:
                            if res in residue_dict[amino]:
                                # print(
                                #     f"The modification is fixed at  {mods_repo[modification].fixed} which is in this fragment!")
                                accum_mz_fx += (mods_repo[modification].mass) / charge
                                accum_neutral_fx += mods_repo[modification].mass
                                mods.append(modification)
                            # else:
                            #     # print(
                            #     #     f"The modification is fixed at  {mods_repo[modification].fixed} which is NOT in this fragment!")
                            #     mods.remove(modification)


                    elif mods_repo[modification].terminal_flag:

                        # if mods_repo[modification].terminal_flag == FragmentSiteObj.term:
                        #     # mods_repo[modification].current_num = 0
                        #     # print(
                        #     #     f"The modification is variable with a max number of {mods_repo[modification].max_num} and it can also be terminal at the {mods_repo[modification].terminal_flag} term, so here! ")
                        #
                        #     mods_repo[modification].current_num += 1
                        #     accum_mz_fx += (mods_repo[modification].mass) / charge
                        #     accum_neutral_fx += mods_repo[modification].mass
                        #     mods.append(f"{mods_repo[modification].terminal_flag}-terminal_{mods_repo[modification].name}")
                        #
                        #     var_mods_dict[f"{modification}"] = []
                        #     var_mods_dict = varmods_processing(FragmentSiteObj, modification, amino, var_mods_dict)
                        #
                        # else:
                        #     # print(
                        #     #     f"The modification is variable with a max number of {mods_repo[modification].max_num} and it can also be terminal at the {mods_repo[modification].terminal_flag} term, not this term")
                        #
                        #     var_mods_dict[f"{modification}"] = []
                        #     var_mods_dict = varmods_processing(FragmentSiteObj, modification, amino,
                        #                                              var_mods_dict)
                        #     # mods.remove(modification)
                        print("terminal flags make no sense for internal fragments!")

                    else:
                        # print(f"The modification is variable with a max number of {mods_repo[modification].max_num}")
                        var_mods_dict[f"{modification}"] = []
                        var_mods_dict = varmods_processing(residue_dict, modification, amino, var_mods_dict)


                else:
                    # print(f"There is no {amino} to put it at!")
                    # mods.remove(modification)
                    pass

        else:
            # # print(f"{modification} is only terminal!")
            # if mods_repo[modification].terminal_flag == FragmentSiteObj.term:
            #     mods.append(
            #         f"{mods_repo[modification].terminal_flag}-terminal_{mods_repo[modification].name}")
            #     mods_repo[modification].current_num += 1
            #     accum_mz_fx += (mods_repo[modification].mass) / charge
            #     accum_neutral_fx += mods_repo[modification].mass
            #
            #
            # else:
            #     # mods.remove(modification)
            #     # print("This terminal modification does not belong in this terminal!")
                pass

    mz_mono += accum_mz_fx
    neutral_mono += accum_neutral_fx

    return mz_mono, neutral_mono, mods, var_mods_dict

def mass_calc(seq, types, maxcharge, dictionary, ss_bonds=None, cysmod_dict = None, cys_num=None, mods_cys = None, cysloc = None, sscount = None, iseqstart=None, iseqend=None, reverse_flag=None, modobj =None, libraryofmods=None):
    """
    Function to calculate masses of specific internal fragments with specific disulfides configurations
    :param seq: str, interfragment sequence
    :param types: ls, ions types
    :param mincharge: Minimun charge
    :param maxcharge: Maximun charge
    :param ss_bonds: bool, disulfide bonds present
    :param cysmod_dict: Dictionary with number of disulfides available for modification (keys), possible modifications (values)
    :param cys_num: Cysteines available for modification
    param cysloc: set, location of cysteines
    :param ss_count: int, number of disulfide bonds
    :param iseqstart:  Index where the internal fragment starts at in the protein sequence
    :param iseqend:  Index where the internal fragment ends at in the protein sequence
    :param reverse_flag:  bool, is this a reverse sequence?
    :param modobj = ls of noncys modification objects
    """

    # print(f"{FragmentSiteObj}")
    # print(f"{FragmentSiteObj.seq}")
    # For each ion type

    residue_census = residue_content(seq,iseqstart)


    for ion_type in types:
        # For each charge
        for charge in range(maxcharge, maxcharge + 1):
            # print(f"charge = {charge}")
            #calculate m/z and neutral
            mz_mono_org = mass.fast_mass2(seq, ion_type=ion_type, charge=charge)
            neutral_mono_org = mass.fast_mass2(seq, ion_type=ion_type, charge=0)

            # print(f"\nOriginal mz = {mz_mono_org}")
            # print(f"Original mass = {neutral_mono_org}")

            # Create a new m/z and mass variable to contain the original mz and neutral mass
            mz_mono = mz_mono_org
            neutral_mono = neutral_mono_org

            # To contain the end results of going thru all the conditions below
            mods = []
            var_mods_dict = {}
            # If disulfides are considered
            if ss_bonds:
                # print(f"ss_bonds = {ss_bonds}")
                # print(f"cys_num = {cys_num}")
                # If the sequence has not free cysteines for modification
                if cys_num == 0:
                    # if besides considering disulfides, possible modifications need to be considered

                    if modobj:
                        mz_mono, neutral_mono, mods, var_mods_dict = modificator(residue_census, modobj, charge,
                                                                                 var_mods_dict, mz_mono, neutral_mono,libraryofmods)

                        # Removing the hydrogens due to the disulfide bonds
                        neutral_mono = neutral_mono + (sscount * (-1.0078 * 2))
                        mz_mono = mz_mono + (sscount * ((-1.0078 / charge) * 2))


                    else:
                        # Do not modified cysteines
                        mods = []

                        # Removing the hydrogens due to the disulfide bonds
                        neutral_mono = neutral_mono + (sscount * (-1.0078 * 2))
                        mz_mono = mz_mono + (sscount * ((-1.0078 / charge) * 2))

                # If there are disulfides to be modified because their disulfide partner is not in this fragment
                else:
                    # If modifications need to be considered
                    if modobj:
                        mz_mono, neutral_mono, mods, var_mods_dict = modificator(residue_census, modobj, charge,
                                                                                 var_mods_dict, mz_mono, neutral_mono,libraryofmods)


                        # Modify neutral mass with the modifications possible
                        # print(f"cysmod_dict = {cysmod_dict}/cys_num = {cys_num}/cysmod_dict[cys_num] = {cysmod_dict[cys_num]}")
                        neutral_mono = neutral_mono + cysmod_dict[cys_num][mods_cys] + (sscount * (-1.0078 * 2))
                        mz_mono = mz_mono + (cysmod_dict[cys_num][mods_cys]/ charge) + (sscount * ((-1.0078 / charge) * 2))

                    else:
                        # If the internal fragment contains cysteines for modifications


                        # Modify neutral mass with the modifications possible
                        neutral_mono = neutral_mono + cysmod_dict[cys_num][mods_cys] + (sscount * (-1.0078 * 2))
                        mz_mono = mz_mono + (cysmod_dict[cys_num][mods_cys] / charge) + (
                                sscount * ((-1.0078 / charge) * 2))


            # If there is no disulfides to be considered, but only modifications
            elif modobj:
                print(f"modobj = {modobj}")
                mz_mono, neutral_mono, mods, var_mods_dict = modificator(residue_census, modobj, charge,
                                                                         var_mods_dict, mz_mono, neutral_mono,libraryofmods)




            # If there were variable modifications to considered
            if var_mods_dict:
                # print(var_mods_dict)
                # print(f"There are {len(var_mods_dict)} var mods!")
                possible_varmodscombos = varmod_combos(var_mods_dict)
                for combo in possible_varmodscombos:
                    # print(combo)
                    combols = list(combo)
                    mass_combo = 0
                    for mod in combols:
                        modsplits = mod.split("_")
                        # mass_combo += mods_repo[mod[:-2]].mass * int(mod[-1])
                        mass_combo += mods_repo[modsplits[0]].mass * int(modsplits[1])
                    # print(f"Mass Combo {mass_combo}")
                    mods.append(combols)
                    mz_mono += mass_combo / charge
                    neutral_mono += mass_combo
                    modscopy = mods.copy()
                    # Create theoretical ion!
                    interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono, modscopy, cysloc,
                                             sscount, mods_cys, iseqstart, iseqend, reverse_flag)

                    # print(f"Interobj {interobj}")

                    # Add theoretical ion to its site in the main dictionary of fragment sites
                    dictionary[interobj.mono_neutral] = interobj

                    mz_mono -= mass_combo / charge
                    neutral_mono -= mass_combo

                    # Remove variable mod after being done with it
                    mods.remove(combols)

            else:

                # Creating theoretical ion
                modscopy = mods.copy()
                interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono, modscopy, cysloc,
                                         sscount, mods_cys, iseqstart, iseqend, reverse_flag)

                dictionary[interobj.mono_neutral] = interobj
                # print(f"Interobj {interobj}")






def ifragments(analysis_name, sequence, types=('b', 'y'), maxcharge=1, maxstart = 4, maxlength=10,
               modbool = None, combo_dict = None, cystine = None, uniprot_offset=None, allow_ssbroken = None,
               reverse_seq=None, foo_ls=None, reduced_cys=None, modsdictio=None, fragmentmode=None):
    """
    :param analysis_name: Str, the name of the current analysis or pass
    :param sequence: a string of the sequence desired to fragment.
    :param types: the ion types like b or y. To pass more than one type of ion put them in parenthesis
    :param maxcharge: the maximun charge for the fragments.
    :param maxstart: int, the smallest internal fragment size (number of residues)
    :param maxlength: int, the largest internal fragment size (number of residues)
    :param modbool: Modifications
    :param maxmods:
    :param combo_dict: A dictionary with all the possible combinations of disulfide bon mods
    :param cystine: Bool, considered disulfie bonds
    :param uniprot_offset: The amount of residues the experimental sequence is offset from the its uniprot sequence and be able to properly locate disulfide bonds
    Important in order to match disulfide bond location
    :param allow_ssbroken: How many disulfide bonds within an internal fragment
    :param reverse_seq: If it is true it means that the protein sequence are reverse
    :param foo_ls: Disulfide ls
    :return: A dictionary with neutral masses as keys and interobj as values
    """

    print('---------Fragmentation starts--------')
    all_start = time.time()


    global_dict = {}
    counter = 0
    #creating internal sequences
    interfrag_list = []
    isequence = sequence[0:-1]

    for i in range(1, len(isequence)):

        #Proline effect
        if isequence[i] == "P" and fragmentmode != "CID":
            continue
        else:

            for n in range(maxstart, maxlength+1):
                interseq = isequence[i:i+n]

                #Eliminating redundancies
                if interseq not in interfrag_list:
                    interfrag_list.append(interseq)
                else:
                    continue

    # print(interfrag_list)

    #For sequences with desired length
    for iseq in interfrag_list:
        if len(iseq) > maxstart:
            # print(f"\n {iseq}")
            #Use regex to locate initial residue location in the whole protein sequence
            interindex = re.compile(iseq, re.I)
            mo = interindex.search(sequence)
            # print(mo.start()+1)

            #Add one to "translate from python indexing to residue number
            #Add uniport offset to match uniport residue numbers and properly id disulfuide bonds
            iseqstart = mo.start() + 1 + uniprot_offset
            # print(iseqstart)
            # print(modbool)
            # print((iseqstart/isequence)*100)
            iseqend = mo.end() + 1 + uniprot_offset
            # print(iseqindex)

            # print(f" Before modbool loop Modbool = {modbool}, lenght = {len(modbool)}")

            if cystine:
                # print("Considering Cysteines!")

                unboundcys = 0

                # disulfide analysis
                Cysteines = re.compile(r'C', re.I)
                #cysteine locations
                cysloc = set()
                iterc = Cysteines.finditer(iseq)
                for cys in iterc:
                    # print(cys.group(), cys.start())
                    cysloc.add(cys.start() + iseqstart)


                # The number of cysteines in the internal fragment to be considered for disulfides and for modification
                cysloclen = len(cysloc)

                # print(f"cysloc = {cysloc}")
                # print(f"foo_ls = {foo_ls}")

                disulfide_counter = 0
                # foos_ls is the list of disulfide bond pairs from uniprot
                for pair in foo_ls:
                    if pair.issubset(cysloc):
                        # print(pair)
                        disulfide_counter += 1

                # print(disulfide_counter)



                # naturally_reducedcys_ls = [58, 476]
                # Carefull taht is a list of str
                # print(f"reduced_cys = {reduced_cys}")

                naturally_reducedcys_ls = reduced_cys

                for cys in cysloc:
                    if int(cys) in naturally_reducedcys_ls:
                        # print(cys)
                        cysloclen -= 1

                # Allow disulfides inside the internal fragment to be broken
                if allow_ssbroken and cysloclen>1:
                    # Unbound cysteines are cysteines that are free, don't have their didulfide bond partner in the current fragment, or
                    # cysteines that were reduced (allow_ssbroken)
                    unboundcys = cysloclen - (disulfide_counter * 2) + (allow_ssbroken * 2)
                else:
                    # Considering disulfides that reached outside the internal fragment to be broken
                    # Cysteines available for modification
                    unboundcys = cysloclen - (disulfide_counter * 2)


                # print(f"disulfide count = {disulfide_counter}")
                # Id the amount of unbound cys is negative it means there is no cys for modification
                if unboundcys < 0:
                    unboundcys = 0

                # print(f" before loop, unboundcys = {unboundcys}")
                if len(modbool) >=1:
                    # print("Considering Cysteines and modifications!")

                    if unboundcys > 0:

                        if combo_dict[unboundcys]:
                            for cysmod in combo_dict[unboundcys]:
                                mass_calc(iseq, types, maxcharge, global_dict, ss_bonds=True,
                                          cysmod_dict=combo_dict, mods_cys=cysmod,
                                          cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter,
                                          iseqstart=iseqstart,
                                          iseqend=iseqend, reverse_flag=reverse_seq, modobj=modbool,libraryofmods=modsdictio)
                                counter += 1
                        else:
                            print("Did you remember to put the cys modifications under mods_array column in Parameters file?")




                    else:

                        mass_calc(iseq, types, maxcharge, global_dict, ss_bonds=True,
                                      cysmod_dict=combo_dict, mods_cys=None,
                                      cys_num=0, cysloc=cysloc, sscount=disulfide_counter, iseqstart=iseqstart,
                                      iseqend=iseqend, reverse_flag=reverse_seq, modobj=modbool,libraryofmods=modsdictio)
                        counter += 1
                else:
                    pass

                # print("Considering Cysteines and no modifications!")

                if unboundcys > 0:

                    for cysmod in combo_dict[unboundcys]:
                        mass_calc(iseq, types, maxcharge, global_dict, ss_bonds=True, cysmod_dict=combo_dict,
                                  mods_cys=cysmod,
                                  cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter, iseqstart=iseqstart,
                                  iseqend=iseqend, reverse_flag=reverse_seq,libraryofmods=modsdictio)
                        counter += 1

                else:

                    mass_calc(iseq, types, maxcharge, global_dict, ss_bonds=True,
                              cysmod_dict=combo_dict, mods_cys=None,
                              cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter, iseqstart=iseqstart,
                              iseqend=iseqend, reverse_flag=reverse_seq,libraryofmods=modsdictio)
                    counter += 1

            else:
                # print(f"modsdictio = {modsdictio}")
                # print(f"modbool = {modbool}")
                if modbool:
                    # print("Considering modifications!")
                    mass_calc(iseq, types, maxcharge, global_dict, ss_bonds=0, cysmod_dict=0,
                              cys_num=0, cysloc=0, sscount=0, iseqstart=iseqstart, iseqend=iseqend, reverse_flag=reverse_seq,
                              modobj=modbool,libraryofmods=modsdictio)
                    counter += 1
                else:
                    # print("Regular primary sequence!")
                    mass_calc(iseq, types, maxcharge, global_dict, ss_bonds=0, cysmod_dict=0,
                              cys_num=0, cysloc=0, sscount=0, iseqstart=iseqstart, iseqend=iseqend,
                              reverse_flag=reverse_seq,libraryofmods=None)
                    counter += 1

    all_end = time.time() - all_start
    print('Total prediction time: {}'.format(round(all_end)))
    print(f"{counter} Internal fragments were produced")
    # print(global_dict)
    return (analysis_name, global_dict)




def fragments_to_picklefile(file_title, fragments_dict):
    """
    Save the fragment function as a binary file
    :return:
    """
    print(f"file_title = {file_title}")
    outname = file_title.replace(".csv",".theofrags")
    with open(outname, 'wb') as picklefile:
        pickle.dump(fragments_dict, picklefile)

def fragments_to_FRAGSfile(file_title, fragments):
    """
    Save the fragment function output as a .txt file (recommended only for small proteins)
    :param file_title: Str, title of the file
    :param fragments: Str, a string with the fragments
    :return: void
    """
    outname = file_title.replace(".csv", ".txt.frags")
    output = open(outname, 'w')
    output.write(fragments)
    output.close()

def unmatched_expions_outfile_writer(exp_ion_list, output_filename, csv = False):
    """

    :param exp_ion_list:
    :param output_filename:
    :param csv:
    :return:
    """
    # Test to get the unid exp-ions out and somehow, even if it manually to determine what are the neutral/radicall losses

    if csv:
        charge_dict = {}
        for ion in exp_ion_list:
            if ion.charge in charge_dict:
                charge_dict[ion.charge].append(ion)
            else:
                charge_dict[ion.charge] = [ion]

        # print(charge_dict)

        outstr = ""
        output = open(output_filename.strip(".csv") + "_Unmatched" + '.csv', 'w')
        for ion in exp_ion_list:
            outstr += f"{ion.exp_mz},{ion.charge},{ion.pkht_cluster}\n"

        output.write(outstr)
        output.close()
    else:
        tfm.print_unmatched(exp_ion_list, output_filename)
        outnamefile = output_filename.strip(".csv")
        tfm.save_fragments(exp_ion_list, outnamefile + ".unmatched")


def intern_multipass_batchparser(batch_file):

    masterlist = []
    workingdir = ''
    with open(batch_file, 'r') as pfile:

        indx = 0
        for line in list(pfile):

            if indx == 0:
                splits = line.rstrip('\n').split(',')
                workingdir = splits[0]
                print(workingdir)
                indx += 1

            else:
                if line.startswith('#'):
                    indx += 1
                else:
                    line_ls = []
                    print(line)
                    splits = line.rstrip('\n').split(',')
                    expionfile = splits[0]
                    expionpath = os.path.join(workingdir,expionfile)
                    line_ls.append(expionpath)
                    ionsfile = splits[1]
                    if ionsfile == "":
                        line_ls.append("")
                    else:
                        ionpath = os.path.join(workingdir, ionsfile)
                        line_ls.append(ionpath)
                    paramsfile = splits[2]
                    if paramsfile == "":
                        line_ls.append("")
                    else:
                        parampath = os.path.join(workingdir, paramsfile)
                        line_ls.append(parampath)

                    masterlist.append(line_ls)

                    indx += 1

    print(masterlist)
    return masterlist


def internal_main_batch_multipass(main_outdir,repoofmods, massresolution, error_ppm):
    """
    Main method to run several analysis with several passes for one experimental ion file
    The ions that are matched will not be considered for matching in consecutive passes
    :return: void
    """

    #Parse Batch File
    batchfile = filedialog.askopenfilename(title='Load Batch File', filetypes=[('CSV', '.csv')])
    filestorun = intern_multipass_batchparser(batchfile)


    # expfiles = filedialog.askopenfilename(title='Load Experimental masses', filetypes=[('CSV', '.csv')])

    for expindx, explist in enumerate(filestorun):
        exp_file = explist[0]
        exp_name = exp_file.split("\\")
        print('Searching file {}, number {} of {}...'.format(exp_name[-1], expindx+1, len(filestorun)))
        #Parse experimental ions
        org_expions, proteinName = expion_parser(exp_file)

        #Use previously created fragments
        # load_ions_bool = messagebox.askyesno('Load Ions File?',
        #                                      'Do you want to load an existing .theofrags file (or calculate theoretical ions fresh from a template file)? Choose YES to load .ions file or NO to open a parameter template')

        load_ions_bool = explist[1]

        # A dict of analysis. Each one has a list of passes
        analysis_dict = {}

        # A dict of analysis. Each one hasdictionary of passes
        # analysis_dict_rev = {}
        paramobj_dict = ''
        param_name = ''
        if load_ions_bool:
            # theofrags_file = filedialog.askopenfilename(title='Choose Ions File', filetypes=[('Theoretical Fragments File', '.theofrags')])
            theofrags_file = load_ions_bool
            with open(theofrags_file, 'rb') as picklefile:
                analysis_dict = pickle.load(picklefile)

            # revtheofrags_file = filedialog.askopenfilename(title='Choose Reverse Ions File',
            #                                             filetypes=[('Theoretical Fragments File', '.theofrags')])
            # with open(revtheofrags_file, 'rb') as revpicklefile:
            #     analysis_dict_rev = pickle.load(revpicklefile)

        else:
            print("No ions file given, therefore moving on to obtaining params file!")
            param_file = explist[2]
            # Load theoretical database parameters
            if param_file:
                # paramfile = filedialog.askopenfilename(title='Load Parameter File', filetypes=[('CSV', '.csv')])
                paramobj_dict = Parameter_Parser.parse_param_template_batch_multipass(param_file)
                param_name = param_file.split("\\")
                print(f"Creating database using file = {param_name[-1]}")


            for analysis in paramobj_dict:
                print(f"Current Analysis = {analysis}")
                #Alwasy reset to this output path so that all passes are saved under one folder
                os.chdir(main_outdir)

                proteinSeq = ""
                #A dictionary of passes per analysis
                passes_dict = {}
                for paramobj in paramobj_dict[analysis]:
                    print(f"PassName = {paramobj.analysisName}")
                    # print(f"before ifragments = {paramobj.combodict_calc()}")

                    #pass sequence so output has original protein sequence
                    proteinSeq = paramobj.seq

                    # It actually returns a tuple with the analysis name in first position then the fragments dictionary
                    frag_dict = ifragments(analysis_name=paramobj.analysisName, sequence=proteinSeq, types=paramobj.iontypes,
                                           maxcharge=paramobj.maxcharge,
                                           maxstart=paramobj.min_len, maxlength=paramobj.max_len,
                                           modbool=paramobj.noncysmods,
                                           combo_dict=paramobj.combodict_calc(),
                                           cystine=paramobj.disulfide_bool,
                                           uniprot_offset=paramobj.uniprot_offset,
                                           allow_ssbroken=paramobj.ss_allowbroken,
                                           reverse_seq=False, foo_ls=paramobj.disulfide_ls,
                                           reduced_cys=paramobj.naturally_redcys, modsdictio=repoofmods, fragmentmode=paramobj.fragmentchem)
                    passes_dict[paramobj.analysisName] = frag_dict

                analysis_dict[analysis] = []
                analysis_dict[analysis].append(proteinSeq)
                analysis_dict[analysis].append(passes_dict)

                # Create pickle file to save theoretical database
                out_folder_str = f"{proteinName}_{analysis}"
                try:
                    os.mkdir(out_folder_str)
                except FileExistsError:
                    out_folder_str = out_folder_str + "_(1)"
                    os.mkdir(out_folder_str)


                print(f"param_name[-1] = {param_name[-1]}")
                fragments_to_picklefile(f"{param_name[-1]}", analysis_dict)
                fragments_to_FRAGSfile(f"{param_name[-1]}", str(analysis_dict))

                os.chdir(out_folder_str)

            # print(f'analysis_dict = {analysis_dict}')

            # print("Calculating Reverse Sequence Theoretical Fragments")
            #
            # for rev_analysis in paramobj_rev_dict:
            #     print(f"Current Analysis for Reverse Sequences = {rev_analysis} of {len(paramobj_dict)}")
            #     passes_rev_dict = {}
            #     for paramobj_rev in paramobj_rev_dict[rev_analysis]:
            #         print(f"PassName = {paramobj_rev.analysisName}")
            #
            #         #It actually returns a tuple with the analysis name in first position then the fragments dictionary
            #         frag_dict_rev = ifragments(analysis_name= paramobj_rev.analysisName, sequence=paramobj_rev.seq, types=paramobj_rev.iontypes, mincharge=paramobj_rev.mincharge,
            #                                maxcharge=paramobj_rev.maxcharge, maxstart=paramobj_rev.min_len,
            #                                maxlength=paramobj_rev.max_len,
            #                                modbool=paramobj_rev.noncysmods, max_mods=0, combo_dict=paramobj_rev.combodict_calc(),
            #                                cystine=paramobj_rev.disulfide_bool, uniprot_offset=paramobj_rev.uniprot_offset,
            #                                allow_ssbroken=paramobj_rev.ss_allowbroken, reverse_seq=True,
            #                                foo_ls=paramobj_rev.disulfide_ls, reduced_cys=paramobj_rev.naturally_redcys)
            #
            #         passes_rev_dict[paramobj_rev.analysisName] = frag_dict_rev
            #
            #     #Analsis dictionary with passes dictionaries as values
            #     analysis_dict_rev[rev_analysis] = passes_rev_dict
            #     out_folder_str_rev = f"{proteinName}_{rev_analysis}"
            #     os.chdir(main_outdir + f"/{out_folder_str_rev}")
            #     # Create pickle file
            #     fragments_to_picklefile(f"{out_folder_str_rev}_Reverse", analysis_dict_rev)
            #
            # # print(f'analysis_dict_rev = {analysis_dict_rev}')
            # print(f'Expions before analysis_num loop = {len(org_expions)}')


        #Matching process begins
        for analysis_num in analysis_dict:
            print(f'Matching currently Analysis # {analysis_num}')
            #Make a copy of the experimental ions list
            analysis_expions = org_expions.copy()
            print(f'Initial Exp ions in {analysis_num} = {len(analysis_expions)}')
            for pass_label in analysis_dict[analysis_num][1]:

                print(f"go back to folder = {proteinName}_{analysis_num}")
                # print(f"analysis_dict[analysis_num] = {analysis_dict[analysis_num]}")
                print(f"pass_label = {pass_label}")



                try:
                    os.chdir(f"{proteinName}_{analysis_num}")
                except FileNotFoundError:
                    os.mkdir(f"{proteinName}_{analysis_num}")
                    os.chdir(f"{proteinName}_{analysis_num}")


                print(f"analysis_dict = {analysis_dict.keys()}")
                matched_expions = matchmaker2_multipass(analysis_dict[analysis_num][1][pass_label], analysis_expions, massresolution, error_ppm,analysis_dict[analysis_num][0])
                #Removing matched ions from futures passes
                for ion in matched_expions:
                    if ion in analysis_expions:
                        analysis_expions.remove(ion)
                print(f'Expions after {pass_label}= {len(analysis_expions)}\n')

            unmatched_expions_outfile_writer(analysis_expions, analysis_num, csv=True)



if __name__ == '__main__':

    # print(f"x-c: {mass.Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1' + 'NH3')}")
    # print(f"x-cdot: {mass.Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1' + 'NH3' + 'H-1')}")
    # print(f"a-zdot: {mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1' + 'ON-1H-1')}")
    # isotope_xtractor()
    internal_main_batch_multipass()
    # main_batch()
    # main()












