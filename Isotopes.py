"""
Isotopes: Palying with pyteomics and isitope envelopes for matching interfrags
CRR
02/27/20
"""
from pyteomics import mass
# from IsoSpecPy import IsoSpecPy
from math import exp
from brainpy import isotopic_variants
from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
import scipy
import pythoms.molecule as pythmole


water = mass.calculate_mass(formula='H[2]2O') # Heavy water

heavy_water = mass.calculate_mass(formula='H[2]HO') # Semiheavy water

print(water)
print(heavy_water)

p = mass.Composition(sequence = 'PEPTIDEPEPTIDEPEPTIDE')
print(p)
theo_iso_env = mass.most_probable_isotopic_composition(sequence="PEPTIDEPEPTIDEPEPTIDE")
print(theo_iso_env)

peptide_dict = {'C': 102, 'H': 159, 'N': 21, 'O': 43, 'S': 0, 'Fe': 0}

print(peptide_dict.keys())

outstr = ''
for key in peptide_dict.keys():
    print(key)
    if peptide_dict[key] == 0:
        pass
    else:
        outstr += f"{key}{peptide_dict[key]}"
print(f"outstr = {outstr}")

#'C': 154, 'H': 246, 'N': 49, 'O': 40, 'S': 1, 'Fe': 1
mol = pythmole.IPMolecule(outstr, charge = 4, resolution = 17000)
# print(pythmole.IPMolecule('L2PdAr+I', charge = 1))
# print(pythmole.IPMolecule('C', charge = 1, resolution=17000))
details = mol.print_details()
print(details)
masss = mass.calculate_mass(composition=peptide_dict, charge=4)
print(f"mass = {masss}")
# print(details[0])
pythmole_mz = mol.bar_isotope_pattern[0]
pythmole_int = mol.bar_isotope_pattern[1]
print(f"pythmole_mz = {pythmole_mz}")
mz_grid = np.arange(pythmole_mz[0] - 1,
                    pythmole_mz[-1] + 1, 0.02)
intensity = np.zeros_like(mz_grid)
sigma = 0.0003
index = 0
for peak in pythmole_int:
    # Add gaussian peak shape centered around each theoretical peak
    intensity += peak * np.exp(-(mz_grid - pythmole_mz[index]) ** 2 / (2 * sigma)) / (np.sqrt(2 * np.pi) * sigma)
    index += 1

# Normalize profile to 0-100
intensity = (intensity / intensity.max())

plt.plot(mz_grid, intensity, label = "pythmole")
plt.show()


#Finally, you can find the most probable isotopic composition for a substance with pyteomics.mass.most_
# probable_isotopic_composition() function. The substance is specified as a formula, a pyteomics.mass.
# Composition object or a modX sequence string.

# pep = mass.most_probable_isotopic_composition(sequence='PEPTIDE')
# print(pep)
#
# comp_one = mass.isotopic_composition_abundance(formula='H2O') # Water with an unspecified isotopic state
# print(comp_one)
# comp_two = mass.isotopic_composition_abundance(formula='H[2]2O')
# print(comp_two)
# comp_three = mass.isotopic_composition_abundance(formula='H[2]H[1]O') # Semiheavy water
# print(comp_three)
# print("Peptide")
# pep_ab_one = mass.isotopic_composition_abundance(formula='C34H53O15N7')
# print(pep_ab_one)
# pep_ab_two = mass.isotopic_composition_abundance(formula='C[12]33C[13]1H53O15N7')
# print(pep_ab_two)
# pep_ab_three = mass.isotopic_composition_abundance(formula='C[12]33C[13]1H[1]52H[2]1O15N7')
# print(pep_ab_three)
# pep_ab_fr = mass.isotopic_composition_abundance(formula='C[12]594C[13]4H594O1141N96S5')
# print(pep_ab_fr)

hsa = mass.most_probable_isotopic_composition(sequence='CAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKH')
print(hsa[0])

# hsa_celeven = mass.most_probable_isotopic_composition(sequence = 'DAHKSEVAHRF')
# print(hsa_celeven)

hsa_celevenw = mass.isotopic_composition_abundance(formula='C[12]55C[13]1H85O171N19')
print(hsa_celevenw)

# #IsoSpec
# i = IsoSpecPy.IsoSpec.IsoFromFormula("H2O1", 0.9)
#
# print(i)
# print(len(i))
#
# confs = i.getConfs()
# print(f"Mass: {confs}")
# print(f"isotopelogues log probabilities: {confs[1]}")
#
#
# hydrogen_probs = (0.9, 0.1)
# oxygen_probs = (0.5, 0.3, 0.2)
# hydrogen_masses = (1.00782503207, 2.0141017778)
# oxygen_masses = (15.99491461956, 16.99913170, 17.9991610)
# atom_counts = (2, 1)
#
# m = IsoSpecPy.IsoSpec(atom_counts, (hydrogen_masses, oxygen_masses), (hydrogen_probs, oxygen_probs), 0.9)
# confsm = m.getConfs()
#
#
# def plot_Isospec(config_obj):
#     mz_ls = config_obj[0]
#     probls = []
#     for prob in config_obj[1]:
#         # print(prob)
#         # print(100*(exp(prob)))
#         probls.append(100 * (exp(prob)))
#
#     maximun = max(probls)
#
#     probls_norm = []
#
#     for prob in probls:
#         probnomr = prob / maximun
#         probls_norm.append(probnomr)
#
#     plt.scatter(mz_ls, probls_norm)
#
#
# plot_Isospec(confs)
# plot_Isospec(confsm)



print("\nBRAIN\n")

#BRAIN
# Generate theoretical isotopic pattern
# peptide = {'H': 571, 'C': 376, 'O': 103, 'N': 102, 'S':6, "Fe":1}
peptide = {'C': 376, 'H': 569, 'N': 101, 'O': 103, 'S': 4, 'Fe': 1}
# {'C': 408, 'H': 634, 'N': 111, 'O': 112, 'S': 3}
peptide_o = {'C': 154, 'H': 246, 'N': 49, 'O': 40, 'S': 1, 'Fe': 1}
#It accepts a charge of 0
masss = mass.calculate_mass(composition=peptide, charge=0)
print(f"mass = {masss}")
massz = mass.calculate_mass(composition=peptide, charge = 3)
print(f"mass = {massz}")

print("Peptide 1")
theoretical_isotopic_cluster = isotopic_variants(peptide, npeaks=15, charge=3)
mz_ls = []
int_ls = []
for peak in theoretical_isotopic_cluster:
    # print(peak.mz, peak.intensity)
    mz_ls.append(peak.mz)
    int_ls.append(peak.intensity)

# produce a theoretical profile using a gaussian peak shape
print(f"Foo = {theoretical_isotopic_cluster[0]}")
print(f"Foo2 = {theoretical_isotopic_cluster[0].mz}")
print(f"Foo3 = {theoretical_isotopic_cluster[-1].mz}")
mz_grid = np.arange(theoretical_isotopic_cluster[0].mz - 1,
                    theoretical_isotopic_cluster[-1].mz + 1, 0.02)
intensity = np.zeros_like(mz_grid)
sigma = 0.003
for peak in theoretical_isotopic_cluster:
    # Add gaussian peak shape centered around each theoretical peak
    intensity += peak.intensity * np.exp(-(mz_grid - peak.mz) ** 2 / (2 * sigma)
            ) / (np.sqrt(2 * np.pi) * sigma)

# Normalize profile to 0-100
intensity = (intensity / intensity.max()) * 100

#It accepts a charge of 0
theoretical_isotopic_cluster_o = isotopic_variants(peptide_o, npeaks=15, charge=4)

mzm_ls = []
intm_ls = []
for peak in theoretical_isotopic_cluster_o:
    # print(peak.mz, peak.intensity)
    mzm_ls.append(peak.mz)
    intm_ls.append(peak.intensity)

# produce a theoretical profile using a gaussian peak shape
import numpy as np
mz_grid_o = np.arange(theoretical_isotopic_cluster_o[0].mz - 1,
                    theoretical_isotopic_cluster_o[-1].mz + 1, 0.02)
intensity_o = np.zeros_like(mz_grid_o)
sigma = 0.003
for peak in theoretical_isotopic_cluster_o:
    # Add gaussian peak shape centered around each theoretical peak
    intensity_o += peak.intensity * np.exp(-(mz_grid_o - peak.mz) ** 2 / (2 * sigma)
            ) / (np.sqrt(2 * np.pi) * sigma)

# Normalize profile to 0-100
intensity_o = (intensity_o / intensity_o.max()) * 100

# draw the profile
# print(np.append(mz_grid, mz_grid_o))
# print(np.append(intensity, intensity_o))
plt.plot(mz_grid, intensity, label = "Fragment + 3 a.m.u")
# plt.plot(mz_grid_o, intensity_o, label = "Fragment")
# plt.plot(np.append(mz_grid, mz_grid_o),np.append(intensity, intensity_o), label = "Adding intensities")
plt.legend(loc='best')
# plt.scatter(mz_ls, int_ls)

plt.xlabel("m/z")
plt.ylabel("Relative intensity")
plt.plot()
plt.show()

# print(int_ls)
# print(intm_ls)
# print(np.corrcoef(int_ls,intm_ls))
# diff = np.asarray(mz_ls) - np.asarray(mzm_ls)
# print(diff)
# print(sum(diff))

# print("Find peaks")
#
# peaks = scipy.signal.find_peaks(intensity)
# print(f"Peaks = {peaks}")
# prominances = scipy.signal.peak_prominences(intensity, peaks[0])
# print(f"Prominances = {prominances}")
# print(prominances[0])
#
# print(f"Intensity = {len(intensity)} and Intensity_o = {len(intensity_o)}")
# corr = np.corrcoef(intensity,intensity_o)
# print(f"Corr = {corr}")
# for x in corr:
#     print(x)
#     for ele in x:
#         print(ele)
# print(f"Corr score = {corr[0][1]}")
#

