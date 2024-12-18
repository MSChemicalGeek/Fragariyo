"""Testing my own peak picking script
https://medium.com/@chrisjpulliam/quickly-finding-peaks-in-mass-spectrometry-data-using-scipy-fcf3999c5057
"""
import os.path
from tkinter import filedialog
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks
import numpy as np
import ms_deisotope
import time


# https://stackoverflow.com/questions/63177236/how-to-calculate-signal-to-noise-ratio-using-python
def signaltonoise(ms_spectrum, axis=0, ddof=0, sectionum = 10, abovethreshold=3):
    a = np.asanyarray(ms_spectrum["int"])
    b = np.asanyarray(ms_spectrum["mz"])
    # Slicing intensity array in three sections. Rounding up because I need integers for slicing.
    listofslices_int= np.array_split(a, sectionum)
    listofslices_mz = np.array_split(b, sectionum)

    results = []
    for indx in range(sectionum):
        # Calculating means
        # print(f"indx = {indx}")
        meancalculation = listofslices_int[indx].mean(axis)
        print(meancalculation)
        results.append([listofslices_int[indx], listofslices_mz[indx], meancalculation*abovethreshold])


    return results

def peaklist_fromrawcsv(numberofsections, timesabovenoise):
    coordinatesfile = filedialog.askopenfilename(title='Select Coordinate files', filetypes=[('CSV', '.csv')])

    ms_data = pandas.read_csv(coordinatesfile, header=1)

    ms_data.columns = ["mz", "int"]
    # ms_data.set_index("mz")
    # print(ms_data)
    # default kind is line
    ms_data.plot(x='mz', y='int', color='blue')
    plt.title('MS spectrum')
    plt.xlabel('m/z')
    plt.ylabel('Intensity')

    mzslices_with_thresholds = signaltonoise(ms_data, axis=0, ddof=0, sectionum = numberofsections, abovethreshold=timesabovenoise)
    # LIst of slices (int,mz) and mean calculations

    selectedmz_arr = np.array([])
    for dataregion in mzslices_with_thresholds:
        # print(f"dataregions = {dataregion}")
        region_mean = dataregion[2]
        region_mz = dataregion[1]
        region_int = dataregion[0]
        # print(len(region_mz), len(region_int))
        peak_idx, _ = find_peaks(region_int,
                                 prominence=100,
                                 height=region_mean,
                                 distance=None)

        # convert region indeces into mzvalues

        result_mz = region_mz[peak_idx]

        selectedmz_arr = np.concatenate((selectedmz_arr, result_mz), axis=0)

    allindeces_arr = []
    for mzvalue in selectedmz_arr:
        index_val = ms_data.index[ms_data['mz'] == mzvalue]
        # From 1D array to scalar
        index_scalar = index_val.item()
        allindeces_arr.append(index_scalar)
    # print(allindeces_arr)
    peak_data = ms_data.iloc[allindeces_arr]
    # print(peak_data)
    sns.scatterplot(x=peak_data["mz"], y=peak_data["int"],
                    color='red', alpha=0.5)

    plt.savefig("peak_picking.png")

    #
    # print(peak_data)
    # print(peak_data.index)
    # print(peak_data.values.tolist())

    # # It seems I don't need to peak picked...it will do it for you...never mind is taking froever! (more than 12 hrs)
    # # need to prepare first into PeakSet Object: https://github.com/mobiusklein/ms_deisotope?tab=readme-ov-file
    start_time = time.time()
    preparedpeaks = ms_deisotope.deconvolution.utils.prepare_peaklist(peak_data.values.tolist())
    print(preparedpeaks)

    print("--- %s seconds preparing ---" % (time.time() - start_time))

    start_time = time.time()
    deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(preparedpeaks, averagine=ms_deisotope.peptide,
                                                           scorer=ms_deisotope.MSDeconVFitter(10.),
                                                           use_quick_charge=False, left_search_limit=2,
                                                           right_search_limit=2)
    print(deconvoluted_peaks)

    print("--- %s seconds deconvoluting ---" % (time.time() - start_time))

    outputstr = "mz,z,intensity,score,envelope\n"
    deisotope_mz = []
    deisotope_int =[]
    for peak in deconvoluted_peaks:
        if peak.score >= 190 and len(peak.envelope) > 2:
            outputstr += f"{peak.mz},{peak.charge},{peak.intensity},{peak.score},{peak.envelope}\n"
            deisotope_mz.append(peak.mz)
            deisotope_int.append(peak.intensity)

    # where to save
    outputfolder = os.path.dirname(coordinatesfile)
    # File path for the CSV file
    csv_file_path = os.path.join(outputfolder, 'deconvolutedpeaklist.csv')

    with open(csv_file_path, 'w') as f:
        # Write some text to the file
        f.write(outputstr)

    sns.scatterplot(x=deisotope_mz, y=deisotope_int,
                    color='green', alpha=0.90)

    plt.savefig("peak_picking.png")


if __name__ == '__main__':

   peaklist_fromrawcsv(10, 3)