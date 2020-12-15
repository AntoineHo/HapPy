#!/usr/bin/python
# -*- coding: utf-8 -*-

# General
import os, sys, math, json
from itertools import combinations, chain # Obtaining peaks combinations

# Happy
from happy.utils import *
from happy.plot import *

# Stats
from scipy.stats import pearsonr
from scipy.optimize import curve_fit # Fitting multiple gaussian model
from scipy.signal import savgol_filter  # Smoothing frequency histogram
from scipy.signal import find_peaks  # Finding peaks
from scipy.signal import peak_widths # Obtaining peaks widths

def auto_estimate_haploidy(
    infile, outfile:str, size: str, min_peak: int, peak_prominence: int,
    window_search: float, score: float, plot = False, debug = False,
):
    """Finds peaks and modality, then computes scores of haploidy"""

    # Get histogram file and genome size estimation
    HIST = check_files(infile)
    SIZE = size_from_string(size)
    outfile = os.path.abspath(os.path.join(os.getcwd(), outfile))
    outdir = os.path.dirname(outfile)
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)

    # DEBUG
    if debug :
        print("-- Debugging --")
        print("Plot is {} and debug is {}".format(plot, debug))
        print("Input:", infile, type(infile))
        print("min_peak:", min_peak, type(min_peak))
        print("prominence:", peak_prominence, type(peak_prominence))
        print("window search:", window_search, type(window_search))
        print("Score threshold:", score, type(score))
        print("Size:", SIZE, type(SIZE))
        print("Output:", outfile, type(outfile))
        print("-- End debugging --")

    print("# Hap.py autoestimate")
    print("Coverage histogram:\t{0}\nOutput file:\t{1}\n".format(HIST, outfile))
    print(
        "===============================================================================\n"
    )

    # Read coverage histogram
    log("Reading histogram!")
    freqs = [int(line.split()[1]) for line in open(HIST, "r")][1:-1]

    # Apply Savitzky-Golay filter to smooth coverage histogram
    log("Analysing curve...")
    smoothed = [s for s in savgol_filter(freqs, 41, 3) if s >= 0] # WARNING might be a bad idea to remove 0.0
    peaks, props = find_peaks(smoothed, height=min_peak, prominence=peak_prominence)
    heights = [smoothed[i] for i in peaks]
    widths = peak_widths(smoothed, peaks)[0]  # Get peak widths
    #peak_props_dc = {pk:{"w":wdh, "h":hgt} for pk, wdh, hgt in zip(peaks, widths, heights)}

    if debug :
        print("Peaks:", peaks)
        debug_smooth_histogram(freqs, smoothed, peaks, heights, widths, outdir)

    # Consider each peak in a row to be the haploid coverage peak
    json_models = {}
    for one_allele_peak in peaks :
        if debug :
            print("Considering haploid peak: {}".format(one_allele_peak))

        # consider each possible number of alleles found in the assembly
        for n_alleles in [1,2,3,4,5,6] : # up to 6 alleles remaining in the assembly
            # Get a curve function (multiple gaussians) according to the number of alleles (up to 6 alleles)
            curve_function = get_curve_function(n_alleles)
            # Find expected peaks and score their multipliers: one_allele_peak * i/n_alleles <-- modifier
            expected_peaks = []
            peak_multipliers = {}
            for i in range(1, n_alleles+1) :
                exp_peak = int((i*one_allele_peak)/n_alleles)
                expected_peaks.append(exp_peak)
                peak_multipliers[exp_peak] = i/n_alleles

            # Get array of integers between 0 and len smoothed
            x = np.arange(len(smoothed))
            # Create a list of expectations for curve fitting (mu, sigma, amplitude)
            p0 = [item for sublist in [(mu,15,np.mean(heights)) for mu in expected_peaks] for item in sublist]
            # Create bounds for the search space :
            # - mu +/- 0.15*mu (mu = average = peak position)
            # - sigma between 0 and 20 <-- before == fourty but it is too high
            # - amplitude between 0 and +inf
            lower_bounds = [item for sublist in [(mu-(int(0.15*mu)+1),0,0) for mu in expected_peaks] for item in sublist]
            upper_bounds = [item for sublist in [(mu+(int(0.15*mu)+1),20,np.inf) for mu in expected_peaks] for item in sublist]

            # Try curve fitting with SciPy
            try :
                popt, pcov = curve_fit(curve_function, x, smoothed, p0=p0, bounds=(lower_bounds, upper_bounds))
                peak_opts = {} # Store each peak detected curve parameters
                for i, peak in enumerate(expected_peaks) :
                    peak_opts[peak] = popt[i*3 : i*3 + 3]

                # Write information to model
                # Model name = one_allele_peak + "_" + n_alleles
                model_name = str(one_allele_peak) + "_" + str(n_alleles)
                json_models[model_name] = {
                                            "peaks_models":{k:list(v) for k,v in peak_opts.items()},
                                            "one_allele_peak":int(one_allele_peak),
                                            "max_ploidy":int(n_alleles),
                                          }

                if debug :
                    print("lower:", lower_bounds)
                    print("actual", peak_opts)
                    print("upper:", upper_bounds)
                    print("")

            # Reports errors but do not break!
            except RuntimeError as e:
                print("RuntimeError: {}".format(e))
                continue
            except ValueError as e :
                print("ValueError: {}".format(e))
                continue
            except  :
                print("UnknownError!")
                continue

            # Computes correlation between model and fitted function
            corr, p = pearsonr(smoothed, curve_function(x, *popt))
            #corr = round(corr, 4)

            # If bad correlation --> skip : bad model
            if corr < 0.5 :
                continue

            # NOTE: find a better way to normalize this!
            # Computes a distance score to the actual peak distribution
            cumulative_distance = 0
            for current_peak in expected_peaks :
                all_peaks = list(peaks)+[0,len(smoothed)]
                closest_peak = min(all_peaks, key=lambda x:abs(x-current_peak))
                diff = abs(current_peak-closest_peak) / len(smoothed) # NOTE: bad normalisation
                cumulative_distance += diff
            distance_score = 1.0 - cumulative_distance

            # Computes area (with multiplier) of each peak in the model
            total_peak_area = 0
            peak_areas = {}
            for peak, multiplier in peak_multipliers.items() :
                peak_curve = gauss(x, *peak_opts[peak]) # Subcurve in the total model
                current_peak_area = sum(peak_curve) # Area of the subcurve
                peak_areas[peak] = current_peak_area # Store area
                total_peak_area += current_peak_area * multiplier # Computes area with modifier

            # Computes various scores and metrics:
            # total area no modifier except the one allele peak
            multi_allele_area = sum(area for peak, area in peak_areas.items() if peak != one_allele_peak)
            # one allele peak area
            one_allele_area =  peak_areas[one_allele_peak]
            # AUC score = 1 - (multi allele peak areas / one_allele_area)
            auc_ratio = 1 - (multi_allele_area / one_allele_area)
            # haploidy score = one_allele_area / (sum area / divider)
            total_haploid_area = 0
            for peak, multiplier in peak_multipliers.items() :
                total_haploid_area += peak_areas[peak] * multiplier
            haploidy_score = one_allele_area / total_haploid_area
            # difference of expected total size WITH MODIFIERS
            # divergence of total sequences found some are added or removed
            divergence = SIZE - total_peak_area
            # Overcollapsed if divergence is bigger than 0 because computed size is BELOW total size
            is_overcollapsed = True if divergence >= 0 else False
            # Percent of over/under collapsing
            percent_size = (divergence/SIZE)*100
            # Size score is proportional to the absolute divergence divided by SIZE
            size_score = 1 - (abs(divergence) / SIZE)

            # merged_model_score = size_score * haploidy_score * corr
            merged_model_score = size_score * haploidy_score * corr
            json_models[model_name]["score"] = merged_model_score

            one_to_one, correspondance = one_to_one_correspondance(
                                         peaks,
                                         expected_peaks,
                                         )
            if one_to_one :
                # get correspondance between expected and actual in reverse
                reversed_correspondance = {}
                for k,v in correspondance.items() :
                    reversed_correspondance[v] = k

                actual_one_allele_peak = correspondance[one_allele_peak]

                # check_limits between peaks
                limits = find_peaks_limits(peaks, smoothed)
                # Computes area (with multiplier) of each peak in the model
                actual_total_peak_area = 0
                actual_peak_areas = {}
                for expected_peak, peak in correspondance.items() :
                    start, end = limits[peak]
                    current_peak_area = sum(smoothed[start:end]) # Area of the subcurve
                    actual_peak_areas[peak] = current_peak_area # Store area
                    # Need expected peak to guess the right multiplier to apply
                    actual_total_peak_area += current_peak_area * peak_multipliers[expected_peak]

                actual_multi_allele_area = sum(area for peak, area in actual_peak_areas.items() if peak != actual_one_allele_peak)
                actual_one_allele_area =  actual_peak_areas[actual_one_allele_peak]
                actual_auc_ratio = 1 - (actual_multi_allele_area / actual_one_allele_area)
                actual_total_haploid_area = 0
                for expected_peak, peak in correspondance.items() :
                    actual_total_haploid_area += actual_peak_areas[peak] * peak_multipliers[expected_peak]
                actual_haploidy_score = actual_one_allele_area / actual_total_haploid_area
                actual_divergence = SIZE - actual_total_peak_area
                actual_is_overcollapsed = True if actual_divergence >= 0 else False
                actual_percent_size = (actual_divergence/SIZE)*100
                actual_size_score = 1 - (abs(actual_divergence) / SIZE)
                actual_data = (correspondance, actual_peak_areas, actual_haploidy_score,
                               actual_size_score, actual_divergence)
            else :
                pass

            # Write stats for this model
            auto_write_stats(
            outfile,
            n_alleles,
            one_allele_peak,
            peak_areas,
            round(distance_score, 4),
            round(corr, 4),
            round(haploidy_score, 4),
            total_peak_area,
            round(size_score, 4),
            is_overcollapsed,
            divergence,
            percent_size,
            one_to_one,
            actual_data
            )

            if plot :
                outname = "One_allele_peak_{}_curves_{}".format(one_allele_peak, n_alleles)
                plot_model(x, smoothed, curve_function, popt, peak_opts,
                               expected_peaks, peaks, peak_areas,
                               peak_multipliers, outdir, outname,
                               )

    # Write models information to JSON
    write_json_peaks(outfile, json_models)

    log("Finished!")
    sys.exit(0)

def get_curve_function(n) :
    """Programmatically create a curve function for fitting n gaussian"""
    funcstring = "def curve_function(x, {}) : ".format(", ".join("mu{i}, sigma{i}, A{i}".format(**{"i":i}) for i in range(1,n+1)))
    funcstring += "\n\treturn np.sum([" + ",".join("A{i}*np.exp(-(x-mu{i})**2/2/sigma{i}**2)".format(**{"i":i}) for i in range(1,n+1)) + "], axis=0)"
    exec(funcstring)
    return locals()['curve_function']

def auto_write_stats(
    outfile,
    n_alleles,
    one_allele_peak,
    peak_areas,
    distance_score,
    corr,
    haploidy_score,
    total_peak_area,
    size_score,
    is_overcollapsed,
    divergence,
    percent_size,
    one_to_one,
    actual_data,
):

# #Curves, OneAllelePeak, PeaksUsed, PeakDistanceScore, ModelPearsonR, OneAlleleBases, OneAlleleScore, SizeScore, Collapsed, SizeDivergence
# n_alleles, one_allele_peak, expected_peaks, distance_score, corr, one_allele_area, haploidy_score, SizeScore, over/under, divergence

    if not os.path.isfile(outfile) :
        # Writing format
        f = open(outfile, "w")

        header_string = "##source=HapPy\n"
        header_string += "##Curves=<Type=Integer,Description=\"Number of gaussian curves used in model. Approx. number of uncollapsed allelic versions in the assembly\">\n"
        header_string += "##OneAllelePeak=<Type=Integer,Description=\"Peak position (coverage) considered as corresponding to one allele (fully collapsed) in the model\">\n"
        header_string += "##ModelPeakCoverage=<Type=[Integer list],Description=\"Peak coverage in the model\">\n"
        header_string += "##ModelPeakAreas=<Type=[Integer list],Description=\"Number of bases under the peak (highest coverage corresponds to the cumulative size of fully collapsed alleles in the assembly if the model is correct)\">\n"
        header_string += "##ModelPeakDistanceScore=<Type=Float,Description=\"Score of how well the expected peaks fit the observed peak pattern. Best: 1.0\">\n"
        header_string += "##PearsonR=<Type=Float,Description=\"Pearson's R correlation between model and observed curves\">\n"
        header_string += "##HaploidyScore=<Type=Float,Description=\"\">\n"
        header_string += "##TotalPeakArea=<Type=Integer,Description=\"Computed genome size according to the model (taking peak multipliers into account)\">\n"
        header_string += "##SizeScore=<Type=Float,Description=\"\">\n"
        header_string += "##Collapsed=<Type=String,Description=\"Over or Under: according to the model: how is the assembly collapsed (over: too few sequences found or under: too many sequences found)\">\n"
        header_string += "##Divergence=<Type=Integer,Description=\"Total size of divergence between actual assembly size and estimated genome size\"\n"
        header_string += "##%Diff=<Type=Float,Description=\"Percent of total size of divergence compared to the given estimated genome size\"\n"
        header_string += "##1To1=<Type=Bool,Description=\"True if there is a one to one correspondance between peaks and expected peaks\"\n"
        header_string += "##ActualPeaks=<Type=[Integer list],Description=\"Actual peaks if there is a one to one correspondance\"\n"
        header_string += "##ActualPeakArea=<Type=[Integer list],Description=\"Actual peaks area if there is a one to one correspondance\"\n"
        header_string += "##ActualHaploidyScore=<Type=Float,Description=\"\"\n"
        header_string += "##ActualSizeScore=<Type=Float,Description=\"\"\n"
        header_string += "##ActualDivergence=<Type=Integer,Description=\"\"\n"
        header_string += "#Curves\tOneAllelePeak\tModelPeakCoverage\tModelPeakArea\tModelPeakDistanceScore\tPearsonR\tHaploidyScore\tTotalPeakArea\tSizeScore\tCollapsed\tDivergence\t%Diff\t1To1\tActualPeaks\tActualPeakAreas\tActualHaploidyScore\tActualSizeScore\tActualDivergence\n"
        f.write(header_string)
        f.close()

    # Build output string
    outstring = ""
    outstring += "{}\t".format(n_alleles)
    outstring += "{}\t".format(one_allele_peak)
    outstring += "{}\t".format(",".join(str(k) for k, v in peak_areas.items()))
    outstring += "{}\t".format(",".join(str(int(v)) for k, v in peak_areas.items()))
    outstring += "{}\t".format(distance_score)
    outstring += "{}\t".format(corr)
    outstring += "{}\t".format(haploidy_score)
    outstring += "{}\t".format(int(total_peak_area))
    outstring += "{}\t".format(size_score)
    outstring += "{}\t".format("Over" if is_overcollapsed else "Under")
    outstring += "{}\t".format(int(divergence))
    outstring += "{}%\t".format(round(percent_size, 3))
    outstring += "{}\t".format(one_to_one)
    if one_to_one :
        outstring += "{}\t".format(",".join(str(v) for k, v in actual_data[0].items())) # get correspondance
        outstring += "{}\t".format(",".join(str(int(v)) for k, v in actual_data[1].items()))
        outstring += "{}\t".format(round(actual_data[2],4))
        outstring += "{}\t".format(round(actual_data[3],4))
        outstring += "{}\n".format(int(actual_data[4]))
    else :
        outstring += "/\t/\t/\t/\t/\n".format()
    # (correspondance, actual_peak_areas, actual_haploidy_score, actual_size_score, actual_divergence)
    #outstring += "{}\n".format()

    f = open(outfile, "a+")
    f.write(outstring)
    f.close()

def one_to_one_correspondance(observed_peaks, expected_peaks) :
    """Check if all expected peaks have one only observed peak corresponding"""
    correspondance = {}
    for pk in expected_peaks :
        closest_peak = min(observed_peaks, key=lambda x:abs(x-pk))
        divergence = abs(pk - closest_peak)
        correspondance[pk] = closest_peak

    # return (True, correspondance) if no duplicated closest peaks
    # else return (False; None)
    if len(set(correspondance.values())) == len(correspondance.values()) :
        return True, correspondance
    else :
        return False, None

def find_peaks_limits(peaks, curve) :
    # Find peak limits in a given curve
    # A peak limit is the lowest value found in the expected curve between two expected peaks (contam or not)
    peaks = sorted(peaks)
    peak_limits = {}
    previous_peak_limit = 0

    for p in range(len(peaks)-1) :
        subcurve = curve[peaks[p]: peaks[p+1]]
        minfreq = min(subcurve)
        for relpos, freq in enumerate(subcurve) :
            if freq == minfreq :
                break
        left_limit_position = peaks[p] + relpos
        peak_limits[peaks[p]] = (previous_peak_limit, left_limit_position)
        previous_peak_limit = left_limit_position

    peak_limits[peaks[-1]] = (previous_peak_limit, len(curve)) # to the end of the distribution

    return peak_limits


def write_json_peaks(outfile, json_models) :
    """Write model information to JSON file"""

    # Debug
    """
    for k,v in json_models.items() :
        for kk, vv in v.items() :
            if kk == "peaks_models" :
                print(type(v[kk]))
    """

    # Remove bad models
    out_models = {}
    for model, data in json_models.items() :
        store = True
        for key in ["score", "one_allele_peak", "max_ploidy"] :
            if key not in data.keys() :
                store = False
                break
        if store :
            out_models[model] = data

    out = os.path.join(os.path.dirname(outfile), os.path.splitext(os.path.basename(outfile))[0] + ".model.json")
    f = open(out, "w")
    json.dump(out_models, f)
    f.close()
