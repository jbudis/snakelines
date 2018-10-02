from __future__ import print_function
from __future__ import division

from pandas import Series
import matplotlib
import collections as col
import numpy as np
import statsmodels.api as sm
from scipy import interpolate
import dill

matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get_bins_from_file(gc, nn):
    """
    Load the GC content and N content of the reference.
    :param gc: str - filename of the GC content pickle
    :param nn: str - filename of the N content pickle
    :return: 2x dict{window:dict{chromsome:ndarray}} - directory of window and chromsomes to GC/N contents of those windows
    """
    with open(gc, "r") as f:
        gc_bins_all = dill.load(f)
    with open(nn, "r") as f:
        nn_bins_all = dill.load(f)

    return gc_bins_all, nn_bins_all


def get_frag_counts_allch(pos, chroms, window, gc_bins):
    """
    Count fragment per bin.
    :param pos: ndarray - positions of the fragments
    :param chroms: ndarray - chromsomes of the fragments
    :param window: int - size of the bin
    :param gc_bins: dict{chromsome:ndarray} - directory of chromsomes to GC contents
    :return: ndarray(#chormosomes, #bins) - fragment count per chromosome and bin
    """
    max_length = max(len(gc_bins[k]) for k in gc_bins)

    # shape: (chromosome, length)
    frag_counts = np.zeros((len(gc_bins), max_length))

    binning = pos // window

    for ch in range(24):
        bins = np.bincount(binning[chroms == ch])[:max_length]
        frag_counts[ch, :bins.size] = bins

    return frag_counts


def get_mean(frag_counts, mean_precision, gc_bins, nn_bins):
    """
    From frag_counts make mean counts. Use mean_precision as precision (if window == mean_precision - take all possible bins)
    :param frag_counts: ndarray(#chormosomes, #bins) - fragment count per chromosome and bin
    :param mean_precision: int - precision for mean counts, maximum is window
    :param gc_bins: dict{chromsome:ndarray} - directory of chromsomes to GC contents
    :param nn_bins: dict{chromsome:ndarray} - directory of chromsomes to N contents
    :return: 2x pd.Series - series with mean count per gc content and number of fragments
    """
    gc_content = col.defaultdict(list)
    gc_res = col.defaultdict(float)

    # print type(indices), type(frag_counts), type(gc_bins)

    chromosomes = range(len(gc_bins))

    for ch in chromosomes:
        indices = nn_bins[ch] == 0
        # print len(frag_counts[ch]), len(gc_bins[ch]), len(nn_bins[ch])
        for f, g in zip(frag_counts[ch][:len(gc_bins[ch])][indices], gc_bins[ch][indices]):
            gc_content[int(g * mean_precision) / float(mean_precision)].append(f)

    for g in gc_content:
        gc_res[g] = np.mean(gc_content[g])

    return Series(gc_res.values(), index=gc_res.keys()), Series([len(gc) for gc in gc_content.values()], index=gc_content.keys())


def get_loess_mean(gc_res, gc_weight=None):
    """
    Calculate loess curve from MEAN fragment counts.
    :param gc_res: pd.Series - mean counts per gc_content
    :param gc_weight: ndarray/None - array of weights to use for weighted lowess, if None do not use weighted lowess
    :return: lowess, float - lowess smoothing object (lists of gc contents and their mean count), average
    """

    if gc_weight is not None:
        x = np.zeros(np.sum(gc_weight))
        y = np.zeros_like(x)

        start = 0
        for r in gc_res.index:
            x[start:start + gc_weight[r]] = r
            y[start:start + gc_weight[r]] = gc_res[r]
            start += gc_weight[r]
    else:
        x = gc_res.index
        y = gc_res

    loess_iter = 4
    loess_frac = 0.33
    loess = sm.nonparametric.lowess(y, x, frac=loess_frac, it=loess_iter)

    return loess, sum(y)/float(len(y))


def calculate_new_loess_allch(chroms, pos, window, gc_bins, nn_bins, mean_precision=None, weight_it=True, figname=None):
    """
    Calculate the lowess smoothing of the whole sample.
    :param chroms: ndarray - chromosome numbers of each fragment
    :param pos: ndarray - position of each fragment
    :param window: int - size of the bin
    :param gc_bins: dict{chromsome:ndarray} - directory of chromsomes to GC contents
    :param nn_bins: dict{chromsome:ndarray} - directory of chromsomes to N contents
    :param mean_precision: int - precision for mean counts, maximum is window
    :param weight_it: bool - whether to weight lowess
    :param figname: str - path to figure, where to save loess figure
    :return: ndarray - lowess weights of each fragment
    """
    if mean_precision is None:
        mean_precision = window

    # get framgent counts
    frag_counts = get_frag_counts_allch(pos, chroms, window, gc_bins)

    # get means
    gc_res, gc_weight = get_mean(frag_counts, mean_precision, gc_bins, nn_bins)

    # compute loess
    loess, average_loess = get_loess_mean(gc_res, gc_weight if weight_it else None)
    loess_x = map(lambda x: x[0], loess)
    loess_y = map(lambda x: x[1], loess)
    loess_f = interpolate.interp1d(loess_x, loess_y)

    # compute the average:
    used_bins = 0
    used_frac = 0

    for ch in range(min(chroms), max(chroms) + 1):
        ind_bins = nn_bins[ch] == 0  # & (frag_counts[ch][:len(nn_bins[ch])] > 0)
        used_bins += np.sum(ind_bins)

        pos_i = pos[(chroms == ch) & (pos < len(nn_bins[ch])*window)]
        binning = pos_i // window
        ind_frac = nn_bins[ch][binning] == 0
        used_frac += np.sum(ind_frac)

    average = used_frac / float(used_bins)

    # plot figure
    if figname is not None:
        plt.figure()
        plt.plot(gc_res.index, gc_res, "rx")
        loess_x = map(lambda x: x[0], loess)
        loess_y = map(lambda x: x[1], loess)
        plt.plot(loess_x, loess_y, "b-", lw=2)
        plt.xlim(0.3, 0.7)
        plt.xlabel("GC Content")
        plt.ylabel("Mean bin count")
        # plt.ylim(60, 140)
        # plt.title(str(sample)+" Window: "+str(window/1000)+"K")
        plt.plot([0, 1], [average, average], "g-")
        # plt.plot([0, 1], [average_loess, average_loess], "k-")
        plt.savefig(figname)

    weights = np.zeros_like(pos, dtype="float")

    # ok we have the base, now we can go chromosome per chromosome:
    for ch in range(min(chroms), max(chroms) + 1):
        indices = (chroms == ch) & (pos < len(nn_bins[ch])*window)
        pos_i = pos[indices]
        binning = pos_i // window

        gc_cont = np.floor(gc_bins[ch][binning] * mean_precision) / float(mean_precision)
        nn_cont = nn_bins[ch][binning]

        ind_frac = nn_cont == 0

        loess_inter = np.zeros_like(gc_cont)
        loess_inter[ind_frac] = loess_f(gc_cont[ind_frac])

        w_i = weights[indices]
        w_i[ind_frac] = (np.array(gc_res[gc_cont[ind_frac]]) - (loess_inter[ind_frac] - average)) / np.array(gc_res[gc_cont[ind_frac]], dtype="float")
        w_i[~ind_frac] = 0.0
        weights[indices] = w_i

    return weights


def get_loess_weights(chroms, pos, figname, gc_file, nn_file):
    """
    Calculate the lowess smoothing of the whole sample. Use default parameters.
    :param chroms: ndarray - chromosome numbers of each fragment
    :param pos: ndarray - position of each fragment
    :param figname: str - path to figure, where to save loess figure
    :param gc_file: str - path to GC file
    :param nn_file: str - path to NN file
    :return: ndarray - lowess weights of each fragment
    """
    if len(chroms) == 0 or len(pos) == 0:
        return np.array([], dtype='float')
    
    gc_bins, nn_bins = get_bins_from_file(gc_file, nn_file)

    window = 20000

    # select only positions on first 24 chromosomes (excluding the M chromosome)
    ind = chroms < 24
    loess = np.ones_like(pos, dtype="float")

    loess[ind] = calculate_new_loess_allch(chroms[ind], pos[ind], window, gc_bins[window], nn_bins[window], None, True, figname)

    return loess
