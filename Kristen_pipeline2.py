import kristen_traversal as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr
from csv import writer as csv_writer
import debug_renders as dbg
import numpy as np
from time import time
# TODO: Delete this pipeline when it is confirmed that it is no longer needed
# Goal of this pipeline
#     1. Detect the number of cells that were properly stained
#     2. For the successfully stained cells, determine how much GFP is located inside the mitochondria


translator = {'C1':0,
              'C3':1,
              'C4':2}

source = uf.Kristen_traverse('/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/20170209/Transfection C', matching_map=translator)

named_source = uf.name_channels(source, ['DAPI','GFP', 'mCherry'])

max_mCherry = cf.max_projection(named_source, in_channel = 'mCherry', out_channel = 'max_mCherry')

max_GFP = cf.max_projection(max_mCherry, in_channel = 'GFP', out_channel = 'max_GFP')

max_DAPI = cf.sum_projection(max_GFP, in_channel = 'DAPI', out_channel = 'sum_DAPI')

# smoothed_DAPI = cf.smooth_2d(max_DAPI, in_channel = 'sum_DAPI', smoothing_px=.5)

binarized_nuclei = cf.robust_binarize(max_DAPI, in_channel = 'sum_DAPI', out_channel = 'nuclei', _dilation=0, heterogeity_size=50, feature_size=50)

segmented_nuclei = cf.label_and_correct(binarized_nuclei, in_channel = ['nuclei', 'sum_DAPI'], out_channel = 'nuclei',min_px_radius=15, min_intensity=20)

stabilized_mCherry = cf.gamma_stabilize(segmented_nuclei, in_channel = 'max_mCherry', min='min', alpha_clean=.5)

smoothed_mCherry = cf.smooth_2d(stabilized_mCherry, in_channel = 'max_mCherry', smoothing_px=.5)

mCherry_o_n_segmented = cf.robust_binarize(smoothed_mCherry,
                                       in_channel='max_mCherry',
                                       out_channel='max_mCherry_binary',
                                       heterogeity_size=10, feature_size=250)

running_render = rdr.Kristen_render(mCherry_o_n_segmented, in_channel=['name pattern',
                                                                'group id',
                                                               'max_mCherry',
                                                               'max_mCherry_binary',
                                                               'GFP',
                                                               'mCherry'],
                                   out_channel='_',
                                    output='Kristen_Transfection_B_and_C_GFP_analysis_results_2.csv',
                                   save=False)

#
# Kristen_summary = rdr.Kristen_summarize_a(running_render, in_channel = ['name pattern', 'q_mean','q_median', 'q_std', 'nq_mean', 'nq_median', 'nq_std', 'slope', 'r2', 'p'],
#                                out_channel='_',
#                                output='Kristen_Transfection_B_and_C_GFP_analysis_results.csv')

for i in enumerate(running_render):
    print 'Analyzed group %s - image %s' % (['group id'], ['name pattern'])