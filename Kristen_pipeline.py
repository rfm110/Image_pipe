import traversals as uf
import core_functions as cf
import render as rdr
from csv import writer as csv_writer
from time import time

# Goal of this pipeline
#     1. Detect the number of cells that were properly stained/transfected
#             quantification only for successfully transfected cells
#     2. For the successfully stained cells, determine how much GFP is located inside the mCHerry-stained mitochondria


translator = {'C1':0,
              'C3':1,
              'C4':2}

source = uf.Kristen_read_file('/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/20170209', matching_map=translator)

# start_csv = uf.csv_writer()
# source = uf.Kristen_yield()

named_source = uf.name_channels(source, ['DAPI','GFP', 'mCherry'])

max_mCherry = cf.max_projection(named_source, in_channel='mCherry', out_channel='max_mCherry')

max_GFP = cf.max_projection(max_mCherry, in_channel='GFP', out_channel='max_GFP')

stabilized_mCherry = cf.gamma_stabilize(max_GFP, in_channel='max_mCherry', min='min', alpha_clean=.5)

smoothed_mCherry = cf.smooth_2d(stabilized_mCherry, in_channel='max_mCherry', smoothing_px=.5)

mCherry_o_n_segmented = cf.robust_binarize(smoothed_mCherry,
                                       in_channel='max_mCherry',
                                       out_channel='max_mCherry_binary',
                                       heterogeity_size=10, feature_size=250)

# running_render = rdr.Kristen_render(mCherry_o_n_segmented, in_channel=['name pattern',
#                                                                 'group id',
#                                                                'max_mCherry',
#                                                                'max_mCherry_binary',
#                                                                'GFP',
#                                                                'mCherry'],
#                                     out_channel='_',
#                                     output='Kristen_Transfection_B_and_C_GFP_analysis_results.csv',
#                                     save=True)
running_render = rdr.Kristen_render(mCherry_o_n_segmented, in_channel=['name pattern',
                                                                'group id',
                                                               'max_mCherry',
                                                               'max_mCherry_binary',
                                                               'GFP',
                                                               'mCherry'],
                                    out_channel=['my_mask', 'cell_label'])

GFP_limited_to_cell_mask = cf. for_each(running_render, cf._3d_stack_2d_filter, 'per_cell', in_channel = ['GFP', 'my_mask'], out_channel = 'GFP_limited_to_cell_mask')
mCherry_limited_to_cell_mask = cf. for_each(running_render, cf._3d_stack_2d_filter, 'per_cell', in_channel = ['mCherry', 'my_mask'], out_channel = 'mCherry_limited_to_cell_mask')
analysis = rdr.Kristen_quantification_analysis(mCherry_limited_to_cell_mask, in_channel = ['name_pattern', 'GFP_limited_to_cell_mask', 'mCherry_limited_to_cell_mask', 'max_mCherry'])

# HOW to replace cell label in for loop with 'per_cell' channel and 'for_each' method
#
# Kristen_summary = rdr.Kristen_summarize_a(running_render, in_channel = ['name pattern', 'q_mean','q_median', 'q_std', 'nq_mean', 'nq_median', 'nq_std', 'slope', 'r2', 'p'],
#                                out_channel='_',
#                                output='Kristen_Transfection_B_and_C_GFP_analysis_results.csv')

for i in enumerate(running_render):
    print 'Analyzed group %s - image %s' % (['group id'], ['name pattern'])