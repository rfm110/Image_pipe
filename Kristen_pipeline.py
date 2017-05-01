import traversals as uf
import core_functions as cf
import render as rdr
from csv import writer as csv_writer
from csv import DictWriter
from time import time

# Goal of this pipeline
#     1. Detect the number of cells that were properly stained/transfected
#             quantification only for successfully transfected cells
#     2. For the successfully stained cells, determine how much GFP is located inside the mCHerry-stained mitochondria


translator = {'C1':0,
              'C3':1,
              'C4':2}



source = uf.Kristen_read_file('/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/20170209', matching_map=translator)
start_csv = uf.csv_writer(source)
kristen_yield = uf.Kristen_yield(start_csv)

named_source = uf.name_channels(kristen_yield, ['DAPI','GFP', 'mCherry'])

max_mCherry = cf.max_projection(named_source, in_channel='mCherry', out_channel='max_mCherry')

max_GFP = cf.max_projection(max_mCherry, in_channel='GFP', out_channel='max_GFP')

stabilized_mCherry = cf.gamma_stabilize(max_GFP, in_channel='max_mCherry', min='min', alpha_clean=.5)

smoothed_mCherry = cf.smooth_2d(stabilized_mCherry, in_channel='max_mCherry', smoothing_px=.5)

mCherry_o_n_segmented = cf.robust_binarize(smoothed_mCherry,
                                       in_channel='max_mCherry',
                                       out_channel='max_mCherry_binary',
                                       heterogeity_size=10, feature_size=250)


# CHANGES FOR RESTRUCTURING made below
running_render = rdr.Kristen_GFP_cutoff(mCherry_o_n_segmented, in_channel=['name pattern',
                                                                           'group id',
                                                                           'max_mCherry',
                                                                           'max_mCherry_binary'],
                                    out_channel=['my_mask', 'labels'])

per_cell_split = cf.splitter(running_render, 'per_cell',
                             sources=['mCherry', 'GFP',
                                      'max_mCherry', 'max_GFP'],
                             mask='my_mask')


analysis = cf.for_each(per_cell_split, rdr.Kristen_quantification_and_stats, 'per_cell',
                                                                in_channel = ['GFP',
                                                                              'mCherry',
                                                                              'max_mCherry'],

                                                                out_channel = ['mCherry_cutoff',
                                                                               'GFP_1d',
                                                                               'mCherry_1d',
                                                                               'sum_qualifying_GFP',
                                                                               'sum_total_GFP',
                                                                               'average_3d_GFP',
                                                                               'median_3d_GFP',
                                                                               'std_3d_GFP',
                                                                               'average_nonqualifying_3d_GFP',
                                                                               'median_nonqualifying_3d_GFP',
                                                                               'std_nonqualifying_3d_GFP',
                                                                               'slope_linregress',
                                                                               'r_value',
                                                                               'p_value'])
# Need to render the following (tile them first)
#     mCherry   mCherry_2   mCherry_cutoff  GFP_1d  mCherry_1d

mCherry_tiled = cf.tile_from_mask(analysis, 'per_cell', 'max_mCherry')
mCherry_2_tiled = cf.tile_from_mask(mCherry_tiled, 'per_cell', 'mCherry_cutoff')
mCherry_cutoff_tiled = cf.tile_from_mask(mCherry_2_tiled, 'per_cell', 'mCherry_cutoff')
max_mCherry_orig = cf.max_projection(mCherry_cutoff_tiled, in_channel='mCherry', out_channel='max_mCherry_orig')
render = rdr.Kristen_image_render(max_mCherry_orig, in_channel = ['name pattern',
                                                          'max_mCherry_binary',
                                                          'max_mCherry_orig',
                                                          'mCherry_cutoff'],
                                            out_channel='_',
                                            save=False,
                                            directory_to_save_to='verification')

# write_csv = rdr.Kristen_write_to_csv(render, in_channel =   ['name pattern',
#                                                              'cell_label',
#                                                              'sum_qualifying_GFP',
#                                                              'sum_total_GFP',
#                                                              'average_3d_GFP',
#                                                              'median_3d_GFP',
#                                                              'std_3d_GFP',
#                                                              'average_nonqualifying_3d_GFP',
#                                                              'median_nonqualifying_3d_GFP',
#                                                              'std_nonqualifying_3d_GFP',
#                                                              'regression_results'],
#                                      out_channel = '_',
#                                      output='Kristen_Transfection_B_and_C_GFP_analysis_results.csv')
write_csv = rdr.Kristen_write_to_csv(render,
                                     output='Kristen_Transfection_B_and_C_GFP_analysis_results.csv')
# HOW to replace cell label in for loop with 'per_cell' channel and 'for_each' method
#
# Kristen_summary = rdr.Kristen_summarize_a(running_render, in_channel = ['name pattern', 'q_mean','q_median', 'q_std', 'nq_mean', 'nq_median', 'nq_std', 'slope', 'r2', 'p'],
#                                out_channel='_',
#                                output='Kristen_Transfection_B_and_C_GFP_analysis_results.csv')

# TODO: add header
with open('Kristen_Transfection_B_and_C_GFP_analysis_results.csv', 'ab') as output_file:
    writer = csv_writer(output_file, delimiter = '\t')
    writer.writerow(['Transfection Label', 'Image Number', 'Cell Number',
                     'Qualifying GFP', 'Total GFP',
                     'Average-Qualifying GFP', 'Median-Qualifying GFP','Standard Deviation-Qualifying GFP',
                     'Average-Nonqualifying GFP','Median-Nonqualifying GFP', 'Standard Deviation-Nonqualifying GFP',
                     'Slope', 'R-Value', 'P-Value'])


for i in enumerate(write_csv):
    print 'Analyzed group %s - image %s' % (['group id'], ['name pattern'])