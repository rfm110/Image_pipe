from wrappers import generator_wrapper
import core_functions as cf


gamma_stabilize = generator_wrapper(in_dims=(None,))(cf.gamma_stabilize)

smooth = generator_wrapper(cf.smooth)

smooth_2d = generator_wrapper(in_dims=(2,))(cf.smooth_2d)

sum_projection = generator_wrapper(in_dims=(3,), out_dims=(2,))(cf.sum_projection)

max_projection = generator_wrapper(in_dims=(3,), out_dims=(2,))(cf.max_projection)

random_walker_binarize = generator_wrapper(in_dims=(2,))(cf.random_walker_binarize)

robust_binarize = generator_wrapper(in_dims=(2,), out_dims=(2,))(cf.robust_binarize)

voronoi_segment_labels = generator_wrapper(in_dims=(2,))(cf.voronoi_segment_labels)

filter_labels = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(cf.filter_labels)

exclude_region = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(cf.exclude_region)

in_contact = generator_wrapper(in_dims=(2, 2))(cf.in_contact)

improved_watershed = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(cf.improved_watershed)

label_and_correct = generator_wrapper(in_dims=(2, 2,), out_dims=(2,))(cf.label_and_correct)

label_based_aq = generator_wrapper(in_dims=(2, 2), out_dims=(1, 2))(cf.label_based_aq)

aq_gfp_per_region = generator_wrapper(in_dims=(2, 2, 2), out_dims=(1, 2))(cf.aq_gfp_per_region)

detect_upper_outliers = generator_wrapper(in_dims=(1,), out_dims=(1, 1, None))(cf.detect_upper_outliers)

paint_mask= generator_wrapper(in_dims=(2, 1), out_dims=(2,))(cf.paint_mask)

mask_filter_2d = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(cf.mask_filter_2d)

clear_based_on_2d_mask = generator_wrapper(in_dims=(3, 2), out_dims=(3,))(cf.clear_based_on_2d_mask)

binarize_3d = generator_wrapper(cf.binarize_3d)

volume_mqvi = generator_wrapper(in_dims=(3, 3), out_dims=(None,))(cf.volume_mqvi)

volume_aqvi = generator_wrapper(in_dims=(3, 3), out_dims=(None,))(cf.volume_aqvi)

_3d_mask_from_2d_mask = generator_wrapper(in_dims=(3, 2), out_dims=(3,))(cf._3d_mask_from_2d_mask)

binarize_2d = generator_wrapper(in_dims=(2,), out_dims=(2,))(cf.binarize_2d)

agreeing_skeletons = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(cf.agreeing_skeletons)

classify_fragmentation_for_mitochondria = generator_wrapper(in_dims=(2, 2), out_dims=(None, 2, 2, 2))(cf.classify_fragmentation_for_mitochondria)

locally_normalize = generator_wrapper(in_dims=(3,), out_dims=(3,))(cf.locally_normalize)
