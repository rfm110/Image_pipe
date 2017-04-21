# the following is the orioinal kristen render to which modifications are being made
@generator_wrapper(in_dims=(None,None, 2, 2, 3, 3), out_dims=(None,))
def Kristen_render(name_pattern,
                   group_id,
                   mCherry,
                   extranuclear_mCherry_pad,
                   GFP_orig,
                   mCherry_orig, output,
                   save=False, directory_to_save_to='verification'):
    labels, _ = ndi.label(extranuclear_mCherry_pad)
    unique_segmented_cells_labels = np.unique(labels)[1:]
    mCherry_cutoff = np.zeros_like(mCherry)
    qualifying_cell_label = []
    qualifying_regression_stats = []

    for cell_label in unique_segmented_cells_labels:
        mCherry_2 = np.zeros_like(mCherry)
        my_mask = labels == cell_label
        average_apply_mask = np.mean(mCherry[my_mask])
        intensity = np.sum(mCherry[my_mask])
        binary_pad = np.zeros_like(mCherry)
        binary_pad[my_mask] = 1
        pixel = np.sum(binary_pad[my_mask])

        if (average_apply_mask > .05 or intensity > 300) and pixel > 4000:

            GFP_limited_to_cell_mask = cf._3d_stack_2d_filter(GFP_orig, my_mask)
            mCherry_limited_to_cell_mask = cf._3d_stack_2d_filter(mCherry_orig, my_mask)

            qualifying_3d_GFP = GFP_limited_to_cell_mask[mCherry_limited_to_cell_mask>50]
            average_3d_GFP = np.mean(qualifying_3d_GFP)
            median_3d_GFP = np.median(qualifying_3d_GFP)
            std_3d_GFP = np.std(qualifying_3d_GFP)
            sum_qualifying_GFP = np.sum(qualifying_3d_GFP)

            nonqualifying_3d_GFP = GFP_limited_to_cell_mask[mCherry_limited_to_cell_mask<=50]
            average_nonqualifying_3d_GFP = np.mean(nonqualifying_3d_GFP)
            median_nonqualifying_3d_GFP = np.median(nonqualifying_3d_GFP)
            std_nonqualifying_3d_GFP = np.std(nonqualifying_3d_GFP)
            sum_nonqualifying_GFP = np.sum(nonqualifying_3d_GFP)

            sum_total_GFP = sum_qualifying_GFP + sum_nonqualifying_GFP
            percent_qualifying_over_total_GFP = sum_qualifying_GFP/sum_total_GFP
            # report the percentage too or sums are sufficient?

            GFP_orig_qualifying = cf._3d_stack_2d_filter(GFP_orig, my_mask)
            mCherry_orig_qualifying = cf._3d_stack_2d_filter(mCherry_orig, my_mask)
            mCherry_1d = mCherry_orig_qualifying[mCherry_orig_qualifying > 50]
            GFP_1d = GFP_orig_qualifying[mCherry_orig_qualifying>50]
            regression_results = stats.linregress(GFP_1d, mCherry_1d)

            mCherry_2[my_mask] = mCherry[my_mask]
            mCherry_cutoff[my_mask] = mCherry[my_mask]
            qualifying_cell_label.append(cell_label)
            qualifying_regression_stats.append((regression_results[0], regression_results[2], regression_results[3]))

            name_pattern_split = name_pattern.split(' - ')
            transfection_label = name_pattern_split[0]
            cell_type = name_pattern_split[1]
            exp_time = name_pattern_split[2]
            image_number = name_pattern_split[4]

            with open(output, 'ab') as output_file:
                writer = csv_writer(output_file, delimiter='\t')
                writer.writerow([transfection_label, cell_type, exp_time, image_number, cell_label, sum_qualifying_GFP, sum_total_GFP, average_3d_GFP, median_3d_GFP, std_3d_GFP, average_nonqualifying_3d_GFP, median_nonqualifying_3d_GFP, std_nonqualifying_3d_GFP, regression_results[0], regression_results[2], regression_results[3]])

            plt.figure(figsize=(26.0, 15.0))
            plt.title('Kristen\'s Data')
            plt.suptitle(name_pattern)

            main_ax = plt.subplot(221)
            plt.subplot(221, sharex=main_ax, sharey=main_ax)
            plt.title('mCherry Binary')
            im = plt.imshow(extranuclear_mCherry_pad, interpolation='nearest', cmap = 'hot')
            plt.colorbar(im)
            plt.subplot(222, sharex=main_ax, sharey=main_ax)
            plt.title('mCherry')
            plt.imshow(mCherry, interpolation='nearest')
            plt.contour(extranuclear_mCherry_pad, [0.5], colors='k')
            plt.subplot(223)
            dplt.better2D_desisty_plot(GFP_1d, mCherry_1d)
            plt.title('mCherry Intensity as a Function of GFP Voxel')
            plt.xlabel('GFP Voxel')
            plt.ylabel('mCherry Intensity')
            plt.subplot(224, sharex=main_ax, sharey=main_ax)
            plt.title('mCherry-cutoff applied')
            plt.imshow(mCherry_2, interpolation='nearest')

            if not save:
                plt.show()

            else:
                name_puck = directory_to_save_to + '/' + 'Kristen-' + name_pattern+ '_cell' + str(cell_label)+ '.png'
                plt.savefig(name_puck)
                plt.close()
    plt.figure(figsize=(26.0, 15.0))
    main_ax = plt.subplot(121)
    plt.subplot(121, sharex=main_ax, sharey=main_ax)
    plt.suptitle('mCherry Before and After Qualifying Cell Cutoff is Applied')
    plt.title('mCherry')
    im = plt.imshow(mCherry, interpolation='nearest')
    plt.colorbar(im)
    plt.subplot(122, sharex=main_ax, sharey=main_ax)
    plt.title('mCherry')
    plt.imshow(mCherry_cutoff, interpolation='nearest')
    if not save:
        plt.show()

    else:
        name_puck = directory_to_save_to + '/' + 'Kristen-' + name_pattern + 'cutoff_app' + '.png'
        plt.savefig(name_puck)
        plt.close()

    return qualifying_regression_stats
