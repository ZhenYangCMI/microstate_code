GroupICAT v3.0a Updates (October 1, 2013):

1. File icatb/icatb_helper_functions/icatb_run_dfnc.m is fixed to avoid averaging of session timecourses when computing dFNC correlations. 

GroupICAT v3.0a Updates (Sep 12, 2013):

1. Feature specific defaults are added in mancova batch file.

2. Added batch script to do dFNC.

3. Option is provided to use low frequency and high frequency values to compute fALFF.

4. Sliding window computation using tukey window is disabled.

5. Added an option to use component network names in plotting FNC correlation matrix in mancova. 

The following files are modified or added:

1. icatb/icatb_batch_files/input_mancovan.m
2. icatb/icatb_batch_files/input_dfnc.m
3. icatb/icatb_display_functions/icatb_component_viewer.m
4. icatb/icatb_helper_functions/icatb_loadComp.m
5. icatb/icatb_helper_functions/icatb_dfnc_batch.m
6. icatb/icatb_helper_functions/icatb_plot_FNC.m
7. icatb/icatb_helper_functions/icatb_dfnc_options.m
8. icatb/icatb_helper_functions/icatb_setup_dfnc.m
9. icatb/icatb_helper_functions/icatb_run_dfnc.m
10. icatb/icatb_mancovan_files/icatb_mancovan_batch.m
11. icatb/icatb_mancovan_files/icatb_setup_mancovan.m
12. icatb/icatb_mancovan_files/icatb_mancovan_feature_options.m
13. icatb/icatb_mancovan_files/icatb_run_mancovan.m
14. icatb/icatb_mancovan_files/icatb_get_spec_stats.m
15. icatb/icatb_mancovan_files/icatb_ggmix.m
16. icatb/icatb_mancovan_files/icatb_display_mancovan.m

GroupICAT v3.0a Updates (July 05, 2013):

1. IVA code is fixed to work on older versions of MATLAB. The following files are modified:
	a. icatb/icatb_analysis_functions/icatb_algorithms/icatb_iva_laplace.m
	b. icatb/icatb_analysis_functions/icatb_algorithms/icatb_iva_second_order.m

GroupICAT v3.0a Updates (June 18, 2013):

1. Copyright information is added in file icatb/icatb_analysis_functions/icatb_algorithms/icatb_gigicar.m.

GroupICAT v3.0a Updates (June 13, 2013):

1. Orthogonal viewer button outside of display GUI in SBM is modified to not use functional data files.

GroupICAT v3.0a Updates (June 11, 2013):

1. Spatial-temporal regression utility is now added in SBM -> Utilities.

GroupICAT v3.0a Updates (June 02, 2013):

1. File icatb/icatb_helper_functions/icatb_rename_4d_file.m is modified to generate spaces at the end of the nifti file numbers. 
