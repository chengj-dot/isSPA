# isSPA
To compile the C codes, you should have installed EMAN 1.9. Also check your g++ version, 4.8 and 5.x works fine to me,but 6.3 occurs errors.	
For example, use "g++-4.8  apply_my_weighting_optimized_v2.c  -o apply_my_weighting_
optimized_v2  -IwhereyouinstallEMAN/include  whereyouinstallEMAN/lib/libEM.so" to complie.

Workflow:
1.	Typical SSNR parameters by fitting power spectrum with two exponential functions are given  --  ssnr.txt.  The first three parameters are the estimation of noise spectrum and the last three parameters are for signal spectrum);
2.	Before apply weighting, you should have achieved an image file in .lst format containing CTF parameters. Raw images listed in this file are then filtered by our optimized weighting function with known SSNR parameters using program “apply_my_weighting_optimized_v2”. Specific usage can be found by typing enter after the program. The size of raw images should be smaller than 1500*1500.
This operation will generate the filtered images and correlated image txt file in .lst format. Run lstfast.py on the output .lst file.
3.	Detect target proteins in the filtered images by using program “do_global_detection_with_
my_weight”. The input file is the output file in .lst of step 2. The parameter --threshold is up to each dataset, you can run this program with threshold 0 on one image, and you can choose a proper value as the threshold according to the distribution of scores. The constant –kk must be the same as you used in step 2.
4.	The output of program presents a large group of locations and orientations, then we need merge the neighbor solutions according to their translations and rotations by doing: “remove_repeat_particles_from_list.py locations_and_orientations_file.lst  center_distance _ threshold_in_pixels  euler_distance_threshold_in angles  output.lst” on each file.
5.	Convert each locations_and_orientations_file from the filtered image to original image by executing: “read_my-output_write_lst.py merged_ locations_and_orientations_file original_ image_file scale_factor output.lst”, if the filtered images have not been re-scaled, the scale factor is 1, otherwise depending on the scaling. Run lstfast.py on the output file.
6.	Extract potential particles from raw images by doing: “e2proc2d.py output_of_step5  ???.spi  --process=xform.centerbyxform –clip=particle_size,particle_size”
7.	Convert output_of_step5 to particle file with “convert_to_particle_file.py ???.spi particle_file.lst boxsize 0”
8.	Convert particle_file in .lst format to relion (.star) format for 3D classification: “lstfile_2_relion_class3d.py particle_file.lst ???.spi boxsize pixel_value particle.star metadata_line 0”, metadata_line is the number of lines before particles in lstfile.
9.	Run 3D classification by skipping alignment.
10.	Select proper particles according to 3D classification to local auto-refinement.
11.	If there exists model bias effects, you need do the following operations.
Firstly, convert relion star file to EMAN lst file by using “relion2lst_dont_write_particles.py refined_data.star –lst refined_data.lst –boxsize boxsize”
Secondly, calculate phase residual between each particle and projection generated by 3D map from local auto-refinement. You can execute jalign or our home-made program “calculate_cc_with_fsc –F FSC_file –a pixel_size –v acceleration voltage –H high_resolution_threshold –L low_resolution_threshold –r 3D_map –i refined_data.lst –f first_image_number –l last_image_number –R mask_radius_in_pixel –u soft_edge –o out.lst –c 0 “.
Thirdly, run “sort_by_score.py out.lst number_of_particle_to_stay metadata_line scored.lst”
Finally, convert lst file to relion star file as step 8 suggests, and do local auto-refinement.
