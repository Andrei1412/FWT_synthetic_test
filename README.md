# FWT_synthetic_test
Full waveform tomography for synthetic test case

1. Create test configuration: Run Test_configuration.py in the Master folder
- Create initial and true models,move the velocity files to Model/Mdels
- create SOURCE.txt and STATION.txt for source and station information, move the files to StaInform

2. Create observed data: Run observed_data_nshots.py in FwdSims/V3.0.7-dc_stf_file 
- Record seismograms from forward emod3d simulation for the list of station at every source, using the same source type, source signature defined in e3d_obs.par.
- Move the seismogram folders to Kernels/Vel_ob

3. Prepare for inversion:
- Set the filter for functions: Read_GP_adj_calc_i_new_displacement_ts.py, Read_GP_obs_calc_err_rwm_new.py,  Read_GP_opt_calc_err_rwm_new.py in Kernels (To calculate the adjoint source, the misfit error at iteration and the misfit error for searching optimal step length).
- Prepare gradient taper matrixs: run Tape_Geometry.py,  Tape_Gradients_nShot_exp.py in Kernels/Iters/iter1 to generate a tape matrix for each source (damping near source/ receiver effects) and a tape matrix for the whole domain (smoothening kernels close to the x,y-boundaries) to Kernels/Iters/iter1/Dump. Recommend to call in the corresponding slurm script since the job may take up to 2 hours for 8 sources and domain of 200x200x75.

4. Invert models: Run Master_iterations_python.py or Master_iterations_matlab.py (Master python script) in Kernels/Iters/iter1 to do inversion (with adjoint source calculated in python and Matlab accordingly). The steps of inversion include:

For each source:

i) Generate the synthetic data for the initial (current) models: Call FWT_emod3d_shot_i_part1.sl in Master python script
  - Output seismograms for the list of stations and the deviatoric strain tensors at every cell in the domain, using the same source type, source signature defined in e3d_mysource_xyz_default.par and merge_P3_xyz.par in folder Kernels/Iters/iter1/FWD/.
  - Move the seismogram folders to Kernels/Vel_es and the strain wavefields to Kernels/Dev_Strain

ii) Calculate the adjoint source:
  - Call adj_test.sl in Master python script or export the seismogram to matlab and import back the adjoint source. The adjoint sources in ascii format for every station are moved to AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii
  
iii) Adjoint (backward) simulation: Call FWT_emod3d_shot_i_part2.sl in Master python script
  - Convert the adjoint sources in ascii format to binary format for adjoint simulation for the initial (current) models using parameters set in set_run+merge_params_new_h.csh in Kernels/Iters/iter1/ADJ
  - Output the deviatoric strain tensors at every cell in the domain and move the files to Kernels/Dev_Strain 

iv) Calculate the kernels for each idividual source:
  - Call kernel_shot_i_iter1.sl to calculate the kernels for Vs and Vp and move the kernels for the corresponding source to Kernels/Iters/iter1/All_shots
  
For all sources:

v) Sum and precondition the kernels for all sources using taper matrices to get the gradients for Vs and Vp. Smoothen (using Gausian filter) and constrain the gradients.
  - Save the gradients for the current iterations to Kernels/Iters/iter1/Dump.
  - Load the gradients from the previous iteration if  conjugate gradient method is applied.
  
vi) Calculate the optimal step length for model updating by call in Step_length_3D_pad2 function in the Master python script:
  - if the optimal step length exists, update and save the model to Kernels/Iters/iter1/Dump for the current iteration.
  - if the optimal step length does not exist, terminate the iteration process or start the current iteration again with higher frequency data.
  





