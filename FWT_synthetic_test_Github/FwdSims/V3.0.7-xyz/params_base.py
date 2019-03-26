run_name='2012p075555'
version='3.0.4'
bin_process_ver='slurm'
global_root='/nesi/project/nesi00213'
tools_dir='/nesi/project/nesi00213/EMOD3D/tools'
user_root='/nesi/project/nesi00213/RunFolder/tdn27'
run_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/'
sim_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/2012p075555'
lf_sim_root_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/2012p075555/LF'
hf_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/2012p075555/HF'
bb_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/2012p075555/BB'
figures_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/2012p075555/Figures'
srf_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/Sources/'
srf_files=['/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/Sources/2012p075555/Srf/2012p075555.srf']
hf_slips=['/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/Sources/2012p075555/Stoch/2012p075555.stoch']
# A single value of hf_kappa and hf_sdrop specified in params.py will be used by default
# If you wish to use a specific hf_kappa and hf_sdrop value for each SRF, uncomment and edit below
#hf_kappa_list=[]
#hf_sdrop_list=[]
vel_mod_dir='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555'
v_mod_1d_dir='/nesi/project/nesi00213/VelocityModel/Mod-1D'
params_vel='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555/params_vel.py'

#FROM VELECITY MODEL
MODEL_LAT = '-42.0332006429'
MODEL_LON = '172.559006104'
MODEL_ROT = '0'

#spatial grid spacing
hh = '0.4' #must be in formate like 0.200 (3decimal places)

#x,y,z grid size (multiple grid spacing
nx = '267'
ny = '269'
nz = '75'
sufx = '_rt01-h0.400'
sim_duration = '62.42'
dt = 0.0200
hf_dt = 0.0050
flo = '0.25'

#dir for vel_mod 
vel_mod_params_dir = '/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555'
GRIDFILE = '/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555/gridfile_rt01-h0.400' #gridout-x used to be referred to as GRIDFILE by gen_ts 
GRIDOUT = '/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555/gridout_rt01-h0.400'
#input for statgrid gen
MODEL_COORDS = '/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555/model_coords_rt01-h0.400'
MODELPARAMS = '/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555/model_params_rt01-h0.400'
MODEL_BOUNDS = '/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Data/VMs/2012p075555/model_bounds_rt01-h0.400'
stat_vs_est = '/nesi/project/nesi00213/StationInfo/geoNet_stats+2018-05-22.vs30'
stat_vs_ref= '/nesi/project/nesi00213/StationInfo/geoNet_stats+2018-05-22.vs30ref'
stat_file='/nesi/project/nesi00213/StationInfo/geoNet_stats+2018-05-22.ll'
STAT_FILES=[stat_file]
stat_coords='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/2012p075555/fd_rt01-h0.400.statcords'
FD_STATLIST='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/Runs/2012p075555/fd_rt01-h0.400.ll'
mgmt_db_location='/nesi/nobackup/nesi00213/RunFolder/tdn27/testrun/'
hf_stat_vs_ref='/nesi/project/nesi00213/StationInfo/geoNet_stats+2018-05-22.hfvs30ref'
