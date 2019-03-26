#!/bin/csh

set VERSION = 3.0.7

set NPROC_X = 8
set NPROC_Y = 8
set NPROC_Z = 5
set NPROC_TOTAL = `echo $NPROC_X $NPROC_Y $NPROC_Z | gawk '{printf "%.0f\n",$1*$2*$3;}'`
echo $NPROC_TOTAL

set NX = 200
set NY = 200
set NZ = 75
set HH = 1.0
set NT = 5000
set DT = 0.02

set ORDER = 2
set ORDER = 4

set RUN_NAME = adj01
set SIM_DIR_ROOT = ../../../AdjSims/V3.0.7-a2a_xyz
#` pwd `
set MAIN_OUTPDIR = ${SIM_DIR_ROOT}/OutBin

set ADJ_FILE = ${SIM_DIR_ROOT}/Adj-SourceBin/my_adjoint_file.e3d
set ADJ_COMPS = "vx,vy,vz"

set TS_XYZ = 0
set TS_XYZ = 1
set DT_TS = 1
set DT_TS = 10
set DX_TS = 1
set DY_TS = 1
set DZ_TS = 1

set NSEIS = 0
set STATCORDS = ../../../StatInfo/fd_rt01-h0.400.statcords

set MODEL_STYLE = 1
#set VMODDIR = ../../Model/Models
#set VMODDIR = ~/scratch/large_files/BB-Valid/SoCal/LongBeach33-2018/Model/Models
#set PMOD = v_sgt-lb01-h0.050.p
#set SMOD = v_sgt-lb01-h0.050.s
#set DMOD = v_sgt-lb01-h0.050.d
set VMODDIR = ../../../Model/Models/
set PMOD = vp3dfile_in.p
set SMOD = vs3dfile_in.s
set DMOD = rho3dfile_in.d

#set MODEL_STYLE = 0
#set VMODDIR = ../../Model/Mod-1D
set FD_VMOD = iters
#nr02-vs850.fd-h0.400

set DEFAULT_PARFILE = ${SIM_DIR_ROOT}/e3d_default_new_h.par
set PARFILE = ${SIM_DIR_ROOT}/e3d.par

set MAXMEM = 12000

set MODEL_LAT = -42.0332006429
set MODEL_LON = 172.559006104
set MODEL_ROT = 0.0

set BFILT = 4
set FLO = 1.0

set DUMP_ITINC = 5000
set DUMP_FLAG = 0
set DUMP_FLAG = 1

set ENABLE_RESTART = 0
set READ_RESTART = 0
set MAIN_RESTARTDIR = ${SIM_DIR_ROOT}/Restart
set RESTART_ITINC = 20000

set SEISDIR = SeismoBin
set TMP_SEISDIR = ${SIM_DIR_ROOT}/${SEISDIR}
set TMP_SEISDIR = ${SIM_DIR_ROOT}/${SEISDIR}
set MAIN_OUTPDIR = ${SIM_DIR_ROOT}/OutBin
#set TMP_SEISDIR = ${HOME}/scratch/${SIM_DIR_ROOT}/${SEISDIR}

\cp $DEFAULT_PARFILE $PARFILE

echo "#"	>> $PARFILE
echo "# Run specific parameters start here"	>> $PARFILE
echo "#"	>> $PARFILE

echo "version=${VERSION}-mpi"		>> $PARFILE
echo "name=${RUN_NAME}"			>> $PARFILE

echo "report=100"			>> $PARFILE
echo "reportId=0"			>> $PARFILE

echo "maxmem=${MAXMEM}"                 >> $PARFILE
echo "order=$ORDER"          >> $PARFILE

echo "nx=$NX"                           >> $PARFILE
echo "ny=$NY"                           >> $PARFILE
echo "nz=$NZ"                           >> $PARFILE
echo "h=$HH"                            >> $PARFILE
echo "nt=$NT"                           >> $PARFILE
echo "dt=$DT"                           >> $PARFILE

echo "nproc_x=$NPROC_X"                           >> $PARFILE
echo "nproc_y=$NPROC_Y"                           >> $PARFILE
echo "nproc_z=$NPROC_Z"                           >> $PARFILE

echo "bfilt=$BFILT"				>> $PARFILE
echo "flo=$FLO"				>> $PARFILE
echo "fhi=0.0"				>> $PARFILE
echo "bforce=0"				>> $PARFILE
echo "dblcpl=0"				>> $PARFILE
echo "pointmt=0"                        >> $PARFILE
echo "ffault=0"                         >> $PARFILE

echo "adjoint=1"                         >> $PARFILE
echo "adjoint_file=${ADJ_FILE}"            >> $PARFILE
echo "adjoint_tfastest=1"            >> $PARFILE
echo "adjoint_comps=${ADJ_COMPS}"            >> $PARFILE

echo "model_style=${MODEL_STYLE}"                    >> $PARFILE
echo "vmoddir=${VMODDIR}"               >> $PARFILE

echo "pmodfile=$PMOD"                   >> $PARFILE
echo "smodfile=$SMOD"                   >> $PARFILE
echo "dmodfile=$DMOD"                   >> $PARFILE

echo "model=$FD_VMOD"                   >> $PARFILE

echo "qpfrac=100"                        >> $PARFILE
echo "qsfrac=50"                       >> $PARFILE
echo "qpqs_factor=2.0"                  >> $PARFILE

echo "fmax=25.0"                        >> $PARFILE
echo "fmin=0.01"                        >> $PARFILE
###echo "fmax=20.0"                        >> $PARFILE
###echo "fmin=0.05"                        >> $PARFILE

echo "modellon=$MODEL_LON"              >> $PARFILE
echo "modellat=$MODEL_LAT"              >> $PARFILE
echo "modelrot=$MODEL_ROT"              >> $PARFILE

echo "enable_output_dump=${DUMP_FLAG}"		>> $PARFILE
echo "dump_itinc=$DUMP_ITINC"		>> $PARFILE
echo "main_dump_dir=${MAIN_OUTPDIR}"    >> $PARFILE

echo "seisdir=${TMP_SEISDIR}"           >> $PARFILE
echo "nseis=$NSEIS"			>> $PARFILE
echo "seiscords=${STATCORDS}"           >> $PARFILE

echo "ts_xyz=$TS_XYZ"			>> $PARFILE
echo "ts_xy=0"			>> $PARFILE
echo "ts_xz=0"			>> $PARFILE
echo "ts_yz=0"			>> $PARFILE
echo "dtts=$DT_TS"                      >> $PARFILE
echo "dxts=$DX_TS"                      >> $PARFILE
echo "dyts=$DY_TS"                      >> $PARFILE
echo "dzts=$DZ_TS"                      >> $PARFILE

echo "enable_restart=${ENABLE_RESTART}" >> $PARFILE
echo "restartdir=${MAIN_RESTARTDIR}"    >> $PARFILE
echo "restart_itinc=${RESTART_ITINC}"   >> $PARFILE
echo "read_restart=${READ_RESTART}"     >> $PARFILE
echo "restartname=${RUN_NAME}"          >> $PARFILE

###echo "elas_only=1"          >> $PARFILE

echo "vpvs_max_global=3.0"          >> $PARFILE

####----------------------------------------------------
#### merge_P3-mpi parameters follow
####

set N_STRAIN_COMPS = 6
set STRAIN_COMPS = "exxp,eyyp,eii3,exyp,exzp,eyzp"

set MERGE_LOGDIR = Mlog
set FILELIST = ${SIM_DIR_ROOT}/merge_filelist_new.txt

set OUTFILE_ROOT = ${RUN_NAME}_xyzts

set MERGE_PARFILE = ${SIM_DIR_ROOT}/merge_P3_new.par
\rm $MERGE_PARFILE

echo "filelist=$FILELIST"		>> $MERGE_PARFILE
echo "n_comps=$N_STRAIN_COMPS"		>> $MERGE_PARFILE
echo "comp_names=$STRAIN_COMPS"		>> $MERGE_PARFILE
echo "outdir=$MAIN_OUTPDIR"		>> $MERGE_PARFILE
echo "outfile_root=$OUTFILE_ROOT"	>> $MERGE_PARFILE
echo "logdir=$MERGE_LOGDIR"		>> $MERGE_PARFILE
echo "logname=$RUN_NAME"		>> $MERGE_PARFILE

\rm $FILELIST
set ip = 0
while ( $ip < $NPROC_TOTAL )

echo $MAIN_OUTPDIR/${RUN_NAME}_xyzts- $ip | gawk '{printf "%s%.5d.e3d\n",$1,$2;}' >> $FILELIST

@ ip ++
end
