#!/bin/csh

set HH = 0.400

set INFILE = nr02-vs850.vmod
set OUTFILE = nr02-vs850.fd-h$HH

set VPMIN = -1
set VSMIN = 0.5
set DNMIN = -1
set QPMIN = -1
set QSMIN = -1

set QSFAC = 50
set QP2QS = 2

gawk -v h=$HH -v vpm=$VPMIN -v vsm=$VSMIN -v dnm=$DNMIN -v qpm=$QPMIN -v qsm=$QSMIN \
-v qsfac=$QSFAC -v qp2qs=$QP2QS \
'BEGIN{printf "DEF HST\n";dd=0.0;nlay=0;}{ \
if(substr($0,1,1)!="#"){ \
nlay++; \
dd=dd+$1; \
zz[nlay]=dd; \
vp[nlay]=$2;if(vp[nlay]<vpm)vp[nlay]=vpm; \
vs[nlay]=$3;if(vs[nlay]<vsm)vs[nlay]=vsm; \
dn[nlay]=$4;if(dn[nlay]<dnm)dn[nlay]=dnm; \
qs[nlay]=$6;if(qsfac>0.0)qs[nlay]=qsfac*vs[nlay];if(qs[nlay]<qsm)qs[nlay]=qsm; \
qp[nlay]=$5;if(qp2qs>0.0)qp[nlay]=qp2qs*qs[nlay];if(qp[nlay]<qpm)qp[nlay]=qpm; \
}} \
END{ \
h2=0.5*h;for(i=1;i<=nlay;i++){zz[i]=int((zz[i]+h2)/h + 0.5)*h - h2;} \
for(i=1;i<nlay;i++){if(zz[i]>0.0 && \
zz[i]<zz[i+1] && \
vp[i]!=vp[i+1] && \
vs[i]!=vs[i+1] && \
dn[i]!=dn[i+1] && \
qp[i]!=qp[i+1] && \
qs[i]!=qs[i+1]) \
printf "%8.4f %8.4f %8.4f %10.2f %10.2f %11.5f\n",vp[i],vs[i],dn[i],qp[i],qs[i],zz[i];} \
printf "%8.4f %8.4f %8.4f %10.2f %10.2f %11.5f\n",vp[nlay],vs[nlay],dn[nlay],qp[nlay],qs[nlay],zz[nlay];}' \
$INFILE > $OUTFILE
