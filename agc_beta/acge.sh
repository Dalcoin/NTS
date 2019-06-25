make
python opti.py
f90 $F90FLAGS -o run -s -w eos_mat.f $LINK_FNL
./run
rm gam.srt
rm run
