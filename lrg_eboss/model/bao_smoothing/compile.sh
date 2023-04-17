gcc main_olin_cf.c  -lgsl -lm -lgslcblas cubature.c  -I/OPT/gsl-1.16/include/ -L//OPT/gsl-1.16/lib/ -lm -o file.out

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/OPT/gsl-1.16/lib/

./file.out 15
