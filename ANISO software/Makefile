F90=ifort
F77=ifort
CC=ifort
AAAFFLAGS= -O3 -mcmodel=large   
F90FLAGS= -O3 -mcmodel=large 
F77FLAGS= -traceback  -O3
SPECIALFLAG= -d2
MODULEFLAG=
CCFLAGS=
LD=ifort -mkl=sequential 
LIBS= 
srclayer=aniso_lmsvars.o aniso_globalpar.o anisolayers_main.o

.SUFFIXES :
.SUFFIXES : .o .c .f .f90 
.f90.o:
	$(F90) $(F90FLAGS) $(MODULEFLAG) -c $<
.f.o:
	$(F77) $(F77FLAGS) -c $<
.c.o:
	$(CC) $(CCFLAGS) -c $<
aniso:	${srclayer}
	${LD} ${F90FLAGS} ${modules} ${lapack} ${srclayer} ${LIBS} -o ./aniso.x
#dsyevd_files.o : 
#	$(F90) $(SPECIALFLAG) -c dsyevd_files.f
#zheevd_files.o :
#	$(F90) $(SPECIALFLAG) -c zheevd_files.f
#lattice_com.o: kinds.o
#parameters_mod.o : kinds.o
src2: ${modules}
clean:
	rm *.o *.mod
