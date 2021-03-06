F90=ifort
F77=ifort
CC=ifort
AAAFFLAGS= -O3   
F90FLAGS= -O  
F77FLAGS= -O3  
SPECIALFLAG= -d2
MODULEFLAG=
CCFLAGS=  
LD=ifort -mkl=sequential 
src0= modules.o CBESSELnew.o fftpack5.1d.o breathingsphere.o
#
.SUFFIXES :
.SUFFIXES : .o .c .f .f90 
.f90.o:
	$(F90) $(F90FLAGS) $(MODULEFLAG) -c $<
.f.o:
	$(F77) $(F77FLAGS) -c $<
.c.o:
	$(CC) $(CCFLAGS) -c $<
radt:	${src0}
	${LD} ${F90FLAGS} ${src0} ${LIBS} -o ./dyntmat.x

clean:
	rm *.o *.mod
