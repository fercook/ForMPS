#
# issue make all
#

ARCH:=$(shell uname)

ifeq ($(ARCH),Darwin)
SYS = MacOSX-x86-64
LAPACK=-framework vecLib
#-L$MKLPATH -I$MKLINCLUDE -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread
BLAS=
FLAGS=-O3 -ffree-line-length-300 -cpp -DTYPEORCLASS="type" -fprofile-arcs -ftest-coverage
RM=rm
FCOMP=gfortran
LINKER=gfortran
LINKFLAGS=
DEBUGFLAG=-g
endif

ifeq ($(ARCH),Linux)
FLAGS = -O3
SYS = Linux-x86-64
LAPACK=
BLAS=
RM=rm
FCOMP=ifort
LINKER=ifort
LINKFLAGS=
endif

BINARIES=_$(SYS)
SOURCES =  constants.f90 error.f90 MatrixHelp.f90 MPSTensor_Class.f90 MatrixProductState_Class.f90
OBJS = $(SOURCES:.f90=.o)
TESTED = MPSTensor_Class MatrixProductState_Class

all: fullmake
obj: object
exec: executable
test: $(TESTED) 
	#testsuite
debug: debug

#--------  HERE START THE USEFUL BITS

object: $(SOURCES)
	$(FCOMP) $(FLAGS) -c $?

fullmake: $(OBJS) main.o
	${LINKER} $(OBJS) main.o ${LINKFLAGS} -L${LIBDIR}  -L${LAPACK} -l${LAPACK} -l${BLAS} -o $@_$(SYS)

%.o: %.f90
	$(FCOMP) $(FLAGS) -c $?

install:
	cp  $(DIR)

clean :
	@ ${RM} -rf *.o *.mod $(BINARIES) *.gcov *.gcda *.gcno
	funit --clean
#-----------------------------------------------------------------------
#
testsuite: 
	
$(TESTED): $(OBJS)
	echo $@
	funit $@ > $@.test 
	gcov $@.f90
	tail -n 5 $@.test
#
