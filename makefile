#
# issue make all or make test
#

ARCH:=$(shell uname)

ifeq ($(ARCH),Darwin)
SYS = MacOSX-x86-64
LAPACK=-framework vecLib
#-L$MKLPATH -I$MKLINCLUDE -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread
BLAS=
FLAGS=$(LAPACK) -Wl,-stack_size -Wl,0x40000000
RM=rm
FCOMP=ifort
LINKER=ifort
LINKFLAGS=
DEBUGFLAG=-profile-functions -profile-loops=all
#-g -debug -save-temps -profile-functions
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
BASICSOURCES = constants.f90 error.f90 Tensor_Class.f90 Operator_Class.f90
MPSSOURCES = MPSTensor_Class.f90 MPOTensor_Class.f90 MPS_Class.f90 MPO_Class.f90 Multiplicator_Class.f90 MPSAlgorithms_Class.f90
PEPSSOURCES = PEPSTensor_Class.f90 PEPOTensor_Class.f90 PEPS_Class.f90 PEPO_Class.f90  Multiplicator2D_Class.f90 PEPSAlgorithms_Class.f90
SOURCES = $(BASICSOURCES) $(MPSSOURCES) $(PEPSSOURCES)
OBJS = $(SOURCES:.f90=.o)
TESTED = Tensor MPSTensor MPOTensor MPS MPO Multiplicator PEPSTensor PEPOTensor PEPS PEPO Multiplicator2D MPSAlgorithms PEPSAlgorithms

all: fullmake
ising: IsingTest
obj: object
exec: executable
test: $(TESTED) 
debug: debug

#--------  HERE START THE USEFUL BITS

object: $(SOURCES)
	$(FCOMP) $(FLAGS) $(DEBUGFLAG) -c $?

fullmake: $(OBJS) main.o
	${LINKER} $(OBJS) main.o ${LINKFLAGS} -L${LIBDIR}  -L${LAPACK} -l${LAPACK} -l${BLAS} -o $@_$(SYS)

%.o: %.f90
	$(FCOMP) $(FLAGS) $(DEBUGFLAG) -c $?

install:
	cp  $(DIR)

IsingTest: $(OBJS) Ising_helper.f90 Ising_tester.f90
	$(FCOMP) $(FLAGS) $(DEBUGFLAG) -c Ising_helper.f90
	$(FCOMP) $(FLAGS) $(DEBUGFLAG) $(OBJS) Ising_helper.o Ising_tester.f90 ${LINKFLAGS} -o IsingTest

clean:
	@ ${RM} -rf *.o *.mod $(BINARIES) *.gcov *.gcda *.gcno *.dyn profiles/*
	funit --clean
#-----------------------------------------------------------------------
#
testsuite: 
	
$(TESTED): $(OBJS)
	$(FCOMP) $(DEBUGFLAG) -c $@.helper.f90
	funit $@_Class > TESTS/$@.test 
	#gcov $@.f90
	#./checkcoverage.sh
	tail -n 5 TESTS/$@.test

coverage:
	./checkcoverage.sh

