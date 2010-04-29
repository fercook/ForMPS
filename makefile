#
# issue make all
#

ARCH:=$(shell uname)

ifeq ($(ARCH),Darwin)
SYS = MacOSX-x86-64
LAPACK=-L$MKLPATH -I$MKLINCLUDE -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread
BLAS=
FLAGS= -O3
RM=rm
FCOMP=ifort
LINKER=ifort
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
SOURCES =  constants.f90 error.f90 MPSTensor_Class.f90 
OBJS = $(SOURCES:.f90=.o)

all: fullmake
obj: object
exec: executable
test: testsuite
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
	@ ${RM} -rf *.o *.mod $(BINARIES)
	funit --clean
#-----------------------------------------------------------------------
#
testsuite: 
	funit MPSTensor_Class > test.results
	tail test.results
#
