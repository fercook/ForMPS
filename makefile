#
# issue make all
#

ARCH:=$(shell uname)

ifeq ($(ARCH),Darwin)
SYS = MacOSX-x86-64
LAPACK=
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
SOURCES =  constants.f90 error.f90 tensor.f90 
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

#-----------------------------------------------------------------------
#
testsuite: $(OBJS) testsuite.f90 testhelp.o
	$(FCOMP) $(FLAGS) $(OBJS) testhelp.o testsuite.f90 -o Testsuite
	./Testsuite > test.output
#
