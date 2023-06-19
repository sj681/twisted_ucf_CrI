#===================================================================
#
#								Makefile for VAMPIRE
#
#===================================================================

# Specify compiler for MPI compilation with openmpi
export OMPI_CXX=g++ -std=c++11

# Compiler
GCC=g++ -std=c++11 -DCOMP='"GNU C++ Compiler"'


LIBS= -lstdc++
#-lm $(FFTLIBS) -L/opt/local/lib/

GCC_CFLAGS=-I./hdr -O0
GCC_LDFLAGS= -lstdc++ -I./hdr
# Objects
OBJECTS= \
obj/main.o \
obj/exchange/exchange.o \
obj/initialise/initialise.o \
obj/positions/positions.o \
obj/read_in_files.o \


EXECUTABLE=main.o

# Set default make target in GNU make > v3.81
.DEFAULT_GOAL := all

# make serial and parallel versions and utilities
all:  $(OBJECTS)
	$(GCC) $(GCC_LDFLAGS)  $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $(OPTIONS) $<

clean:
	@rm -f obj/*.o
	@rm -f obj/*/*.o

purge:
	@rm -f obj/*.o
	@rm -f obj/*/*.o
	@rm -f vampire-*

tidy:
	@rm -f *~
	@rm -f hdr/*~
	@rm -f src/*~
	@rm -f src/*/*~

tests:
	$(MAKE) -C test/integration/
	$(MAKE) -C test/unit/
