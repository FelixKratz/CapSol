CXX = g++
CFLAGS = -std=c++11 -g -Wall
LINKER_FLAGS = 
LINKING_GSL = -lm -lgslcblas -lgsl
LINKING_ARMADILLO = -larmadillo

# MacOS libs and include paths. May vary for other system or may not be needed.
LIBS=-L/opt/homebrew/lib
INCLUDE=-I/opt/homebrew/include -I/usr/local/include

FOLDERS = laplace hooke kelvin io misc contact bending
LOCAL_INCLUDES = $(addprefix -I, $(FOLDERS))
ODIR = bin
ODIR_SUB = $(patsubst %, $(ODIR)/%, $(FOLDERS))

HOOKE_OBJ = hooke/hooke.o hooke/hookeSolver.o hooke/shapeEquations.o hooke/pressureShooting.o
KELVIN_OBJ = kelvin/kelvin.o kelvin/kelvinFit.o kelvin/shapeEquations.o kelvin/pressureShooting.o kelvin/kelvinSequence.o
LAPLACE_OBJ = laplace/laplace.o
CONTACT_OBJ = contact/contact.o contact/shapeEquations.o contact/pressureShooting.o contact/volumeShooting.o contact/angleShooting.o
BENDING_OBJ = bending/bending.o bending/shapeEquations.o
INTERFACE_OBJ = io/interface.o io/test.o io/api.o io/main.o
UTILITY_OBJ = misc/pointSet.o

_OBJ = $(HOOKE_OBJ) $(KELVIN_OBJ) $(LAPLACE_OBJ) $(CONTACT_OBJ) $(BENDING_OBJ) $(INTERFACE_OBJ) $(UTILITY_OBJ)
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ)) 

all: clean
	$(MAKE) api -j10

debug: CFLAGS += -DDEBUG -DWATCH_LEAST_SQUARES -DWATCH_KELVIN_FITTING
debug: $(ODIR)/CapSol
	mv $(ODIR)/CapSol $(ODIR)/debug

asan: CFLAGS += -DTEST -DDEBUG -fsanitize=address -fsanitize=undefined
asan: $(ODIR)/CapSol

main: CFLAGS += # -O3
main: $(ODIR)/CapSol
	cd $(ODIR) && ./CapSol

test: CFLAGS += -O3 -DTEST
test: $(ODIR)/CapSol
	cd $(ODIR) && ./CapSol

api: CFLAGS += -DAPI -O3 -fPIC
api: LINKER_FLAGS += -shared -undefined dynamic_lookup 
api: INCLUDE += $(shell python3 -m pybind11 --includes)
api: SUFFIX_PYBIND = $(shell python3-config --extension-suffix)
api: $(ODIR)/CapSol
	mv $(ODIR)/CapSol frontend/CapSol$(SUFFIX_PYBIND)

$(OBJ): $(ODIR)/%.o: %.cpp %.hpp
	@-mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CFLAGS) $(INCLUDE) $(LOCAL_INCLUDES)

$(ODIR)/CapSol: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINKER_FLAGS) $(LIBS) $(LINKING_GSL) $(LINKING_ARMADILLO)

.PHONY: clean
clean:
	@-rm -r $(ODIR)
