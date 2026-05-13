# Compiler
FC = ifx

# Flags - add -g for debugging if needed
FFLAGS = -c -g
LDFLAGS =

# Executable name
TARGET = disp

# Source files
MAIN_SRC = /home/lfrisoli/MuST/project/testSupercell.f90
MODULE_SRC = /home/lfrisoli/MuST/MST/src/SupercellModule.f90

# Object files (dependencies that are pre-compiled libraries)
OBJ = KindParamModule.o \
      MathParamModule.o \
      PhysParamModule.o \
      ErrorHandlerModule.o \
      ChemElementModule.o \
      breakLine.o

# All object files including our compiled sources
ALLOBJ = $(OBJ) SupercellModule.o testSupercell.o

# Default target
all: $(TARGET)

# Link step - include all object files
$(TARGET): $(ALLOBJ)
	$(FC) -o $(TARGET) $(ALLOBJ)

# Compile the module
SupercellModule.o: $(MODULE_SRC)
	$(FC) $(FFLAGS) $(MODULE_SRC)

# Compile the main program
testSupercell.o: $(MAIN_SRC) SupercellModule.o
	$(FC) $(FFLAGS) $(MAIN_SRC)

