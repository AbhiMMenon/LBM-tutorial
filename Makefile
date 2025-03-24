# Compiler and Flags
FC = gfortran
FFLAGS = -Ofast -Wall -Warray-temporaries -Wfrontend-loop-interchange 

# Directories
SRCDIR = src
OBJDIR = obj
BINDIR = .

# Source files and object files
SRC = $(wildcard $(SRCDIR)/*.f90)
OBJ = $(SRC:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
EXEC = $(BINDIR)/LBMsolver

# Default target
all: $(EXEC)

# Link the program
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJ)

# Compile the source files into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean object and binary files
clean:
	rm -f $(OBJDIR)/*.o $(EXEC)

movie:
	echo "Making movie"
	gnuplot src/plot.gnu | ffmpeg -framerate 120 -f  png_pipe -i pipe: -vf negate -y movie.mp4


# Phony targets (not actual files)
.PHONY: all clean

