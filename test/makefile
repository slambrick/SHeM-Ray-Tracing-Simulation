CC = gcc
CFLAGS = -Wall -pedantic -std=c99
LIBS = -lgsl -lgslcblas -lm
SDIR = src
EDIR = ../mexFiles
BDIR = bin
ODIR = obj
EODIR = extobj
MDIR = mex
INC = -I$(EDIR) -I$(HOME)/.local/share/MATLAB/R2019b/extern/include

# source file names of the tests
SRC = $(shell find $(SDIR) -type f -regex '.*\.c')

# external sources -- source files being tested
_ESRC = small_functions3D.c distributions.c ray_tracing_structs3D.c common_helpers.c extract_inputs.c
ESRC = $(patsubst %,$(EDIR)/%,$(_ESRC))

# external headers
EH = $(patsubst %.c,%.h,$(ESRC))

# objects made from external sources
EOBJ = $(patsubst %.c,$(EODIR)/%.o,$(_ESRC))

.PHONY: all remove clean debug
all: $(BDIR)/simple

# create the necessary directories
directories:
	mkdir -p $(ODIR)
	mkdir -p $(EODIR)
	mkdir -p $(BDIR)

# # compile mex files
# %.mexa64: $(MDIR)/%.c $(ESRC)
# 	mex -R2018a $(LIBS) $(INC) CFLAGS="$$CFLAGS $(CFLAGS)" $< $(ESRC)

# compile external objects
$(EODIR)/%.o: $(EDIR)/%.c $(EDIR)/%.h
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

# compile the test objects
$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

# link into executable
$(BDIR)/simple: $(ODIR)/simple.o $(EOBJ)
	$(CC) $(CFLAGS) $(INC) $^ $(LIBS) -o $@

# link to executable
# $(BDIR)/%: $(ODIR)/%.o $(OBJ)
#         $(CC) $(CFLAGS) $^ $(LIBS) -o $@

remove:
	rm -f $(ODIR)/*.o
	rm -f $(BDIR)/*

clean:
	rm -f $(ODIR)/*.o
