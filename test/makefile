CC = gcc
INC = -I../mtwister -I../mexFiles
CFLAGS = -Wall -pedantic -std=c99
SDIR = src
EDIR = ../mtwister
EDIR2 = ../mexFiles
BDIR = bin
EODIR = extobj
ODIR = obj
LIBS = -lm
DEPS = ../mtwister/mtwister.h ../mexFiles/common_helpers.h

# Make all the directories
$(shell mkdir -p $(ODIR) $(BDIR) $(EODIR))

# Put the mtwister object into an external object directory
_EOBJ = mtwister.o small_functions3D.o common_helpers.o
EOBJ = $(patsubst %,$(EODIR)/%,$(_EOBJ))
$(EODIR)/%.o: $(EDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(INC)

$(EODIR)/%.o: $(EDIR2)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(INC)

# Put the object files in their own directory
_OBJ = simple.o random.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(INC)

# Link into an executable
$(BDIR)/simple: obj/simple.o $(EOBJ) $(DEPS)
	$(CC) -o $@ $^ $(LIBS) $(INC)

$(BDIR)/random: obj/random.o $(EOBJ) $(DEPS)
	$(CC) -o $@ $^ $(LIBS) $(INC) 

# Makes sure that make doesn't do anything with a file called 'clean'
.PHONY: clean

# Clean up the files that make produces
clean:
	rm -f $(ODIR)/*.o $(BDIR)/* $(EODIR)/*.o
