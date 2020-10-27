# Makefile for the atom ray tracing library

CC = gcc  # C compiler
LIBS = -lm
INC = -I../mtwister
CFLAGS = -c -Wall -Wextra -pedantic -O3 -ffast-math # C flags
RM = rm -f   # rm command
TARGET = ../obj/atom_ray_tracing3D.o # target lib
SRCS = atom_ray_tracing3D.c # source files

$(TARGET): $(SRCS)
	$(CC) ${CFLAGS} ${INC} ${LIBS} -o ${TARGET} ${SRCS}

.PHONY: clean
clean:
	-${RM} ${TARGET}
