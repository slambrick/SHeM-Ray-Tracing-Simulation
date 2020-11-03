# Makefile for the pure C tests of ray tracing

CC = gcc
INC = -I../mtwister -I../atom_ray_tracing_library
CFLAGS = -Wall -pedantic -Wextra -std=c99
LIBS = -lm
RM = rm -f
TARGET = bin/single_ray
SRCS = src/tests/src/single_experiment_test.c ../obj/atom_ray_tracing3D.o ../obj/mtwister.o

$(TARGET): $(SRCS)
	$(CC) $(INC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LIBS)

clean:
	rm bin/single_experiment_test
