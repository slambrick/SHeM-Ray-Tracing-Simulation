# Makefile for source tracing

CC = gcc
INC = -I../mtwister -I../atom_ray_tracing_library
CFLAGS = -Wall -pedantic -Wextra -std=c99
LIBS = -lm
RM = rm -f
TARGET = bin/source_tracing
SRCS = src/source_tracing.c ../obj/atom_ray_tracing3D.o ../obj/mtwister.o

$(TARGET): $(SRCS)
	$(CC) $(INC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LIBS)

clean:
	rm bin/source_tracing
