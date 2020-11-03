# Makefile for the mtwister library

CC = gcc  # C compiler
CFLAGS = -c -Wall -Wextra -pedantic -O2 -std=c99  # C flags
RM = rm -f   # rm command
TARGET = ../obj/mtwister.o  # target lib
SRCS = mtwister.c  # source files

$(TARGET): $(SRCS)
	$(CC) ${CFLAGS} -o ${TARGET} ${SRCS}

.PHONY: clean
clean:
	-${RM} ${TARGET}
