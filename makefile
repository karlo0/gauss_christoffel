CC=gcc
SRC=$(wildcard *.c)
OBJ=$(SRC:.c=.o)
COMPILERFLAGS= -O2 -funroll-all-loops -fomit-frame-pointer -finline-functions -fno-strict-aliasing --param max-inline-insns-single=1800

.PHONY: all, clean

all: $(OBJ)

%.o: %.c
	$(CC) $(COMPILERFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)
