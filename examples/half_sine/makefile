NAME=half_sine
CC=gcc
OP=../../
DEPS=$(wildcard $(OP)*.o)
HPATH=-I$(OP)

LDFLAGS= -lm -llapack -lopenblas
COMPILERFLAGS= -O2 -funroll-all-loops -fomit-frame-pointer -finline-functions -fno-strict-aliasing --param max-inline-insns-single=1800

.PHONY: all, clean_obj, clean_bin, echo_src, echo_obj, clean_out
all:
	$(MAKE) -C $(OP)
	$(MAKE) $(NAME)

$(NAME): $(NAME).o $(DEPS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(COMPILERFLAGS) $(HPATH) -c $< -o $@

echo_src:
	@echo *.c

echo_obj:
	@echo $(DEPS)

clean_obj: 
	rm -rf *.o

clean_bin:
	rm -rf  $(NAME)

clean_out:
	rm -rf *.txt *.log
