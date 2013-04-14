# Copyright 2012 Shan An <anshan.tju@gmail.com>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.

# C source code

CC = gcc
CXX = g++

CSRC	= io_png/io_png.c colortransfer.c

# all source code
SRC	= $(CSRC) 

# C objects
COBJ	= $(CSRC:.c=.o)
# all objects
OBJ	= $(COBJ) 

# binary target
BIN	= colortransfer

default	: $(BIN)

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops #delete for gdb debugging 

# C compilation flags 
CFLAGS	= $(COPT) -Wall -Wextra -Werror \
	-Wno-write-strings -ansi #-ggdb #for gdb debugging

# link flags
LDFLAGS	= -lpng -lm

# use local embedded libraries with `make LOCAL_LIBS=1`
ifdef LOCAL_LIBS
# library location
LIBDIR = io_png/libs/build/lib
INCDIR = io_png/libs/build/include 
# libpng is required
LIBDEPS += libpng
# compile options to use the local libpng header
CFLAGS 	+= -I$(INCDIR) -D_LOCAL_LIBS
# link options to use the local libraries
LDFLAGS = $(LIBDIR)/libpng.a $(LIBDIR)/libz.a -lm
# io_png.o needs png.h
io_png/io_png.o	:  $(LIBDEPS)
endif

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
LDFLAGS += -lgomp
else
CFLAGS	+= -Wno-unknown-pragmas
endif

# build the local png library
.PHONY	: libpng
libpng	:
	$(MAKE) -C io_png/libs libpng

# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS) -I/opt/local/include/ -I/usr/local/include/   

# link all the opject code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS) -L/opt/local/lib/ -L/usr/local/lib/  -lpng


# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
	$(RM) $(BIN)
distclean	: clean
	$(RM) $(BIN)
