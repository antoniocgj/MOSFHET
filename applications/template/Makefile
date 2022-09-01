MAKEFILE := $(shell git rev-parse --show-toplevel)/Makefile.def
include $(MAKEFILE)
INCLUDE_FLAGS = $(addprefix -I../../, $(INCLUDE_DIRS))
FLAGS += $(INCLUDE_FLAGS)

all: main

main: $(addprefix ../../src/, $(SRC)) main.c
	$(CC) -g -o main $^ $(OPT_FLAGS) $(LIBS)
