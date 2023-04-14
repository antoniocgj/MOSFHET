MOSFHET_DIR = ../../
include $(MOSFHET_DIR)/Makefile.def

all: main

main: $(SRC_MOSFHET) main.c
	$(CC) -g -o main $^ $(OPT_FLAGS) $(LIBS)
