MOSFHET_DIR = ../../
include $(MOSFHET_DIR)/Makefile.def

all: main

main: $(SRC_MOSFHET) main.c
	$(CC) -g -o main $^ $(OPT_FLAGS) $(LIBS)

vertical_packing: $(SRC_MOSFHET) vertical_packing.c
	$(CC) -g -o vertical_packing $^ $(OPT_FLAGS) $(LIBS)

clean: 
	rm --f main vertical_packing