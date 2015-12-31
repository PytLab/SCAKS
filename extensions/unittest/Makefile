CC = gcc
INC = -I /usr/local/include/CUnit -I ../include/
LIB = -L /usr/local/lib
FUNCDIR = ../src

test.x: $(FUNCDIR)/kmc_functions.c test_kmc_functions.c
	$(CC) -o $@ $(INC) $(LIB) $^ -lcunit

clean:
	rm -f test.x