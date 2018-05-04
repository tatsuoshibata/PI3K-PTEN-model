# Makefile
CC = gcc
CFLAGS = -O3 -I$(HOME)/include
LDFLAGS =   -L$(HOME)/lib -lnrutil

OBJS =  PI3KPTENModel.o

all: PI3KPTENModel.out

PI3KPTENModel.out: $(OBJS)
	$(CC) -o $@ $< $(LDFLAGS)

PI3KPTENModel.o:PI3KPTENModel.c cell.h
	$(CC) $(CFLAGS) -c PI3KPTENModel.c -o $@

.c.o:
	$(CC) $(CFLAGS) -c $<

gc.o: cell.h
# Disp.o:

.PHONY: clean
clean:
	rm -f *.out  *.o
