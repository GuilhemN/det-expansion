
CC = gcc
CFLAGS = -W -Wall -O2 -march=native -I/usr/include/flint/ -DSQUIRRELS_LEVEL=$(NIST_LEVEL)
LD = gcc
LDFLAGS = 
LIBS = -lm -lgmp -lflint

OBJ1 = build/main.o 

all: build build/main

build:
	-mkdir build

clean:
	-rm -f build/main $(OBJ1)

build/main: $(OBJ1)
	$(LD) $(LDFLAGS) $(LIBS) -o build/main $(OBJ1)

build/main.o: main.c
	$(CC) $(CFLAGS) -c -o build/main.o main.c
