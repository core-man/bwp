CC    = gcc
CFLAG = -Wall

BIN = .

all: bwp clean

bwp: bwp.o sacio.o
	$(CC) ${CFLAGS} -o $(BIN)/$@ $^ -lm

clean:
	rm *.o
