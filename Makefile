
CC=g++

CFLAGS=-std=c++11 -O3

objects=LinearCDSfold

linearcdsfold:

		$(CC) $(CFLAGS) src/LinearCDSfold.cpp -o LinearCDSfold 

clean:
	-rm $(objects)
