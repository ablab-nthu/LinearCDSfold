
CC=g++

CFLAGS=-std=c++11 -O3

objects=main

linearcdsfold:

		$(CC) $(CFLAGS) main.cpp -o main

clean:
	-rm $(objects)
