# $@ - Target name
# $^ - Target dependencies

CC = gcc
CFLAGS = -std=c99 -Wall -Wpedantic -Wextra

2-1: 2-1.c
	$(CC) $^ $(CFLAGS) -fsanitize=address,undefined -O3 -o $@.o

2-3: 2-3.c
	$(CC) $^ $(CFLAGS) -fsanitize=address,undefined -O3 -o $@.o

2-1_val: 2-1.c
	$(CC) $^ $(CFLAGS) -O3 -o $@.o

2-3_val: 2-3.c
	$(CC) $^ $(CFLAGS) -O3 -o $@.o

debug_2-1: 2-1.c
	$(CC) $^ $(CFLAGS) -fsanitize=address,undefined -Og -g -o $@.o

debug_2-3: 2-3.c
	$(CC) $^ $(CFLAGS) -fsanitize=address,undefined -Og -g -o $@.o
