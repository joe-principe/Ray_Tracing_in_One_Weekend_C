# $@ - Target name
# $^ - Target dependencies

CC = gcc
CFLAGS = -std=c99 -Wall -Wpedantic -Wextra

main: main.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -O3 -o $@.o

val: main.c
	$(CC) $^ $(CFLAGS) -lm -O3 -o $@.o

debug: main.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -Og -g -o $@.o
