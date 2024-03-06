# $@ - Target name
# $^ - Target dependencies

CC = gcc
CFLAGS = -std=c99 -Wall -Wpedantic -Wextra

main_6-1: 6-1.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -O3 -o $@.o

val_6-1: 6-1.c
	$(CC) $^ $(CFLAGS) -lm -O3 -o $@.o

debug_6-1: 6-1.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -Og -g -o $@.o

main_6-2: 6-2.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -O3 -o $@.o

val_6-2: 6-2.c
	$(CC) $^ $(CFLAGS) -lm -O3 -o $@.o

debug_6-2: 6-2.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -Og -g -o $@.o

main_6-3: 6-3.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -O3 -o $@.o

val_6-3: 6-3.c
	$(CC) $^ $(CFLAGS) -lm -O3 -o $@.o

debug_6-3: 6-3.c
	$(CC) $^ $(CFLAGS) -lm -fsanitize=address,undefined -Og -g -o $@.o
