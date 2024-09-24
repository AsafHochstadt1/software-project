CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm

spkmeans: 
	$(CC) spkmeans.c -D_GNU_SOURCE -o spkmeans $(CFLAGS)