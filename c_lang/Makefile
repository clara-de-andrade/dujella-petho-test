CC = gcc
CFLAGS = -Wall -Wextra 
INCLUDES = -I./
LIBS = -lm -lgmp -lmpfr

SRCS = base.c \
       tests.c \
       main.c
OBJS = $(SRCS:.c=.o)
EXEC = main

.PHONY: build clean

$(EXEC): $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) $(SRCS) -o $(EXEC) $(LIBS)

run: $(EXEC)
	./$(EXEC)

clean:
	rm -f $(EXEC)
