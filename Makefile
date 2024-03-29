TARGET = run_tests
LIBS = -lm -lgmp
CC = gcc
CFLAGS = -g -Wall -std=gnu99

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

test: $(TARGET)
	./run_tests

clean:
	-rm -f *.o
	-rm -f $(TARGET)
