
TARGET := faves
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)
CURRENT_DIR := $(shell pwd)

CC := gcc
CFLAGS := -O3 -mavx2 -msse4.1 -Wall -Wextra -Wpedantic
LDFLAGS := -lm -pthread -lz

# ========================================
#  Build Rules
# ========================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(CXXLIBS) $(LDFLAGS) -L lib -lblend 
	rm -f $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# ========================================
#  Utility Targets
# ========================================

clean: 
	@echo "Cleaning"
	rm -f $(OBJS)
	rm -f $(TARGET)

blend: 
	@mkdir -p lib
	gcc -O3 -mavx2 -msse4.1 -c sketch/blend.c -I blend
	ar rcs lib/libblend.a blend.o
	@rm -f blend.o

install: clean blend $(TARGET)
