CC = gcc
HDF5_INCLUDE = /usr/include/hdf5/serial
CFLAGS = -Wall -Werror -I$(HDF5_INCLUDE) -fPIC
SRCS = bshuf_h5filter.c bitshuffle.c bitshuffle_core.c iochain.c lz4.c
OBJS = $(SRCS:.c=.o)
TARGET = libbitshuffle.a

all: $(TARGET)

$(TARGET): $(OBJS)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
