HDF5_INCLUDE = /usr/include/hdf5/serial
HDF5_LIB = /usr/lib/x86_64-linux-gnu/hdf5/serial
CFLAGS = -Wall -Werror -Isrc -I$(HDF5_INCLUDE) -fPIC
LDFLAGS = -Lsrc -lbitshuffle -L$(HDF5_LIB) -lhdf5_serial -lhdf5_serial_hl -lm -Wl,-rpath,$(HDF5_LIB)

SRCS = main.c
TARGET = sporta
LIB_TARGET = sporta.so

all: $(TARGET) $(LIB_TARGET)

$(TARGET): $(SRCS)
	$(MAKE) -C src
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

$(LIB_TARGET): $(SRCS)
	$(MAKE) -C src
	$(CC) $(CFLAGS) -shared -o $(LIB_TARGET) $(SRCS) $(LDFLAGS)

clean:
	$(MAKE) -C src clean
	rm -f $(TARGET) $(LIB_TARGET)
