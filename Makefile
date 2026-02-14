HDF5_HOME = /home/denis/gsas2main
CFLAGS = -Wall -Werror -Isrc -I$(HDF5_HOME)/include -fPIC
LDFLAGS = -Lsrc -lbitshuffle -L$(HDF5_HOME)/lib -lhdf5 -lhdf5_hl -lm -Wl,-rpath,$(HDF5_HOME)/lib

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
