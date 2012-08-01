
CUDA_INSTALL_PATH ?= /usr/local/cuda

NVCC := $(CUDA_INSTALL_PATH)/bin/nvcc
CFLAGS = -Xcompiler "-Wall -Wextra -fpic"
SOURCES = smartlmsearch.cu
OBJECTS = smartlmsearch.o cuseful.o

# Includes
INCLUDES += -I. -I/apps/lib64/R/include -I../gputools/src

smartlmsearch.so: $(OBJECTS)
	$(NVCC) -shared -lcublas $(OBJECTS) -o smartlmsearch.so

smartlmsearch.o:
	$(NVCC) $(CFLAGS) $(INCLUDES) -c $(SOURCES)


clean:
	rm smartlmsearch.o smartlmsearch.so

