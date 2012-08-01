
CUDA_INSTALL_PATH ?= /usr/local/cuda

NVCC := $(CUDA_INSTALL_PATH)/bin/nvcc
CFLAGS = -Xcompiler "-Wall -Wextra -fpic"
SOURCES = lmsearch.cu
OBJECTS = lmsearch.o cuseful.o lsfit.o qrdecomp.o

# Includes
INCLUDES += -I. -I/apps/lib64/R/include -I../gputools/src

lmsearch.so: $(OBJECTS)
	$(NVCC) -shared -lcublas $(OBJECTS) -o lmsearch.so

lmsearch.o:
	$(NVCC) $(CFLAGS) $(INCLUDES) -c $(SOURCES)

clean:
	rm lmsearch.o lmsearch.so

