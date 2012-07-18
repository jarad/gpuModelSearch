
NVCC = nvcc
CFLAGS = -g -G -Xcompiler "-Wall -Wextra -fpic"
SOURCES = lmsearch.cu
OBJECTS = cuseful.o lsfit.o qrdecomp.o

# Includes
INCLUDES += -I. -I/apps/lib64/R/include -I../gputools/src

lmsearch.so: lmsearch.o
	$(NVCC) -shared -lcublas $(OBJECTS)

lmsearch.o:
	$(NVCC) $(CFLAGS) $(INCLUDES) -c $(SOURCES)

clean:
	rm lmsearch.o lmsearch.so

