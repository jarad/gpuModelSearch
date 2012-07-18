
NVCC = nvcc
SOURCES = lmsearch.cu
OBJECTS = lmsearch.o cuseful.o lsfit.o qrdecomp.o


all:
	$(NVCC) -g -G -I/apps/lib64/R/include -I../gputools/src -Xcompiler "-Wall -Wextra -fpic" -c $(SOURCES) -o lmsearch.o
	$(NVCC) -shared -lcublas $(OBJECTS) -o lmsearch.so

clean:
	rm lmsearch.o

