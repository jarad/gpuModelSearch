
NVCC = nvcc
CFLAGS = -g -G -I/apps/lib64/R/include -I../gputools/src -Xcompiler "-Wall -Wextra -fpic"
SOURCES = lmsearch.cu
OBJECTS = cuseful.o lsfit.o qrdecomp.o


lmsearch.so: lmsearch.o
	$(NVCC) -shared -lcublas $(OBJECTS)

lmsearch.o:
	$(NVCC) $(CFLAGS) -c $(SOURCES)

clean:
	rm lmsearch.o lmsearch.so

