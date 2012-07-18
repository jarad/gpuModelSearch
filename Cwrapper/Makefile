
NVCC = nvcc
SOURCES = lmsearch.cu

all:
	$(NVCC) -g -G -I/apps/lib64/R/include -I../gputools/src -Xcompiler "-Wall -Wextra -fpic" -c $(SOURCES) -o lmsearch.o
	$(NVCC) -shared -lcublas lmsearch.o cuseful.o lsfit.o qrdecomp.o -o lmsearch.so

clean:
	rm lmsearch.o

