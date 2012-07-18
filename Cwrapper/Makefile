all:
	nvcc -g -G -I/apps/lib64/R/include -I../gputools/src -Xcompiler "-Wall -Wextra -fpic" -c lmsearch.cu -o lmsearch.o
	nvcc -shared -lcublas lmsearch.o cuseful.o lsfit.o qrdecomp.o -o lmsearch.so

