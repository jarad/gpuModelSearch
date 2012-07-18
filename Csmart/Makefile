all:
	nvcc -g -G -I/apps/lib64/R/include -I../gputools/src -Xcompiler "-Wall -Wextra -fpic" -c smartlmsearch.cu -o smartlmsearch.o
	nvcc -shared -lcublas smartlmsearch.o cuseful.o  -o smartlmsearch.so

clean:
	rm smartlmsearch.o

