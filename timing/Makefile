
CUDA_INSTALL_PATH ?= /usr/local/cuda

NVCC := $(CUDA_INSTALL_PATH)/bin/nvcc
CFLAGS = -Xcompiler "-Wall -Wextra -fpic"
SOURCES = devreset.cu \
	  timing.r
OBJECTS = devreset.o \
	  devreset.so
RCMD = R CMD BATCH --vanilla

# Includes
INCLUDES += -I. -I/apps/lib64/R/include

timing.r.Rout: $(OBJECTS)
	nohup $(RCMD) timing.r &

devreset.so: devreset.o
	$(NVCC) -shared devreset.o -o devreset.so

devreset.o:
	$(NVCC) $(CFLAGS) $(INCLUDES) -c devreset.cu

clean:
	rm devreset.o devreset.so timing.r.Rout nohup.out
