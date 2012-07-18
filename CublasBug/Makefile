all:
	nvcc -I/apps/lib64/R/include -Xcompiler "-Wall -Wextra -fpic" -c RcublasTest.cu -o RcublasTest.o
	nvcc -shared -lcublas RcublasTest.o -o RcublasTest.so
	nvcc -I/apps/lib64/R/include -Xcompiler "-Wall -Wextra -fpic" -c RcudaTest.cu -o RcudaTest.o
	nvcc -shared RcudaTest.o -o RcudaTest.so
	nvcc -I/apps/lib64/R/include -Xcompiler "-Wall -Wextra -fpic" -c RcublasV2Test.cu -o RcublasV2Test.o
	nvcc -shared -lcublas RcublasV2Test.o -o RcublasV2Test.so
	nvcc -Xcompiler "-Wall -Wextra" -lcublas cublasV2Test.cu -o cublasV2Test.out
	nvcc -Xcompiler "-Wall -Wextra" -lcublas cublasV2nomemTest.cu -o cublasV2nomemTest.out

clean:
	rm RcublasTest.o RcudaTest.o RcublasV2Test.o RcublasTest.so RcudaTest.so RcublasV2Test.so cublasV2Test.out cublasV2nomemTest.out

