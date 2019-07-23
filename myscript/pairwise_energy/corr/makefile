FC = gfortran

TARGET = libcorrfort.so
SRC = corr.f90

all:libcorrfort.so


#object
libcorrfort.so:corr.f90
	$(FC) -fPIC -shared $< -lm -o $@ -O3
clean:
	rm -f libcorrfort.so
