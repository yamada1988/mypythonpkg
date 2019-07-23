FC = gfortran

TARGET = libljfort.so
SRC = lj.f90

all:libljfort.so


#object
libljfort.so:lj.f90
	$(FC) -fPIC -shared $< -lm -o $@ -O3
clean:
	rm -f libljfort.so
