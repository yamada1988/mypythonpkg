FC = gfortran

TARGET = liburealfort.so
SRC = real.f90

all:liburealfort.so


#object
liburealfort.so:real.f90
	$(FC) -fPIC -fopenmp -shared $< -lm -o $@ -O3
clean:
	rm -f liburealfort.so
