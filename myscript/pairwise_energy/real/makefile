FC = gfortran

TARGET = librealfort.so
SRC = real.f90

all:librealfort.so


#object
librealfort.so:real.f90
	$(FC) -fPIC -shared $< -lm -o $@ -O3
clean:
	rm -f librealfort.so
