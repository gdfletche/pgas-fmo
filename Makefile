#  gdf 8 sep 2020    PGAS-FMO proxy  

 FC=mpifort     # my mac
#FC=mpif90      # summit
#FC=mpif90      # cooley
#F77=mpif77     # cooley
#FC=ftn         # theta
#F77=ftn        # theta
 FLAGS=-O3 -w

all: pgas-fmo.x  

pgas-fmo.x: pgas-fmo.0.0.f90  
	$(FC) $(FLAGS) -o pgas-fmo.x pgas-fmo.0.0.f90 

clean:
	rm *.o; rm *.mod


