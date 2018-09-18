# Build:
all:
	gfortran\
		glob_var.f95\
		eos.f95\
		rhs_flux.f95\
		rhs_flux_mm6eq.f95\
		preprocess.f95\
		time_integration.f95\
		main.f95\
		-o multimat1d.exe\
		-ffpe-trap=zero,invalid,overflow\
		-freal-4-real-8\
		-fcheck=all -g

# Remove .mod's and .exe's from the parent directory:
clean:
	rm *.mod *.exe
	rm -r *.dSYM

# Remove data from the output directory:
rmout:
	rm -f ./output/*
