# Build:
all:
	gfortran\
		glob_var.f03\
		gquadrature.f03\
		eos.f03\
		reconstruction_mm6eq.f03 reconstruction_soldyn.f03\
		rhs_flux_mm6eq.f03 rhs_soldyn.f03\
		preprocess.f03\
		time_integration.f03\
		main.f03\
		-o multimat1d.exe\
		-ffpe-trap=zero,invalid,overflow\
		-freal-4-real-8\
		-fcheck=all -g

# Remove .mod's and .exe's from the parent directory:
clean:
	rm *.mod *.exe
	rm -r *.dSYM
