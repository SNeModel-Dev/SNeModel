all: ShockHeatingPlusAccelv3 specuvot3 SpectToMag
ShockHeatingPlusAccelv3: ShockHeatingPlusAccelv3.f90
	gfortran -g ShockHeatingPlusAccelv3.f90 -o ShockHeatingPlusAccelv3.exe
specuvot3: specuvot3.f90
	gfortran -O2 specuvot3.f90 -o specuvot3.exe
SpectToMag: SpectToMag.f90
	gfortran -O2 SpectToMag.f90 -o SpectToMag.exe
