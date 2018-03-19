# This is an commentary line in a makefile
# Start of the makefile
# COMPILER=gfortran -fopenmp

NUMBER=-DNMAX=5000

MESH=-DCASE1=1 -DCASE2=-65 -DCASE3=-64 -DCASE4=-63 -DCASE5=-4096 -DCASE6=-4095 -DCASE7=-4033 -DCASE8=-4032 -DCASE9=-4031 -DCASE10=-4097 -DCASE11=-4159 -DCASE12=-4160 -DCASE13=-4161

CASE=-DALL
#CASE=-DPROJECTILE

#MODEL=-DHertzMindlinVisco
MODEL=-DHertzMindlinResti

COMPILER=ifort -openmp -O2 -fpp $(NUMBER) $(MESH) $(CASE) $(MODEL) 

DEMBody: global.o loadFile.o main.o start.o zero.o loadRestartData.o input.o attitudeBondedWallsQM.o periodic.o meshGenerate.o force.o forceMirror.o forceContactWalls.o forceBondedWalls.o forceFunnelWalls.o planet.o output.o intgrt.o attitudeBondedWalls.o attitude.o
	$(COMPILER)  -o DEMBody global.o loadFile.o main.o start.o zero.o loadRestartData.o input.o attitudeBondedWallsQM.o periodic.o meshGenerate.o force.o forceMirror.o forceContactWalls.o forceBondedWalls.o forceFunnelWalls.o planet.o output.o intgrt.o attitudeBondedWalls.o attitude.o

global.mod: global.o global.f90
	$(COMPILER)  -c global.f90
global.o: global.f90
	$(COMPILER)  -c global.f90
loadfile.mod: global.mod loadFile.o loadFile.f90
	$(COMPILER)  -c loadFile.f90
loadFile.o: global.mod loadFile.f90
	$(COMPILER)  -c loadFile.f90
main.o: global.mod main.f90
	$(COMPILER)  -c main.f90
start.o: global.mod start.f90
	$(COMPILER)  -c start.f90
zero.o: global.mod zero.f90
	$(COMPILER)  -c zero.f90
loadRestartData.o: global.mod loadfile.mod loadRestartData.f90
	$(COMPILER)  -c loadRestartData.f90
input.o: global.mod loadfile.mod input.f90
	$(COMPILER)  -c input.f90
attitudeBondedWallsQM.o: global.mod attitudeBondedWallsQM.f90
	$(COMPILER)  -c attitudeBondedWallsQM.f90
periodic.o: global.mod periodic.f90
	$(COMPILER)  -c periodic.f90
meshGenerate.o: global.mod meshGenerate.f90
	$(COMPILER)  -c meshGenerate.f90
force.o: global.mod force.f90
	$(COMPILER)  -c force.f90
forceMirror.o: global.mod forceMirror.f90
	$(COMPILER)  -c forceMirror.f90
forceContactWalls.o: global.mod forceContactWalls.f90
	$(COMPILER)  -c forceContactWalls.f90
forceBondedWalls.o: global.mod forceBondedWalls.f90
	$(COMPILER)  -c forceBondedWalls.f90
forceFunnelWalls.o: global.mod forceFunnelWalls.f90
	$(COMPILER)  -c forceFunnelWalls.f90
planet.o: global.mod planet.f90
	$(COMPILER)  -c planet.f90
output.o: global.mod output.f90
	$(COMPILER)  -c output.f90
intgrt.o: global.mod intgrt.f90
	$(COMPILER)  -c intgrt.f90
attitudeBondedWalls.o: global.mod attitudeBondedWalls.f90
	$(COMPILER)  -c attitudeBondedWalls.f90
attitude.o: global.mod attitude.f90
	$(COMPILER)  -c attitude.f90

clean:
	rm global.mod loadfile.mod global.o loadFile.o main.o start.o zero.o loadRestartData.o input.o attitudeBondedWallsQM.o periodic.o meshGenerate.o force.o forceMirror.o forceContactWalls.o forceBondedWalls.o forceFunnelWalls.o planet.o output.o intgrt.o attitudeBondedWalls.o attitude.o
# End of the makefile
