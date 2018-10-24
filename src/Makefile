# This is an commentary line in a makefile
# Start of the makefile
# COMPILER=gfortran -fopenmp

NUMBER=-DNMAX=1610000

#MODEL=-DHertzMindlinVisco
MODEL=-DHertzMindlinResti

#SEARCH=-DTraverseSearch
SEARCH=-DLatticeSearch -DParticleLattice
#SEARCH=-DLatticeSearch

REGION=-DConfined
#REGION=-DnConfined

#PARALLEL=-Dnopenmp
PARALLEL=-Dopenmp

PARTICLE=-DPP=1

OUTPUTOPTION=-DhistoryOutput

COMPILER=ifort -openmp -O2 -fpp $(NUMBER) $(MODEL) $(SEARCH) $(REGION) $(PARALLEL) $(PARTICLE) $(OUTPUTOPTION)

DEMBody: global.o loadFile.o main.o start.o initialSystem.o attitudeBondedWallsQM.o loadRestartData.o loadData.o periodic.o meshGenerate.o latticeGenerate.o force.o forceLattice.o forceMirror.o forceContactWalls.o forceBondedWalls.o forceTriMeshWalls.o forceFunnelWalls.o forceGravBody.o planet.o forceRotSystem.o forceGravTriMesh.o forceParticleLattice.o forceTraverse.o output.o intgrt.o attitudeBondedWalls.o attitudeGravBody.o attitude.o
	$(COMPILER)  -o DEMBody global.o loadFile.o main.o start.o initialSystem.o attitudeBondedWallsQM.o loadRestartData.o loadData.o periodic.o meshGenerate.o latticeGenerate.o force.o forceLattice.o forceMirror.o forceContactWalls.o forceBondedWalls.o forceTriMeshWalls.o forceFunnelWalls.o forceGravBody.o planet.o forceRotSystem.o forceGravTriMesh.o forceParticleLattice.o forceTraverse.o output.o intgrt.o attitudeBondedWalls.o attitudeGravBody.o attitude.o

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
initialSystem.o: global.mod loadfile.mod initialSystem.f90
	$(COMPILER)  -c initialSystem.f90
attitudeBondedWallsQM.o: global.mod attitudeBondedWallsQM.f90
	$(COMPILER)  -c attitudeBondedWallsQM.f90
loadRestartData.o: global.mod loadfile.mod loadRestartData.f90
	$(COMPILER)  -c loadRestartData.f90
loadData.o: global.mod loadData.f90
	$(COMPILER)  -c loadData.f90
periodic.o: global.mod periodic.f90
	$(COMPILER)  -c periodic.f90
meshGenerate.o: global.mod meshGenerate.f90
	$(COMPILER)  -c meshGenerate.f90
latticeGenerate.o: global.mod latticeGenerate.f90
	$(COMPILER)  -c latticeGenerate.f90
force.o: global.mod force.f90
	$(COMPILER)  -c force.f90
forceLattice.o: global.mod forceLattice.f90
	$(COMPILER)  -c forceLattice.f90
forceMirror.o: global.mod forceMirror.f90
	$(COMPILER)  -c forceMirror.f90
forceContactWalls.o: global.mod forceContactWalls.f90
	$(COMPILER)  -c forceContactWalls.f90
forceBondedWalls.o: global.mod forceBondedWalls.f90
	$(COMPILER)  -c forceBondedWalls.f90
forceTriMeshWalls.o: global.mod forceTriMeshWalls.f90
	$(COMPILER)  -c forceTriMeshWalls.f90
forceFunnelWalls.o: global.mod forceFunnelWalls.f90
	$(COMPILER)  -c forceFunnelWalls.f90
forceGravBody.o: global.mod forceGravBody.f90
	$(COMPILER)  -c forceGravBody.f90
planet.o: global.mod planet.f90
	$(COMPILER)  -c planet.f90
forceRotSystem.o: global.mod forceRotSystem.f90
	$(COMPILER)  -c forceRotSystem.f90
forceGravTriMesh.o: global.mod forceGravTriMesh.f90
	$(COMPILER)  -c forceGravTriMesh.f90
forceParticleLattice.o: global.mod forceParticleLattice.f90
	$(COMPILER)  -c forceParticleLattice.f90
forceTraverse.o: global.mod forceTraverse.f90
	$(COMPILER)  -c forceTraverse.f90
output.o: global.mod output.f90
	$(COMPILER)  -c output.f90
intgrt.o: global.mod intgrt.f90
	$(COMPILER)  -c intgrt.f90
attitudeBondedWalls.o: global.mod attitudeBondedWalls.f90
	$(COMPILER)  -c attitudeBondedWalls.f90
attitudeGravBody.o: global.mod attitudeGravBody.f90
	$(COMPILER)  -c attitudeGravBody.f90
attitude.o: global.mod attitude.f90
	$(COMPILER)  -c attitude.f90

clean:
	rm global.mod loadfile.mod global.o loadFile.o main.o start.o initialSystem.o attitudeBondedWallsQM.o loadRestartData.o loadData.o periodic.o meshGenerate.o latticeGenerate.o force.o forceLattice.o forceMirror.o forceContactWalls.o forceBondedWalls.o forceTriMeshWalls.o forceFunnelWalls.o forceGravBody.o planet.o forceRotSystem.o forceGravTriMesh.o forceParticleLattice.o forceTraverse.o output.o intgrt.o attitudeBondedWalls.o attitudeGravBody.o attitude.o
# End of the makefile
