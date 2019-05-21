# This is an commentary line in a makefile
# Start of the makefile
# COMPILER=gfortran -fopenmp

PLATFORM=-DLinux
#PLATFORM=-DWindows

NUMBER=-DNMAX=1300000

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

OUTPUTOPTION=-DnhistoryOutput

COMPILER=ifort -openmp -O2 -fpp $(PLATFORM) $(NUMBER) $(MODEL) $(SEARCH) $(REGION) $(PARALLEL) $(PARTICLE) $(OUTPUTOPTION)

DEMBody: global.o loadFile.o main.o start.o initialSystem.o attitudeQ2M.o initialTriMesh.o initialBondedTriMesh.o \
	initialBiDisperse.o initialLattice.o latticeXY2ID.o loadRestartData.o loadData.o periodic.o periodicBiDisperse.o meshGenerate.o latticeGenerate.o latticeGenerateBiDisperse.o \
	force.o forceLattice.o forcePBC.o forceContactWalls.o forceMovingWalls.o forceBondedWalls.o forceTriMeshWalls.o forceBondedTriMeshWalls.o \
	forceFunnelWalls.o forceGravBody.o forceSphereBody.o forceBiDisperse.o forceBiDisperseMirror.o forceBiDisperseInnerMirror.o forceBiDisperseTraverse.o forceBiDisperseContactWalls.o forceBiDisperseMovingWalls.o forcePlanet.o forceRotSystem.o forceGravTriMesh.o forceParticleLattice.o forceTraverse.o \
	output.o intgrt.o attitude.o intgrtBondedWalls.o intgrtGravBody.o intgrtSphereBody.o intgrtBiDisperse.o
	$(COMPILER)  -o DEMBody global.o loadFile.o main.o start.o initialSystem.o attitudeQ2M.o initialTriMesh.o initialBondedTriMesh.o \
	initialBiDisperse.o initialLattice.o latticeXY2ID.o loadRestartData.o loadData.o periodic.o periodicBiDisperse.o meshGenerate.o latticeGenerate.o latticeGenerateBiDisperse.o \
	force.o forceLattice.o forcePBC.o forceContactWalls.o forceMovingWalls.o forceBondedWalls.o forceTriMeshWalls.o forceBondedTriMeshWalls.o \
	forceFunnelWalls.o forceGravBody.o forceSphereBody.o forceBiDisperse.o forceBiDisperseMirror.o forceBiDisperseInnerMirror.o forceBiDisperseTraverse.o forceBiDisperseContactWalls.o forceBiDisperseMovingWalls.o forcePlanet.o forceRotSystem.o forceGravTriMesh.o forceParticleLattice.o forceTraverse.o \
	output.o intgrt.o attitude.o intgrtBondedWalls.o intgrtGravBody.o intgrtSphereBody.o intgrtBiDisperse.o

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
attitudeQ2M.o: attitudeQ2M.f90
	$(COMPILER)  -c attitudeQ2M.f90
initialTriMesh.o: global.mod initialTriMesh.f90
	$(COMPILER)  -c initialTriMesh.f90
initialBondedTriMesh.o: global.mod initialBondedTriMesh.f90
	$(COMPILER)  -c initialBondedTriMesh.f90
initialBiDisperse.o: global.mod initialBiDisperse.f90
	$(COMPILER)  -c initialBiDisperse.f90
initialLattice.o: global.mod initialLattice.f90
	$(COMPILER)  -c initialLattice.f90
latticeXY2ID.o: latticeXY2ID.f90
	$(COMPILER)  -c latticeXY2ID.f90
loadRestartData.o: global.mod loadfile.mod loadRestartData.f90
	$(COMPILER)  -c loadRestartData.f90
loadData.o: global.mod loadData.f90
	$(COMPILER)  -c loadData.f90
periodic.o: global.mod periodic.f90
	$(COMPILER)  -c periodic.f90
periodicBiDisperse.o: global.mod periodicBiDisperse.f90
	$(COMPILER)  -c periodicBiDisperse.f90
meshGenerate.o: global.mod meshGenerate.f90
	$(COMPILER)  -c meshGenerate.f90
latticeGenerate.o: global.mod latticeGenerate.f90
	$(COMPILER)  -c latticeGenerate.f90
latticeGenerateBiDisperse.o: global.mod latticeGenerateBiDisperse.f90
	$(COMPILER)  -c latticeGenerateBiDisperse.f90
force.o: global.mod force.f90
	$(COMPILER)  -c force.f90
forceLattice.o: global.mod forceLattice.f90
	$(COMPILER)  -c forceLattice.f90
forcePBC.o: global.mod forcePBC.f90
	$(COMPILER)  -c forcePBC.f90
forceContactWalls.o: global.mod forceContactWalls.f90
	$(COMPILER)  -c forceContactWalls.f90
forceMovingWalls.o: global.mod forceMovingWalls.f90
	$(COMPILER)  -c forceMovingWalls.f90
forceBondedWalls.o: global.mod forceBondedWalls.f90
	$(COMPILER)  -c forceBondedWalls.f90
forceTriMeshWalls.o: global.mod forceTriMeshWalls.f90
	$(COMPILER)  -c forceTriMeshWalls.f90
forceBondedTriMeshWalls.o: global.mod forceBondedTriMeshWalls.f90
	$(COMPILER)  -c forceBondedTriMeshWalls.f90
forceFunnelWalls.o: global.mod forceFunnelWalls.f90
	$(COMPILER)  -c forceFunnelWalls.f90
forceGravBody.o: global.mod forceGravBody.f90
	$(COMPILER)  -c forceGravBody.f90
forceSphereBody.o: global.mod forceSphereBody.f90
	$(COMPILER)  -c forceSphereBody.f90
forceBiDisperse.o: global.mod forceBiDisperse.f90
	$(COMPILER)  -c forceBiDisperse.f90
forceBiDisperseMirror.o: global.mod forceBiDisperseMirror.f90
	$(COMPILER)  -c forceBiDisperseMirror.f90
forceBiDisperseInnerMirror.o: global.mod forceBiDisperseInnerMirror.f90
	$(COMPILER)  -c forceBiDisperseInnerMirror.f90
forceBiDisperseTraverse.o: global.mod forceBiDisperseTraverse.f90
	$(COMPILER)  -c forceBiDisperseTraverse.f90
forceBiDisperseContactWalls.o: global.mod forceBiDisperseContactWalls.f90
	$(COMPILER)  -c forceBiDisperseContactWalls.f90
forceBiDisperseMovingWalls.o: global.mod forceBiDisperseMovingWalls.f90
	$(COMPILER)  -c forceBiDisperseMovingWalls.f90
forcePlanet.o: global.mod forcePlanet.f90
	$(COMPILER)  -c forcePlanet.f90
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
attitude.o: global.mod attitude.f90
	$(COMPILER)  -c attitude.f90
intgrtBondedWalls.o: global.mod intgrtBondedWalls.f90
	$(COMPILER)  -c intgrtBondedWalls.f90
intgrtGravBody.o: global.mod intgrtGravBody.f90
	$(COMPILER)  -c intgrtGravBody.f90
intgrtSphereBody.o: global.mod intgrtSphereBody.f90
	$(COMPILER)  -c intgrtSphereBody.f90
intgrtBiDisperse.o: global.mod intgrtBiDisperse.f90
	$(COMPILER)  -c intgrtBiDisperse.f90

clean:
	rm global.mod loadfile.mod global.o loadFile.o main.o start.o initialSystem.o attitudeQ2M.o initialTriMesh.o initialBondedTriMesh.o \
	initialBiDisperse.o initialLattice.o latticeXY2ID.o loadRestartData.o loadData.o periodic.o periodicBiDisperse.o meshGenerate.o latticeGenerate.o latticeGenerateBiDisperse.o \
	force.o forceLattice.o forcePBC.o forceContactWalls.o forceMovingWalls.o forceBondedWalls.o forceTriMeshWalls.o forceBondedTriMeshWalls.o \
	forceFunnelWalls.o forceGravBody.o forceSphereBody.o forceBiDisperse.o forceBiDisperseMirror.o forceBiDisperseInnerMirror.o forceBiDisperseTraverse.o forceBiDisperseContactWalls.o forceBiDisperseMovingWalls.o forcePlanet.o forceRotSystem.o forceGravTriMesh.o forceParticleLattice.o forceTraverse.o \
	output.o intgrt.o attitude.o intgrtBondedWalls.o intgrtGravBody.o intgrtSphereBody.o intgrtBiDisperse.o
# End of the makefile
