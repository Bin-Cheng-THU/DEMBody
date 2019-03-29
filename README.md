# DEMBody

DEMBody is Open Source software for Granular Dynamics in Astronomy. The software focus on investigation of geological features, surface evolution and in-situ exploration on celestrial bodies, e.g., size sorting/segregation on asteroid, mass creep/wasting during geological processes and locomotion dynamics of asteroid lander on granular layer.

The software implements the Discrete Element Method (DEM), coupled with NBody gravity simulation. DEMBody is designed for execution on parallel supercomputers, clusters or multi-core PCs running a Linux (or Windows)-based operating system. 

## Papers based on DEMBody

* [Numerical simulations of the controlled motion of a hopping asteroid lander on the regolith surface](https://doi.org/10.1093/mnras/stz633) - *MNRAS*, 2019
* [Collision-based understanding of the force law in granular impact dynamics](https://doi.org/10.1103/PhysRevE.98.012901) - *PRE*, 2018
* [Asteroid surface impact sampling: dependence of the cavity morphology and collected mass on projectile shape](https://doi.org/10.1038/s41598-017-10681-8) - *Scientific Reports*, 2017

## Getting Started

These instructions will get you a copy of the project for development and testing purposes.

### Structure
  `Data`: storage for data file including point data, wall data and biDisperse data.
  
  `src`: source code of DEMBody.
  
  `Input`: input files including systemControl.dembody, input_points.txt, bondedWallPoint.vtk, trimeshWall.mesh, bondedTriMeshWall.mesh, largeParticles.bidisperse, gravTriMesh.force.

### Pre-/Post-process

We developed [DEMTool](https://github.com/Bin-Cheng-THU/DEMTool.git) for pre- and post- process for DEMBody, e.g., point/mesh file generation and data rending based on [POV-Ray](http://www.povray.org/). Please read the instruction for that.

## Authors

* **Bin Cheng** - *Initial work* - [Bin-Cheng-THU](https://github.com/Bin-Cheng-THU)

See also the list of [contributors](AUTHORS.md) who participated in this project.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Prof. Baoyin and colleagues in LAD
* My girlfriend Fanbing Zeng
* My parents, brother and whole family