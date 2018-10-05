#include <stdio.h>
#include <wrap/io_trimesh/export_ply.h>
#include "Definitions.h"
#include "Volume.h"
#include "Walker.h"
#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/extended_marching_cubes.h>

int main(int argc, char *argv[])
{
	BoundingBox	bbox(vcg::Point3i(-20, -20, -20), vcg::Point3i(20, 20, 20));
	vcg::Point3i			resolution(40, 40, 40);
	
	Volume	volume;
	Walker	walker(bbox, resolution);

	typedef vcg::tri::MarchingCubes<Mesh, Walker>					MarchingCubes;
	typedef vcg::tri::ExtendedMarchingCubes<Mesh, Walker> ExtendedMarchingCubes;

	
	
	// MARCHING CUBES
	Mesh		mc_mesh;
	printf("[MARCHING CUBES] Building mesh...");
	MarchingCubes					mc(mc_mesh, walker);
	walker.BuildMesh<MarchingCubes>(mc_mesh, volume, mc);
	vcg::tri::io::ExporterPLY<Mesh>::Save( mc_mesh, "marching_cubes.ply");
	printf("OK!\n");

	// EXTENDED MARCHING CUBES
	Mesh		emc_mesh;
	printf("[EXTENDED MARCHING CUBES] Building mesh...");
	ExtendedMarchingCubes emc(emc_mesh, walker, 30);
	walker.BuildMesh<ExtendedMarchingCubes>(emc_mesh, volume, emc);
	vcg::tri::io::ExporterPLY<Mesh>::Save( emc_mesh, "extended_marching_cubes.ply");
	printf("OK!\n");
};