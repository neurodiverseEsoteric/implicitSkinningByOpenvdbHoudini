#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/GridOperators.h>

int main()
{
    openvdb::initialize();
	
    openvdb::ScalarGrid::Ptr grid =
        openvdb::tools::createLevelSetSphere<openvdb::ScalarGrid>(
            /*radius=*/1, /*center=*/openvdb::Vec3f(0, 0, 0),
            /*voxel size=*/0.2, /*width=*/4.0);

    // Name the grid "LevelSetSphere".
    grid->setName("valueGrid");
	
	openvdb::tools::Gradient<openvdb::ScalarGrid> myGrid(*grid);
	openvdb::VectorGrid::Ptr resultGrid = openvdb::VectorGrid::create();
	resultGrid = myGrid.process();
	resultGrid->setName("gradGrid");
    // Create a VDB file object.
    openvdb::io::File file("sphere.vdb");
    // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;
	grids.push_back(grid);
    grids.push_back(resultGrid);

    // Write out the contents of the container.
    file.write(grids);
    file.close();
}