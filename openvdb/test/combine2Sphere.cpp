#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/Grid.h>
#include <openvdb/Types.h>
#include <openvdb/math/Vec3.h>

#define PI 3.1415926535897931

struct Local {
	static inline void dot(const openvdb::Vec3f& a, const openvdb::Vec3f& b, openvdb::Vec3f& result) {
		//result.x() =  acos( a.dot(b) ) / PI * 180;
		result.x() =  openvdb::math::angle(a, b) / PI * 180;
	}
};

int main()
{
    openvdb::initialize();
	
    openvdb::FloatGrid::Ptr grid1 =
        openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
            /*radius=*/2, /*center=*/openvdb::Vec3f(1, 0, 0),
            /*voxel size=*/0.5, /*width=*/4.0);

	openvdb::VectorGrid::Ptr outGrad1 = openvdb::tools::gradient(*grid1);
    outGrad1->setName("gradient1");
	
    openvdb::FloatGrid::Ptr grid2 =
	openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
            /*radius=*/2, /*center=*/openvdb::Vec3f(-1, 0, 0),
            /*voxel size=*/0.5, /*width=*/4.0);
	openvdb::VectorGrid::Ptr outGrad2 = openvdb::tools::gradient(*grid2);
    // Name the grid "LevelSetSphere".
    outGrad2->setName("gradient2");

	
	openvdb::VectorGrid::Ptr result = openvdb::VectorGrid::create();
	
	result->tree().combine2(outGrad1->tree(), outGrad2->tree(), Local::dot, /*prune=*/false); 
	
	openvdb::VectorGrid::Ptr result2 = openvdb::gridPtrCast<openvdb::VectorGrid>(result->deepCopyGrid());
	
    // Create a VDB file object.
    openvdb::io::File file("gradient.vdb");
    // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(result);

    // Write out the contents of the container.
    file.write(grids);
    file.close();
	
	/*
	std::cout << "gradient1 values are : " << std::endl;
	for (openvdb::VectorGrid::ValueOnIter iter = outGrad1->beginValueOn(); iter; ++iter) {
		openvdb::Vec3f value = iter.getValue();
		std::cout << "outGrad1" << iter.getCoord() << " = " << value << std::endl;
	}
	
	std::cout << "gradient2 values are : " << std::endl;
	for (openvdb::VectorGrid::ValueOnIter iter = outGrad2->beginValueOn(); iter; ++iter) {
		openvdb::Vec3f value = iter.getValue();
		std::cout << "outGrad2" << iter.getCoord() << " = " << value << std::endl;
	}
	
	std::cout << "gradient dot values are : " << std::endl;
	for (openvdb::VectorGrid::ValueOnIter iter = result->beginValueOn(); iter; ++iter) {
		openvdb::Vec3f value = iter.getValue();
		std::cout << "Grid" << iter.getCoord() << " = " << value.x() << std::endl;
	}
	*/
  
}