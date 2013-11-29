#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/GridOperators.h>

using namespace openvdb;

float remap(float value, float r)
{
	float ratio1, ratio2;
    if(value < -r) {
		return 1;
    }else if(value > r) {
		return 0;
	}else{
		ratio1 = value / r;
		ratio2 = ratio1 * ratio1;
		return -3.0f / 16 * ratio2 * ratio2 * ratio1 + 5.0f / 8 * ratio2 * ratio1 - 15.0f / 16 * ratio1 + 0.5f;
	}
}

void makeSphere(FloatGrid::Ptr grid, float radius, const CoordBBox& indexBB, double h, float range)
{
  typename FloatGrid::Accessor accessor = grid->getAccessor();
  range = range * h;
  for (Int32 i = indexBB.min().x(); i <= indexBB.max().x(); ++i) {
    for (Int32 j = indexBB.min().y(); j <= indexBB.max().y(); ++j) {
      for (Int32 k = indexBB.min().z(); k <= indexBB.max().z(); ++k) {
        // transform point (i, j, k) of index space into world space
        Vec3d p(i * h, j * h, k * h);
        // compute level set function value
        float distance = sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z()) - radius;
		
		distance = remap(distance, range);
		
        accessor.setValue(Coord(i, j, k), distance);
      }
    }
  }	

  grid->setTransform(openvdb::math::Transform::createLinearTransform(h));
}

int main()
{
  openvdb::initialize();

  openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);

  float radius = 1;
  float voxelSize = 0.2;
  float bboxSize = 4 / voxelSize;
  CoordBBox indexBB(Coord(-bboxSize, -bboxSize, -bboxSize), Coord(bboxSize, bboxSize, bboxSize));
  float range = 10;
  makeSphere(grid, radius, indexBB, voxelSize, range);

  // specify dataset name
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
  
//   const float outside = grid->background();
//   std::cout << "background value is : " << outside << std::endl;
//   
//   std::cout << "grid values are : " << std::endl;
//   for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter) {
//     std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
//   }
//   
//   openvdb::tools::Gradient<openvdb::ScalarGrid> myGrid(*grid);
//   openvdb::VectorGrid::Ptr resultGrid = openvdb::VectorGrid::create();
//   resultGrid = myGrid.process();
//   
//   openvdb::math::Vec3<float> outgra = resultGrid->background();
//   std::cout << "background value is : " << outgra << std::endl;
//   
//   std::cout << "gradient values are : " << std::endl;
//   for (openvdb::VectorGrid::ValueOnIter iter = resultGrid->beginValueOn(); iter; ++iter) {
//     std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
//   }
  
}
