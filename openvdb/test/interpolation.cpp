#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <iostream>


using namespace openvdb;


void makeBox(FloatGrid::Ptr grid, const CoordBBox& indexBB, double step)
{
  typename FloatGrid::Accessor accessor = grid->getAccessor();

  for (Int32 i = indexBB.min().x(); i <= indexBB.max().x(); ++i) {
    for (Int32 j = indexBB.min().y(); j <= indexBB.max().y(); ++j) {
      for (Int32 k = indexBB.min().z(); k <= indexBB.max().z(); ++k) {

        float value = i + j + k;

        accessor.setValue(Coord(i, j, k), value);
      }
    }
  }	

  grid->setTransform(openvdb::math::Transform::createLinearTransform(step));
}

void makeBoxVec3(VectorGrid::Ptr grid, const CoordBBox& indexBB, double step)
{
  typename VectorGrid::Accessor accessor = grid->getAccessor();

  for (Int32 i = indexBB.min().x(); i <= indexBB.max().x(); ++i) {
    for (Int32 j = indexBB.min().y(); j <= indexBB.max().y(); ++j) {
      for (Int32 k = indexBB.min().z(); k <= indexBB.max().z(); ++k) {

        Vec3f value = Vec3f(i, j, k);

        accessor.setValue(Coord(i, j, k), value);
      }
    }
  }	

  grid->setTransform(openvdb::math::Transform::createLinearTransform(step));
}

int main()
{
	openvdb::initialize();

	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

	openvdb::VectorGrid::Ptr gridVec = openvdb::VectorGrid::create();
	openvdb::VectorGrid::Accessor accessorVec = gridVec->getAccessor();
	
	float size = 1, step = 1;
	CoordBBox indexBB(Coord(-size, -size, -size), Coord(size, size, size));
	makeBox(grid, indexBB, step);
	makeBoxVec3(gridVec, indexBB, step);
	
	
	const openvdb::math::Transform &	xform = grid->transform();
	
	
	float i,j,k;
	float value;
	Vec3f vec;
	std::cout << "Please enter coordinate value: ";
	while(std::cin >> i >> j >> k){
		
	    openvdb::Vec3f vpos = openvdb::Vec3f(i, j, k);
		
		vpos = xform.worldToIndex(vpos);

		openvdb::tools::BoxSampler::sample(grid->tree(), vpos, value);
		std::cout << "Grid[" << i << ", " <<  j << ", " << k << "] " << " = " << value << std::endl;
		
		openvdb::tools::BoxSampler::sample(gridVec->tree(), vpos, vec);
		std::cout << "Grid[" << i << ", " <<  j << ", " << k << "] " << " = " << vec << std::endl;
		
		std::cout << "Please enter coordinate value: ";
		
	}
	

}
