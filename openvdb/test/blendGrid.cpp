#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/GridOperators.h>

#include <algorithm>    // std::max
#include <math.h>       /* tan */		
#include <iostream>

using namespace openvdb;

float distance(float x, float y, float ox, float oy){
	return sqrt( (x - ox) * (x - ox) + (y - oy) * (y - oy) );
}


float calculateGValue(float x, float y, float tan_theta, float voxelSize){
	float r1, r2, k;
	
	// r1
	if ( x > 0.5f )
		k = 0.5f * ( 7 - 5 * tan_theta ) * x + 1.75f * (tan_theta - 1);
	else 
		k = tan_theta * x;
	
	r1 = distance(x, y, k, k) - x + k;
	if(r1 < 0)
		std::cout << "r1 < 0" << std::endl;
	if(r1 == 0)
		return x;
	
	// r2
	bool negative = false;
	float ci1 = x;
	do {
		ci1 += voxelSize;
		if ( ci1 > 0.5f )
			k = 0.5f * ( 7 - 5 * tan_theta ) * ci1 + 1.75f * (tan_theta - 1);
		else 
			k = tan_theta * ci1;
		
		r2 = ci1 - k - distance(x, y, k, k);
		if(r2 == 0)
			return x;
		if(r2 < 0)
			negative = true;
	}while( r2 < 0 && ci1 < 1.0f);
	if ( ci1 >= 1.0f && x > 0.65f ){
// 		std::cout << "error! C(i+1) > 1 then take Ci=" << x << std::endl;
		return x;
	}
// 	if (negative)
// 		std::cout << "r1= " << r1 << " r2=" << r2 << " when x=" << x << " y=" << y << " tan_theta=" << tan_theta << std::endl;
	return ( r2 * x + r1 * ci1 ) / (r1 + r2);
}

	
void makeBox(FloatGrid::Ptr grid, const CoordBBox& indexBB, float voxelSize)
{
	typename FloatGrid::Accessor accessor = grid->getAccessor();

	float x, y, tan_theta, value;
	for (Int32 k = indexBB.min().z(); k <= indexBB.max().z(); ++k) { // tan theta
		float tan_theta = k * voxelSize;
		for (Int32 i = indexBB.min().x(); i <= indexBB.max().x(); ++i) {
			float x = i * voxelSize; // Ci = x;
			Int32 j = indexBB.min().y();
			float y = j * voxelSize;
			while ( y <= x ) {
				if ( x > 0.7f ){
					value = x;
				} else if( x <= 0.5f ){
					if(y < tan_theta * x)
						value = x;
					else {
						// calculate G value
						value = calculateGValue(x, y, tan_theta, voxelSize);
					}
				} else{
					if ( y < ( 0.5f * ( 7 - 5 * tan_theta ) * x + 1.75f * (tan_theta - 1) ) )
						value = x;
					else {
						// calculate G value
						value = calculateGValue(x, y, tan_theta, voxelSize);
					}
				}
// 				if ( value < 0.0f )
// 					std::cout << "error! value=" << value << " < 0! " << std::endl;
				accessor.setValue(Coord(i, j, k), value);
				accessor.setValue(Coord(j, i, k), value);
				j++;
				y += voxelSize;
			}
		}
	}
	grid->setTransform(openvdb::math::Transform::createLinearTransform(voxelSize));
}


int main()
{
	openvdb::initialize();

	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);

	float voxelSize = 1.0f/128;
	float bboxSize = 1.0 / voxelSize;
	CoordBBox indexBB(Coord(0, 0, 0), Coord(bboxSize, bboxSize, bboxSize));
	makeBox(grid, indexBB, voxelSize);

	grid->setName("blendGrid");
	
// 	openvdb::math::Transform& gridform = grid->transform();
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
	for ( int i = 0; i < bboxSize; i++){
		for ( int j = 0; j < bboxSize; j++){
			for ( int k = 0; k < bboxSize; k++){
				openvdb::Coord xyz(i, j, k);
				std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
			}
		}
	}
	
    openvdb::io::File file("blendGrid.vdb");

    openvdb::GridPtrVec grids;
	grids.push_back(grid);

    file.write(grids);
    file.close();
  
}
