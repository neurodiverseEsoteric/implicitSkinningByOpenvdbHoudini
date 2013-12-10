///////////////////////////////////////////////////////////////////////////
//
/// @file SOP_OpenVDB_Implicit_Surf.cc
///
/// @author maphysart
///
/// @brief Implicit Surface node

////////////////////////////////////////


#include <houdini_utils/ParmFactory.h>
#include <openvdb_houdini/Utils.h>
#include <openvdb_houdini/SOP_NodeVDB.h>
#include <openvdb_houdini/GEO_PrimVDB.h>
#include <openvdb_houdini/GU_PrimVDB.h>
#include <UT/UT_Interrupt.h>

#include <openvdb/tools/GridOperators.h>

#include <iostream>     // std::cout
#include <sstream>
#include <limits>       // std::numeric_limits

#include "hrbf_core.h"
#include "hrbf_phi_funcs.h"

namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;

class SOP_OpenVDB_Implicit_Surf : public hvdb::SOP_NodeVDB
{
public:
	SOP_OpenVDB_Implicit_Surf(OP_Network *net, const char *name, OP_Operator *op);
    
	virtual ~SOP_OpenVDB_Implicit_Surf(){};
	
    virtual void getDescriptiveParmName(UT_String& s) const { s = "file_name"; }
    
	static OP_Node* factory(OP_Network*, const char* name, OP_Operator*);
	
	static int writeNowCallback(void* data, int index, float now, const PRM_Template*);
	
protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);
	void fillVDB(openvdb::FloatGrid::Ptr grid, const openvdb::CoordBBox& bbox, const HRBF_fit<float, 3, Rbf_pow3<float> >& hrbf, float voxelSize, float radius);
	float remap(float value, float r);
private:

};


////////////////////////////////////////


OP_Node* 
SOP_OpenVDB_Implicit_Surf::factory(OP_Network* net,
    const char* name, OP_Operator* op)
{
    return new SOP_OpenVDB_Implicit_Surf(net, name, op);
}

SOP_OpenVDB_Implicit_Surf::SOP_OpenVDB_Implicit_Surf(OP_Network* net,
    const char* name, OP_Operator* op) : hvdb::SOP_NodeVDB(net, name, op)
{
}


int
SOP_OpenVDB_Implicit_Surf::writeNowCallback(
    void *data, int /*index*/, float /*now*/,
    const PRM_Template*)
{
    if (SOP_OpenVDB_Implicit_Surf* self = static_cast<SOP_OpenVDB_Implicit_Surf*>(data)) {
        self->forceRecook();
        return 1;
    }
    return 0;
}


float
SOP_OpenVDB_Implicit_Surf::remap(float value, float r)
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


// Build UI and register this operator.
void
newSopOperator(OP_OperatorTable* table)
{
    if (table == NULL) return;

    hutil::ParmList parms;

    // Output VDB File name
    parms.add(hutil::ParmFactory(PRM_FILE, "vdbFileName", "Output VDB File")
        .setDefault(0, "./output.vdb")
        .setHelpText("Path name for the output VDB file"));
	
    parms.add(hutil::ParmFactory(PRM_STRING, "vdbName", "vdb name")
        .setHelpText("vdb name."));
	
    parms.add(hutil::ParmFactory(PRM_STRING, "vdbGradName", "vdb gradient name")
        .setHelpText("vdb gradient name."));

	// Voxel Size
	parms.add(hutil::ParmFactory(PRM_FLT_J, "voxelSize", "Voxel Size")
        .setDefault(PRMpointOneDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 2));
	
	// Affect Radius
	parms.add(hutil::ParmFactory(PRM_INT, "radius", "Affect Radius")
        .setDefault(PRMfiveDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 20));

    // "Write Now" button
    parms.add(hutil::ParmFactory(PRM_CALLBACK, "write", "Write Now")
        .setCallbackFunc(&SOP_OpenVDB_Implicit_Surf::writeNowCallback)
        .setHelpText("Click to write the output file."));
	
    // Register this operator.
    hvdb::OpenVDBOpFactory("OpenVDB Implicit Surf", SOP_OpenVDB_Implicit_Surf::factory, parms, *table)
        .addInput("Input group of points");
	
}


void
SOP_OpenVDB_Implicit_Surf::fillVDB(openvdb::FloatGrid::Ptr grid, const openvdb::CoordBBox& bbox, const HRBF_fit<float, 3, Rbf_pow3<float> >& hrbf, float voxelSize, float radius)
{
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
	
// 	std::cout << "Fill vdb :" << std::endl;
	
	openvdb::Int32 maxX = static_cast<openvdb::Int32>(bbox.max().x() / voxelSize), minX = static_cast<openvdb::Int32>(bbox.min().x() / voxelSize);
	openvdb::Int32 maxY = static_cast<openvdb::Int32>(bbox.max().y() / voxelSize), minY = static_cast<openvdb::Int32>(bbox.min().y() / voxelSize);
	openvdb::Int32 maxZ = static_cast<openvdb::Int32>(bbox.max().z() / voxelSize), minZ = static_cast<openvdb::Int32>(bbox.min().z() / voxelSize);
	
	float r = voxelSize * radius, value;
	for (openvdb::Int32 i = minX; i <= maxX; ++i) {
		for (openvdb::Int32 j = minY; j <= maxY; ++j) {
			for (openvdb::Int32 k = minZ; k <= maxZ; ++k) {
				// transform point (i, j, k) of index space into world space
				Eigen::Matrix<float, 3, 1> pos;
				pos(0) = i * voxelSize; pos(1) = j * voxelSize; pos(2) = k * voxelSize;
				
// 				std::cout << "pos: " << i << " " << j << " " << k << " " << std::endl;
// 				std::cout << pos << std::endl;
				// compute level set function value
				value = hrbf.eval(pos);
				value = remap(value, r);
// 				std::cout << "value: " << value << std::endl;
				accessor.setValue(openvdb::Coord(i, j, k), value);
			}
		}
	}	
	
	grid->setTransform(openvdb::math::Transform::createLinearTransform(voxelSize));
}

OP_ERROR
SOP_OpenVDB_Implicit_Surf::cookMySop(OP_Context &context)
{
	try {
		hutil::ScopedInputLock lock(*this, context);
		hvdb::Interrupter boss("Converting geometry to volume");
		gdp->clearAndDestroy();

        // Get params
		const fpreal time = context.getTime();
		float voxelSize	= evalFloat("voxelSize", 0, time);
		int	 radius 	= evalInt("radius", 0, time);
	
		UT_String fileNameStr;
		evalString(fileNameStr, "vdbFileName", 0, time);
		const std::string filename = fileNameStr.toStdString();
		if (filename.empty()) {
			addWarning(SOP_MESSAGE, "no name given for the output file");
			return error();
		}

		UT_String vdbNameStr;
		evalString(vdbNameStr, "vdbName", 0, time);

		UT_String vdbGradNameStr;
		evalString(vdbGradNameStr, "vdbGradName", 0, time);
		
		/*
		 * Make the grid and its gradient grid
		 */
		// calculate the hrbf func
		//typedef Eigen::Matrix<float, 3, 1>	Vector;
		std::vector< Eigen::Matrix<float, 3, 1> > points;
		std::vector< Eigen::Matrix<float, 3, 1> > normals;
		
		const GU_Detail* gdp = inputGeo(0, context);
		
		// get the handle of points' pos and normals
		GA_ROHandleV3 pHandle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
		GA_ROHandleV3 nHandle(gdp->findAttribute(GA_ATTRIB_POINT, "N"));
		
		GA_Offset ptoff;
		Eigen::Matrix<float, 3, 1> tempPos, tempN;
		float maxX = std::numeric_limits<float>::min(), minX = std::numeric_limits<float>::max(); 
		float maxY = std::numeric_limits<float>::min(), minY = std::numeric_limits<float>::max(); 
		float maxZ = std::numeric_limits<float>::min(), minZ = std::numeric_limits<float>::max();
		GA_FOR_ALL_PTOFF(gdp, ptoff)
		{
			UT_Vector3 pValue = pHandle.get(ptoff);
			UT_Vector3 nValue = nHandle.get(ptoff);
			tempPos(0) = pValue.x();	tempPos(1) = pValue.y();	tempPos(2) = pValue.z();
			
			if ( minX > pValue.x() ) minX = pValue.x();
			if ( maxX < pValue.x() ) maxX = pValue.x();
			if ( minY > pValue.y() ) minY = pValue.y();
			if ( maxY < pValue.y() ) maxY = pValue.y();
			if ( minZ > pValue.z() ) minZ = pValue.z();
			if ( maxZ < pValue.z() ) maxZ = pValue.z();
			
			tempN(0) = nValue.x();	tempN(1) = nValue.y();	tempN(2) = nValue.z();
			points.push_back(tempPos);
			normals.push_back(tempN);
		}
		
		// get the hrbf alphas and betas
		HRBF_fit< float, 3, Rbf_pow3<float> > hrbf_surf;
		hrbf_surf.hermite_fit(points, normals);
		
// 		std::cout << "fit surface " << std::endl;
// 		std::cout << "alphas: " << std::endl;
// 		std::cout << hrbf_surf._alphas << std::endl;
// 		std::cout << "betas: " << std::endl;
// 		std::cout << hrbf_surf._betas << std::endl;
		
		// make the vdb grid and gradient grid by hrbf
		openvdb::GridPtrSet outGrids;
		
		openvdb::FloatGrid::Ptr outGrid = openvdb::FloatGrid::create(0.0f);
		openvdb::CoordBBox bbox(openvdb::Coord(minX - voxelSize * radius, minY - voxelSize * radius, minZ - voxelSize * radius), \
			openvdb::Coord(maxX + voxelSize * radius, maxY + voxelSize * radius, maxZ + voxelSize * radius));	
		
// 		std::cout << "bbox: " << std::endl;
// 		std::cout << coordAsString(bbox.min(),",") << std::endl;
// 		std::cout << coordAsString(bbox.max(),",") << std::endl;
		
		fillVDB(outGrid, bbox, hrbf_surf, voxelSize, radius);
		outGrid->setName(vdbNameStr.toStdString());
		//outGrid->tree().prune();
		outGrid->signedFloodFill();
		
		openvdb::VectorGrid::Ptr outGradGrid = openvdb::tools::gradient(*outGrid);
		outGradGrid->setName(vdbGradNameStr.toStdString());
// 		outGradGrid->tree().prune();
		
		outGrids.insert(outGrid);
		outGrids.insert(outGradGrid);
		
		// Create a VDB file object.
		openvdb::io::File file(filename);
		file.write(outGrids);
		file.close();
		
	} catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
    }
    
    return error();
	
}
