///////////////////////////////////////////////////////////////////////////
//
/// @file SOP_OpenVDB_Vector_Angle.cc
///
/// @author maphysart
///
/// @brief vector angle node

////////////////////////////////////////

#include <houdini_utils/ParmFactory.h>
#include <openvdb_houdini/Utils.h>
#include <openvdb_houdini/SOP_NodeVDB.h>
#include <openvdb_houdini/GEO_PrimVDB.h>
#include <openvdb_houdini/GU_PrimVDB.h>
#include <openvdb_houdini/UT_VDBTools.h>

#include <openvdb/tools/ValueTransformer.h>

#include <UT/UT_Interrupt.h>
#include <GEO/GEO_Point.h>

#include <limits>   		// std::numeric_limits
#include <string>		// std::string::compare

#define PI 3.1415926535897931

namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;

class SOP_OpenVDB_Vector_Angle : public hvdb::SOP_NodeVDB
{
public:	
	
	SOP_OpenVDB_Vector_Angle(OP_Network *net, const char *name, OP_Operator *op);
    
	virtual ~SOP_OpenVDB_Vector_Angle(){};
    
	static OP_Node* factory(OP_Network*, const char* name, OP_Operator*);
	
protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);
};


////////////////////////////////////////

OP_Node* 
SOP_OpenVDB_Vector_Angle::factory(OP_Network* net,
    const char* name, OP_Operator* op)
{
    return new SOP_OpenVDB_Vector_Angle(net, name, op);
}

SOP_OpenVDB_Vector_Angle::SOP_OpenVDB_Vector_Angle(OP_Network* net,
    const char* name, OP_Operator* op) : hvdb::SOP_NodeVDB(net, name, op)
{
}


// Build UI and register this operator.
void
newSopOperator(OP_OperatorTable* table)
{
    if (table == NULL) return;

    hutil::ParmList parms;
	
	parms.add(hutil::ParmFactory(PRM_STRING, "Group",  "Group")
		.setDefault(0, ""));
		
	parms.add(hutil::ParmFactory(PRM_STRING, "groupA",  "Group A")
		.setChoiceList(&hutil::PrimGroupMenu));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "groupB",  "Group B")
		.setChoiceList(&hutil::PrimGroupMenu));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "gradName", "Gradient Grid Name")
		.setDefault(0, "gradGrid"));
	
    hvdb::OpenVDBOpFactory("OpenVDB Vector Angle", SOP_OpenVDB_Vector_Angle::factory, parms, *table)
		.addInput("VDB")
		.addInput("VDBA")
		.addInput("VDBB");
}


OP_ERROR
SOP_OpenVDB_Vector_Angle::cookMySop(OP_Context &context)
{
	try {
		hutil::ScopedInputLock lock(*this, context);
		//duplicatePointSource(0, context);
		const fpreal time = context.getTime();
		duplicateSource(0, context);
		
		/*
		 * Get params
		 */	
		const GU_Detail* vdbA = inputGeo(1, context);
		const GU_Detail* vdbB = inputGeo(2, context);

		UT_String gradNameStr;
		evalString(gradNameStr, "gradName", 0, time);
		
		UT_String groupStr;
		evalString(groupStr, "group", 0, time);
	    const GA_PrimitiveGroup *group = matchGroup(*gdp, groupStr.toStdString());
		
		UT_String groupAStr;	
		evalString(groupAStr, "groupA", 0, time);
		const GA_PrimitiveGroup *groupA = matchGroup(const_cast<GU_Detail&>(*vdbA), groupAStr.toStdString());
		
		UT_String groupBStr;
		evalString(groupBStr, "groupB", 0, time);
		const GA_PrimitiveGroup *groupB = matchGroup(const_cast<GU_Detail&>(*vdbB), groupBStr.toStdString());
		
		openvdb::FloatGrid *  grid;
		for (hvdb::VdbPrimIterator it(gdp, group); it; ++it) {
			GU_PrimVDB* vdb = *it;
			if (vdb->getGrid().isType<openvdb::FloatGrid>())
				grid = static_cast<openvdb::FloatGrid *>(&vdb->getGrid());
			else 
				gdp->destroyPrimitive(*vdb, true);
		}	
				
				
		//std::cout << "Required gradient name : " << gradNameStr.toStdString() << std::endl;
		const hvdb::GU_PrimVDB *pgradA = NULL;
		for (hvdb::VdbPrimCIterator it(vdbA, groupA); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
			if ( gridName.compare( gradNameStr.toStdString() ) == 0 )
				pgradA = *it;
		}
		
				//std::cout << "Required gradient name : " << gradNameStr.toStdString() << std::endl;
		const hvdb::GU_PrimVDB *pgradB = NULL;
		for (hvdb::VdbPrimCIterator it(vdbB, groupB); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
			if ( gridName.compare( gradNameStr.toStdString() ) == 0 )
				pgradB = *it;
		}
		
		const openvdb::VectorGrid *  gridA;
		if ( pgradA != NULL )
			gridA = static_cast<const openvdb::VectorGrid *>(&pgradA->getGrid());
		else 
			throw std::runtime_error("Cannot find gradient grid.");
		
		const openvdb::VectorGrid *  gridB;
		if ( pgradB != NULL )
			gridB = static_cast<const openvdb::VectorGrid *>(&pgradB->getGrid());
		else 
			throw std::runtime_error("Cannot find gradient grid.");
		
		
		openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
		
		const openvdb::math::Transform& gridAXform = gridA->constTransform();
		const openvdb::math::Transform& gridBXform = gridB->constTransform();
		
		for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter.test(); ++iter) {
			
			hvdb::Interrupter progress("Calculate vector angle based on two gradient vector grid");
					
			if (iter.isVoxelValue()) { // set a single voxel
				openvdb::Vec3f pos = grid->constTransform().indexToWorld(iter.getCoord());
				openvdb::Vec3f gradA, gradB;
				openvdb::tools::BoxSampler::sample(gridA->tree(), gridAXform.worldToIndex(pos), gradA);
				openvdb::tools::BoxSampler::sample(gridB->tree(), gridBXform.worldToIndex(pos), gradB);
				float temp = openvdb::math::angle(gradA, gradB) / PI * 180;
// 				float temp = acos( gradA.dot(gradB) ) / PI * 180;
				if ( gradA == openvdb::Vec3f(0.0f, 0.0f, 0.0f) || gradB == openvdb::Vec3f(0.0f, 0.0f, 0.0f) )
					temp = 0.0f;
				iter.setValue(temp);
			} else { // fill an entire tile
				openvdb::CoordBBox bbox;
				iter.getBoundingBox(bbox);
				accessor.getTree()->fill(bbox, 0.0f);
			}
		}
		
		grid->setName("angleGrid");
		
// 		std::cout << "gradient dot angle values are : " << std::endl;
// 		for (openvdb::FloatGrid::ValueOnCIter iter = outGrid->cbeginValueOn(); iter; ++iter) {
// 			float value = iter.getValue();
// // 			if ( value.x() > 0 ) 
// 			std::cout << "Grid world" << outGrid->constTransform().indexToWorld(iter.getCoord()) << " index" << iter.getCoord() << " = " << value << std::endl;
// 		}
		
	} catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
    }
    
    return error();
	
}
