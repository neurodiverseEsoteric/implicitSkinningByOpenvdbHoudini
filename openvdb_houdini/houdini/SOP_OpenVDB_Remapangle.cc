///////////////////////////////////////////////////////////////////////////
//
/// @file SOP_OpenVDB_Remapangle.cc
///
/// @author maphysart
///
/// @brief remap angle node

////////////////////////////////////////

#include <houdini_utils/ParmFactory.h>
#include <openvdb_houdini/Utils.h>
#include <openvdb_houdini/SOP_NodeVDB.h>
#include <openvdb_houdini/GEO_PrimVDB.h>
#include <openvdb_houdini/GU_PrimVDB.h>
#include <openvdb_houdini/UT_VDBTools.h>
#include <openvdb_houdini/UT_VDBUtils.h>

#include <UT/UT_Interrupt.h>
#include <GEO/GEO_Point.h>

#include <limits>   		// std::numeric_limits
#include <string>		// std::string::compare

#include <math.h>       /* exp */

#define PI 3.1415926535897931
#define PI4 0.785398163
#define PI180 0.017453293

namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;

class SOP_OpenVDB_Remapangle : public hvdb::SOP_NodeVDB
{
public:	
	SOP_OpenVDB_Remapangle(OP_Network *net, const char *name, OP_Operator *op);
    
	virtual ~SOP_OpenVDB_Remapangle(){};
	
    //virtual void getDescriptiveParmName(UT_String& s) const { s = "file_name"; }
    
	static OP_Node* factory(OP_Network*, const char* name, OP_Operator*);
	
protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);

private:
	double kaba(double ratio){ return 1.0 - exp(1.0 	- 1.0 / ( 1 - exp(1 - 1.0/ratio))); };
};

////////////////////////////////////////


OP_Node* 
SOP_OpenVDB_Remapangle::factory(OP_Network* net,
    const char* name, OP_Operator* op)
{
    return new SOP_OpenVDB_Remapangle(net, name, op);
}

SOP_OpenVDB_Remapangle::SOP_OpenVDB_Remapangle(OP_Network* net,
    const char* name, OP_Operator* op) : hvdb::SOP_NodeVDB(net, name, op)
{
}


// Build UI and register this operator.
void
newSopOperator(OP_OperatorTable* table)
{
    if (table == NULL) return;

    hutil::ParmList parms;
	
	parms.add(hutil::ParmFactory(PRM_STRING, "group",  "Group")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));

	parms.add(hutil::ParmFactory(PRM_FLT_J, "alpha1", "Start Angle")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, PI));

	parms.add(hutil::ParmFactory(PRM_FLT_J, "alpha2", "Middle Angle")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, PI));

	parms.add(hutil::ParmFactory(PRM_FLT_J, "alpha3", "End Angle")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, PI));
	
	parms.add(hutil::ParmFactory(PRM_FLT_J, "theta1", "Start theta Angle")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, PI4));
	
	parms.add(hutil::ParmFactory(PRM_FLT_J, "theta2", "Middle theta Angle")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, PI4));
	
	parms.add(hutil::ParmFactory(PRM_FLT_J, "theta3", "End theta Angle")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, PI4));

	parms.add(hutil::ParmFactory(PRM_FLT_J, "w1", "Steep1 ")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0.01, PRM_RANGE_UI, 10));

	parms.add(hutil::ParmFactory(PRM_FLT_J, "w2", "Steep2")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0.01, PRM_RANGE_UI, 10));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "outGroup",  "output Group")
		.setDefault(0, "thetaGrid")
		.setHelpText("Output theta grid name."));
	
	// Register this operator.
    hvdb::OpenVDBOpFactory("OpenVDB Remapangle", SOP_OpenVDB_Remapangle::factory, parms, *table)
		.addInput("VDBs");
}


OP_ERROR
SOP_OpenVDB_Remapangle::cookMySop(OP_Context& context)
{
    try {
        hutil::ScopedInputLock lock(*this, context);
        const fpreal time = context.getTime();
		gdp->clearAndDestroy();

		const GU_Detail* vdbA = inputGeo(0, context);
		
		float alpha1 = evalFloat("alpha1", 0, time);
		float alpha2 = evalFloat("alpha2", 0, time);
		float alpha3 = evalFloat("alpha3", 0, time);
		float theta1 = evalFloat("theta1", 0, time);
		float theta2 = evalFloat("theta2", 0, time);
		float theta3 = evalFloat("theta3", 0, time);
		float w1 = evalFloat("w1", 0, time);		
		float w2 = evalFloat("w2", 0, time);		
		
		UT_String groupStr;
		evalString(groupStr, "group", 0, time);
	    const GA_PrimitiveGroup *group = matchGroup(const_cast<GU_Detail&>(*vdbA), groupStr.toStdString());
		
		
		// Create a output grid and group
        GA_PrimitiveGroup* outGroup = NULL;
        UT_String outGroupStr;
        evalString(outGroupStr, "outGroup", 0, time);
        if(outGroupStr.isstring()) {
            outGroup = gdp->newPrimitiveGroup(outGroupStr.buffer());
        }
		openvdb::FloatGrid::Ptr outGrid = openvdb::FloatGrid::create();
		openvdb::FloatGrid::Accessor accessor = outGrid->getAccessor();	
		
		
		hvdb::Interrupter progress("Remap angle according to remap function.");
		
		for (hvdb::VdbPrimCIterator it(vdbA, group); it; ++it) {
			if (progress.wasInterrupted()) 
				throw std::runtime_error("was interrupted");	
			
			const GU_PrimVDB* vdb = *it;
			if (vdb->getGrid().isType<openvdb::FloatGrid>()){
				const openvdb::FloatGrid *  grid = static_cast<const openvdb::FloatGrid *>(&vdb->getGrid());
// 				int success = UTvdbProcessTypedGrid(UTvdbGetGridType(vdb->getGrid()), &vdb->getGrid(), Op);
// 				int j = 0;
				outGrid->setTransform(grid->transform().copy());
				for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter.test(); ++iter){
					const float value = iter.getValue() * PI180;
					float temp;
					
					// remap angle to theta 
					if (value <= alpha1)
						temp = theta1;
					else if (value >= alpha3)
						temp = theta3;
					else {
						if ( value < alpha2 ) {
							temp = pow(kaba( (value - alpha2) / (alpha1 - alpha2 ) ), w1) * ( theta1 - theta2 ) + theta2;
						}
						else {
							temp = pow(kaba( (value - alpha2) / (alpha3 - alpha2 ) ), w2) * ( theta3 - theta2 ) + theta2;
						}	
					}
// 					iter.setValue(temp);
// 					// set theta grid value
					if (iter.isVoxelValue()) { // set a single voxel
						accessor.setValue(iter.getCoord(), temp);
					} else { // fill an entire tile
						openvdb::CoordBBox bbox;
						iter.getBoundingBox(bbox);
						accessor.getTree()->fill(bbox, temp);
					}
				}
				
// 				std::cout << "remap grid value are : " << std::endl;
// 				for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter.test(); ++iter){
// 					float value = iter.getValue();
// 		// 			if ( value.x() > 0 ) 
// 					std::cout << "Grid world" << grid->constTransform().indexToWorld(iter.getCoord()) << " index" << iter.getCoord() << " = " << value << std::endl;
// 				}
			}
			break;
		}
		
		GEO_PrimVDB* vda = hvdb::createVdbPrimitive(*gdp, outGrid, outGroupStr.toStdString().c_str());		
		if (outGroup) outGroup->add(vda);
		
	} catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
    }
    
    return error();
}