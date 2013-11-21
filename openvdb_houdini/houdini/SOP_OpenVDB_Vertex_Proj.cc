///////////////////////////////////////////////////////////////////////////
//
/// @file SOP_OpenVDB_Vertex_Proj.cc
///
/// @author maphysart
///
/// @brief vertex proj node

////////////////////////////////////////


#include <houdini_utils/ParmFactory.h>
#include <openvdb_houdini/Utils.h>
#include <openvdb_houdini/SOP_NodeVDB.h>
#include <openvdb_houdini/GEO_PrimVDB.h>
#include <openvdb_houdini/GU_PrimVDB.h>
#include <openvdb_houdini/UT_VDBTools.h>

#include <UT/UT_Interrupt.h>
#include <GEO/GEO_Point.h>

#include <limits>   		// std::numeric_limits
#include <string>		// std::string::compare

namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;

////////////////////////////////////////


class SOP_OpenVDB_Vertex_Proj : public hvdb::SOP_NodeVDB
{
public:
	SOP_OpenVDB_Vertex_Proj(OP_Network *net, const char *name, OP_Operator *op);
    
	virtual ~SOP_OpenVDB_Vertex_Proj(){};
	
    //virtual void getDescriptiveParmName(UT_String& s) const { s = "file_name"; }
    
	static OP_Node* factory(OP_Network*, const char* name, OP_Operator*);
	
protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);
	void relaxVertex(openvdb::FloatGrid::Ptr grid);

private:
	fpreal 	mError;
	uint		mIterations;
	fpreal 	mStopAngle;
    static const float omega;
};


////////////////////////////////////////

const float SOP_OpenVDB_Vertex_Proj::omega = 0.35;


OP_Node* 
SOP_OpenVDB_Vertex_Proj::factory(OP_Network* net,
    const char* name, OP_Operator* op)
{
    return new SOP_OpenVDB_Vertex_Proj(net, name, op);
}

SOP_OpenVDB_Vertex_Proj::SOP_OpenVDB_Vertex_Proj(OP_Network* net,
    const char* name, OP_Operator* op) : hvdb::SOP_NodeVDB(net, name, op), mError(.0), mIterations(5), mStopAngle(20)
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
			
	parms.add(hutil::ParmFactory(PRM_STRING, "gridName", "Value Grid Name")
		.setDefault(0, "valueGrid")
        .setHelpText("value grid name."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "gradName", "Gradient Grid Name")
		.setDefault(0, "gradGrid")
        .setHelpText("gradient grid name."));
	
	parms.add(hutil::ParmFactory(PRM_FLT_J, "error", "Error Bound")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 0.2));
	
	parms.add(hutil::ParmFactory(PRM_INT, "iter", "Iterations")
        .setDefault(PRMfiveDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 10));

	parms.add(hutil::ParmFactory(PRM_FLT_J, "angle", "Stop Angle")
        .setDefault(PRM20Defaults)
        .setRange(PRM_RANGE_RESTRICTED, 10, PRM_RANGE_UI, 100));
	
    // Register this operator.
    hvdb::OpenVDBOpFactory("OpenVDB Vertex Proj", SOP_OpenVDB_Vertex_Proj::factory, parms, *table)
		.addInput("Polygon Mesh")
		.addInput("VDBs");
}


void
SOP_OpenVDB_Vertex_Proj::relaxVertex(openvdb::FloatGrid::Ptr grid)
{

}

OP_ERROR
SOP_OpenVDB_Vertex_Proj::cookMySop(OP_Context &context)
{
	try {
        hutil::ScopedInputLock lock(*this, context);
        duplicatePointSource(0, context);
        const fpreal time = context.getTime();
		
		
		/*
		 * Get params
		 */
		const GU_Detail* vdbGdp = inputGeo(1, context);
        
		mError = evalFloat("error", 0, time);
		mIterations = evalInt("iter", 0, time);
		mStopAngle = evalFloat("angle", 0, time);
		
		UT_String gridNameStr;
		evalString(gridNameStr, "gridName", 0, time);
		UT_String gradNameStr;
		evalString(gradNameStr, "gradName", 0, time);
		
		UT_String groupStr;
        evalString(groupStr, "group", 0, time);
		const GA_PrimitiveGroup *group = matchGroup(const_cast<GU_Detail&>(*vdbGdp), groupStr.toStdString());
		
		std::cout << "Required grid name : " << gridNameStr.toStdString() << std::endl;
		std::cout << "Required gradient name : " << gradNameStr.toStdString() << std::endl;
		const hvdb::GU_PrimVDB *pgrid = NULL, *pgrad = NULL;
		for (hvdb::VdbPrimCIterator it(vdbGdp, group); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
			std::cout << "Grid name : " << gridName << std::endl;
			if ( gridName.compare( gridNameStr.toStdString() ) == 0 )
				pgrid = *it;
			if ( gridName.compare( gradNameStr.toStdString() ) == 0 )
				pgrad = *it;
		}
		
		openvdb::FloatGrid::Ptr grid;
		openvdb::VectorGrid::Ptr grad;
		if ( pgrid != NULL && pgrad != NULL ) {
			grid = openvdb::gridPtrCast<openvdb::FloatGrid>( pgrid->getGrid().deepCopyGrid() );
			grad = openvdb::gridPtrCast<openvdb::VectorGrid>( pgrad->getGrid().deepCopyGrid() );	
		}
		else{
			throw std::runtime_error("Cannot find value or gradient grid.");
		}
		// gradient and grid have same transform
		openvdb::math::Transform& gridform = grid->transform();

		
		/*
		 * Project 
		 */
		hvdb::Interrupter progress("Projecting vertex based on vdb grid");
		fpreal error = std::numeric_limits<double>::max();
		
		//pHandle(gdp->findAttribute(GA_ATTRIB_POINT, "nearList"));
		GA_ROHandleR restVHandle(gdp->findAttribute(GA_ATTRIB_POINT, "restValue"));
		GA_RWHandleV3 pHandle(gdp->getP());
		
		for ( uint i = 0; i < mIterations && error > mError; i++ ){
				
			if (progress.wasInterrupted()) {
                throw std::runtime_error("Projection was interrupted");
			}
            
            	float value;
			openvdb::Vec3f vec, delta;
			GA_Offset ptoff;
			
			int total = 0;
			for (GA_Iterator it(gdp->getPointRange()); !it.atEnd(); it.advance())
			{
				total++;
				ptoff = it.getOffset();
				UT_Vector3 pos = pHandle.get(ptoff);
				fpreal restv = restVHandle.get(ptoff);
				
				openvdb::Vec3f vpos = openvdb::Vec3f(pos.x(), pos.y(), pos.z());
				vpos = gridform.worldToIndex(vpos);
				
				openvdb::tools::BoxSampler::sample(grid->tree(), vpos, value);
				openvdb::tools::BoxSampler::sample(grad->tree(), vpos, vec);
			
				// find the point attr neighboutList and valueRest
				delta = vec * ( omega * (value - restv) / vec.lengthSqr() );
				pos.x() += delta.x();
				pos.y() += delta.y();
				pos.z() += delta.z();
				
				pHandle.set(ptoff, pos);
			}
			std::cout << "Total points : " << total << std::endl;
            
		}	
		
	} catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
    }
    
    return error();
	
}
