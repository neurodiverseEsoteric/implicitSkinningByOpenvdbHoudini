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

#include <algorithm>		// std::max
#include <math.h>    	// fabs 
#include <boost/dynamic_bitset.hpp>

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
	fpreal 	mBias;
	uint 	mRelaxFreq;
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
    const char* name, OP_Operator* op) : hvdb::SOP_NodeVDB(net, name, op), mError(.0), mIterations(5), mStopAngle(55), mBias(1.0f), mRelaxFreq(3)
{
}

// Build UI and register this operator.
void
newSopOperator(OP_OperatorTable* table)
{
    if (table == NULL) return;

    hutil::ParmList parms;
	
	// TODO only deal with some group of points
// 	parms.add(hutil::ParmFactory(PRM_STRING, "groupPts",  "Group Points")
// 		.setChoiceList(&hutil::PrimGroupMenu)
// 		.setHelpText("Choose only a subset of points."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "groupA",  "Group A")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));
			
	parms.add(hutil::ParmFactory(PRM_STRING, "groupB",  "Group B")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "gridName", "Value Grid Name")
		.setDefault(0, "valueGrid")
        .setHelpText("value grid name."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "gradName", "Gradient Grid Name")
		.setDefault(0, "gradGrid")
        .setHelpText("gradient grid name."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "angleName", "Angle Grid Name")
		.setDefault(0, "angleGrid")
        .setHelpText("Angle grid name."));
	
	parms.add(hutil::ParmFactory(PRM_FLT_J, "error", "Error Bound")
        .setDefault(PRMzeroDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 0.1));
	
	parms.add(hutil::ParmFactory(PRM_INT, "iter", "Iterations")
        .setDefault(PRMfiveDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 30));

	parms.add(hutil::ParmFactory(PRM_FLT_J, "angle", "Stop Angle")
        .setDefault(PRM20Defaults)
        .setRange(PRM_RANGE_RESTRICTED, 10, PRM_RANGE_UI, 100));
	
	parms.add(hutil::ParmFactory(PRM_FLT_J, "bias", "Bias")
        .setDefault(PRMoneDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 1.5));
	
	parms.add(hutil::ParmFactory(PRM_INT, "relaxFreq", "Relax Frequency")
        .setDefault(PRMthreeDefaults)
        .setRange(PRM_RANGE_RESTRICTED, 3, PRM_RANGE_UI, 20));
	
    parms.add(hutil::ParmFactory(PRM_TOGGLE, "finalRelax", "Final relax")
        .setDefault(PRMoneDefaults));
	
    // Register this operator.
    hvdb::OpenVDBOpFactory("OpenVDB Vertex Proj", SOP_OpenVDB_Vertex_Proj::factory, parms, *table)
		.addInput("Polygon Mesh")
		.addInput("VDBs")
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
		const GU_Detail* vdbAGdp = inputGeo(1, context);
		const GU_Detail* vdbBGdp = inputGeo(2, context);
		
		mError = evalFloat("error", 0, time);
		mIterations = evalInt("iter", 0, time);
		mStopAngle = evalFloat("angle", 0, time);
		mBias = evalFloat("bias", 0, time);
		mRelaxFreq = evalInt("relaxFreq", 0, time);
		const bool doFinalRelax = evalInt("finalRelax", 0, time);
		
		UT_String gridNameStr;
		evalString(gridNameStr, "gridName", 0, time);
		
		UT_String gradNameStr;
		evalString(gradNameStr, "gradName", 0, time);
		
		UT_String angleNameStr;
		evalString(angleNameStr, "angleName", 0, time);
		
		// TODO only deal with some group of points
// 		UT_String groupPtsStr;
// 		evalString(groupPtsStr, "groupPts", 0, time);
// 		const GA_PrimitiveGroup *groupA = matchGroup(const_cast<GU_Detail&>(*vdbAGdp), groupPtsStr.toStdString());
		
		UT_String groupAStr;
		evalString(groupAStr, "groupA", 0, time);
		const GA_PrimitiveGroup *groupA = matchGroup(const_cast<GU_Detail&>(*vdbAGdp), groupAStr.toStdString());
		
		UT_String groupBStr;
		evalString(groupBStr, "groupB", 0, time);
		const GA_PrimitiveGroup *groupB = matchGroup(const_cast<GU_Detail&>(*vdbBGdp), groupBStr.toStdString());
		
// 		std::cout << "Required grid name : " << gridNameStr.toStdString() << std::endl;
// 		std::cout << "Required gradient name : " << gradNameStr.toStdString() << std::endl;
		const hvdb::GU_PrimVDB *pgrid = NULL, *pgrad = NULL;
		for (hvdb::VdbPrimCIterator it(vdbAGdp, groupA); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
// 			std::cout << "Grid name : " << gridName << std::endl;
			if ( gridName.compare( gridNameStr.toStdString() ) == 0 )
				pgrid = *it;
			if ( gridName.compare( gradNameStr.toStdString() ) == 0 )
				pgrad = *it;
		}
		
		const hvdb::GU_PrimVDB *pangle = NULL;
		for (hvdb::VdbPrimCIterator it(vdbBGdp, groupB); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
// 			std::cout << "Grid name : " << gridName << std::endl;
			if ( gridName.compare( angleNameStr.toStdString() ) == 0 )
				pangle = *it;
		}
		
		openvdb::FloatGrid::Ptr grid, angle;
		openvdb::VectorGrid::Ptr grad;
		if (pgrid != NULL && pgrad != NULL && pangle != NULL) {
			grid = openvdb::gridPtrCast<openvdb::FloatGrid>( pgrid->getGrid().deepCopyGrid() );
			grad = openvdb::gridPtrCast<openvdb::VectorGrid>( pgrad->getGrid().deepCopyGrid() );
			angle = openvdb::gridPtrCast<openvdb::FloatGrid>( pangle->getGrid().deepCopyGrid() );
		}
		else{
			throw std::runtime_error("Cannot find value or gradient grid.");
		}

		GA_ROHandleS	 neighbours(gdp->findStringTuple(GA_ATTRIB_POINT,"Neighbours",1));
		if (!neighbours.isValid())
			throw std::runtime_error("Cannot find attribute Neighbours");
		
		openvdb::math::Transform& gridform = grid->transform(); 		// gradient and grid have same transform
		openvdb::math::Transform& angleform = angle->transform();

		/*
		 * Project 
		 */
		hvdb::Interrupter progress("Projecting vertex based on vdb grid");
		
		//pHandle(gdp->findAttribute(GA_ATTRIB_POINT, "nearList"));
		GA_ROHandleR restVHandle(gdp->findAttribute(GA_ATTRIB_POINT, "restValue"));
		GA_RWHandleV3 pHandle(gdp->getP());	
		
		
		GA_Offset pointCount = 0;
		for (GA_Iterator it(gdp->getPointRange()); !it.atEnd(); it.advance()){
			pointCount++;
		}
			
			
		boost::dynamic_bitset<> negProj(pointCount);
		boost::dynamic_bitset<> contactRegion(pointCount);
		
		int step = 1;
		while( step <= mIterations ){
			if ( step % mRelaxFreq != 0 ) {
				bool allIsRightPos = true;
				
				for (GA_Iterator it(gdp->getPointRange()); !it.atEnd(); it.advance())
				{
					if(progress.wasInterrupted()) {
						throw std::runtime_error("Projection was interrupted");
					}
					
					openvdb::Vec3f vec, delta;
					float value, angleValue;
					UT_Vector3 pos;
					
					GA_Offset ptoff = it.getOffset();
					if (negProj[ptoff]) continue;		// have right position, so does not project any more
					
					pos = pHandle.get(ptoff);
					fpreal restv = restVHandle.get(ptoff);
					
					openvdb::Vec3f vpos = openvdb::Vec3f(pos.x(), pos.y(), pos.z());
					
					openvdb::Vec3f pos2 = angleform.worldToIndex(vpos);
					openvdb::tools::BoxSampler::sample(angle->tree(), pos2, angleValue);
					if ( angleValue > mStopAngle ) {
						contactRegion[ptoff] = 1;
						continue;
					}
					
					openvdb::Vec3f pos1 = gridform.worldToIndex(vpos);
					openvdb::tools::BoxSampler::sample(grid->tree(), pos1, value);
					openvdb::tools::BoxSampler::sample(grad->tree(), pos1, vec);
					if (fabs(value - restv) < mError){
						negProj[ptoff] = 1;
						continue;							
					}
					// find the point attr neighboutList and valueRest
					delta = vec * ( omega * (restv - value) / vec.length() );
					if (delta.length() < mError){
						negProj[ptoff] = 1;
						continue;							
					}
					
					pos.x() += delta.x();
					pos.y() += delta.y();
					pos.z() += delta.z();
					
					pHandle.set(ptoff, pos);	
					allIsRightPos = false;
				}
				if ( allIsRightPos ) 
					break;
			} // Projection step
			else {
				// relax step
				GA_Offset ptoff;
				UT_Vector3 finalPos;
				std::vector<UT_Vector3> finalPositions;
				float error, value, weight;
				
				for (GA_Iterator it(gdp->getPointRange()); !it.atEnd(); it.advance())
				{
					if (progress.wasInterrupted()) {
						throw std::runtime_error("Relax was interrupted");
					}
					
					ptoff = it.getOffset();
					
					const char *string_value = neighbours.get(*it);

					std::stringstream ss(string_value);
					std::vector<int> neighbourPts;
					int j;
					while (ss >> j) {
						neighbourPts.push_back( j );
						if (ss.peek() == ',')
							ss.ignore();
					}
					
					finalPos = UT_Vector3(0, 0, 0);
					int numberOfNeighbours = neighbourPts.size();
					
					for (unsigned int i = 0; i < numberOfNeighbours; i++){
						GA_Offset pt = neighbourPts[i];
						finalPos += gdp->getPos3(pt);
					}
					finalPos /= numberOfNeighbours;
					finalPositions.push_back(finalPos);
				}
				
				for (GA_Iterator it(gdp->getPointRange()); !it.atEnd(); it.advance())
				{
					if (progress.wasInterrupted()) {
						throw std::runtime_error("Relax was interrupted");
					}
					
					ptoff = it.getOffset();

					if (negProj[ptoff]) continue;		// have right position, so does not relax any more
					
					UT_Vector3 pos = pHandle.get(ptoff);
					openvdb::Vec3f vpos = openvdb::Vec3f(pos.x(), pos.y(), pos.z());
					openvdb::tools::BoxSampler::sample(grid->tree(), gridform.worldToIndex(vpos), value);
					
					error = value - restVHandle.get(ptoff);
					weight = 1- (error - mBias) * (error - mBias) * (error - mBias) * (error - mBias);
					weight = std::max(0.0f, weight);
					
					finalPos = finalPositions[ptoff];
// 					finalPos = ( 1 - weight ) * pos + weight * finalPos;
					finalPos = ( 1 - weight ) * pos + weight * finalPos;
					gdp->setPos3(ptoff, finalPos);
				}
					
			}// relax
			step++;
		}
		
		if ( doFinalRelax ) {
			// final relax , only work at contact region
			std::map<GA_Offset, UT_Vector3> contactPos;
			for (GA_Iterator it(gdp->getPointRange()); !it.atEnd(); it.advance())
			{
				if (progress.wasInterrupted()) {
					throw std::runtime_error("final relax was interrupted");
				}
						
				GA_Offset ptoff = it.getOffset();
				if (contactRegion[ptoff]) continue;		// not in contact region, just ignore

				const char *string_value = neighbours.get(*it);

				std::stringstream ss(string_value);
				std::vector<int> neighbourPts;
				int j;
				while (ss >> j) {
					neighbourPts.push_back( j );
					if (ss.peek() == ',')
						ss.ignore();
				}
				
				UT_Vector3 finalPos = UT_Vector3(0, 0, 0);
				int numberOfNeighbours = neighbourPts.size();
				
				for (unsigned int i = 0; i < numberOfNeighbours; i++){
					GA_Offset pt = neighbourPts[i];
					finalPos += gdp->getPos3(pt);
				}
				finalPos /= numberOfNeighbours;
				contactPos[ptoff] = finalPos;
			}
			
			for(std::map<GA_Offset, UT_Vector3>::iterator iter = contactPos.begin(); iter != contactPos.end(); ++iter)
			{
				GA_Offset ptoff =  iter->first;	
				gdp->setPos3(ptoff, iter->second);
			}
			
			contactPos.clear();
		}
		
	} catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
    }
    
    return error();
	
}
