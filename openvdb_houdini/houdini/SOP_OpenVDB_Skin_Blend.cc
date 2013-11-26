///////////////////////////////////////////////////////////////////////////
//
/// @file SOP_OpenVDB_Skin_Blend.cc
///
/// @author maphysart
///
/// @brief skin blend two field based on gradient angle

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

#include <algorithm>    // std::max
#include <math.h>       /* fabs */

namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;

class SOP_OpenVDB_Skin_Blend : public hvdb::SOP_NodeVDB
{
public:	
	enum ResampleMode {
        RESAMPLE_OFF,    // don't auto-resample grids
        RESAMPLE_A,      // resample B to match A
        RESAMPLE_B     	// resample A to match B
    };
	
	SOP_OpenVDB_Skin_Blend(OP_Network *net, const char *name, OP_Operator *op);
    
	virtual ~SOP_OpenVDB_Skin_Blend(){};
	
    //virtual void getDescriptiveParmName(UT_String& s) const { s = "file_name"; }
    
	static OP_Node* factory(OP_Network*, const char* name, OP_Operator*);
	
protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);
	
	template<typename GridT>
	typename GridT::Ptr resampleToMatch(const GridT& src, const GridT& ref, int order);

	float blend( float a, float b, float t); //look for blend grid table
};


////////////////////////////////////////

namespace {
	struct Local {
		static inline void max(const float& a, const float& b, float& result) {
			//result.x() =  acos( a.dot(b) ) / PI * 180;
			//result.x() =  openvdb::math::angle(a, b) / PI * 180;
			result =  std::max(fabs(a), fabs(b));
		}
	};
}

OP_Node* 
SOP_OpenVDB_Skin_Blend::factory(OP_Network* net,
    const char* name, OP_Operator* op)
{
    return new SOP_OpenVDB_Skin_Blend(net, name, op);
}

SOP_OpenVDB_Skin_Blend::SOP_OpenVDB_Skin_Blend(OP_Network* net,
    const char* name, OP_Operator* op) : hvdb::SOP_NodeVDB(net, name, op)
{
}


// Build UI and register this operator.
void
newSopOperator(OP_OperatorTable* table)
{
    if (table == NULL) return;

    hutil::ParmList parms;
	
	
	parms.add(hutil::ParmFactory(PRM_STRING, "groupA",  "Group A")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "groupB",  "Group B")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "groupC",  "Group C")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "valueName", "Value Grid Name")
		.setDefault(0, "valueGrid")
        .setHelpText("value grid name."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "thetaName", "Theta Grid Name")
		.setDefault(0, "angleGrid")
        .setHelpText("theta grid name."));
	
	parms.add(hutil::ParmFactory(PRM_FILE, "vdbFileName", "Blend Grid VDB")
        .setDefault(0, "./output.vdb")
        .setHelpText("Path name for the output VDB file"));

	parms.add(hutil::ParmFactory(PRM_STRING, "blendName", "blend Grid Name")
		.setDefault(0, "blendGrid")
        .setHelpText("blend grid name."));
		
	parms.add(hutil::ParmFactory(PRM_STRING, "outGroup",  "output Group")
		.setDefault(0, "fieldGrid")
		.setHelpText("Output vector angle grid name."));
	
	// Register this operator.
    hvdb::OpenVDBOpFactory("OpenVDB Vector Angle", SOP_OpenVDB_Skin_Blend::factory, parms, *table)
		.addInput("VDBA")
		.addInput("VDBB")
		.addInput("VDBC");
}


template<typename GridT>
typename GridT::Ptr SOP_OpenVDB_Skin_Blend::resampleToMatch(const GridT& src, const GridT& ref, int order)
{
	hvdb::Interrupter interrupt = hvdb::Interrupter();
	
	typedef typename GridT::ValueType ValueT;

	const openvdb::math::Transform& refXform = ref.constTransform();

	typename GridT::Ptr dest;
	if (src.getGridClass() == openvdb::GRID_LEVEL_SET) {
		// For level set grids, use the level set rebuild tool to both resample the
		// source grid to match the reference grid and to rebuild the resulting level set.
		const ValueT halfWidth = ((ref.getGridClass() == openvdb::GRID_LEVEL_SET)
			? ValueT(ref.background() * (1.0 / ref.voxelSize()[0]))
			: ValueT(src.background() * (1.0 / src.voxelSize()[0])));
		try {
			dest = openvdb::tools::doLevelSetRebuild(src, /*iso=*/openvdb::zeroVal<ValueT>(),
				/*exWidth=*/halfWidth, /*inWidth=*/halfWidth, &refXform, &interrupt);
		} catch (openvdb::TypeError&) {
			this->addWarning(SOP_MESSAGE, ("skipped rebuild of level set grid "
				+ src.getName() + " of type " + src.type()).c_str());
			dest.reset();
		}
	}
	if (!dest && src.constTransform() != refXform) {
		// For non-level set grids or if level set rebuild failed due to an unsupported
		// grid type, use the grid transformer tool to resample the source grid to match
		// the reference grid.
		dest = src.copy(openvdb::CP_NEW);
		dest->setTransform(refXform.copy());
		using namespace openvdb;
		switch (order) {
		case 0: tools::resampleToMatch<tools::PointSampler>(src, *dest, interrupt); break;
		case 1: tools::resampleToMatch<tools::BoxSampler>(src, *dest, interrupt); break;
		case 2: tools::resampleToMatch<tools::QuadraticSampler>(src, *dest, interrupt); break;
		}
	}
	return dest;
}
    

OP_ERROR
SOP_OpenVDB_Skin_Blend::cookMySop(OP_Context &context)
{
	try {
        hutil::ScopedInputLock lock(*this, context);
		//duplicatePointSource(0, context);
        const fpreal time = context.getTime();
		gdp->clearAndDestroy();
		
		/*
		 * Get params
		 */
		const GU_Detail* vdbA = inputGeo(0, context);
		const GU_Detail* vdbB = inputGeo(1, context);
		const GU_Detail* vdbC = inputGeo(2, context);
		
		UT_String valueNameStr;
		evalString(valueNameStr, "valueName", 0, time);

		UT_String thetaNameStr;
		evalString(thetaNameStr, "thetaName", 0, time);
		
		UT_String groupAStr;
		evalString(groupAStr, "groupA", 0, time);
		const GA_PrimitiveGroup *groupA = matchGroup(const_cast<GU_Detail&>(*vdbA), groupAStr.toStdString());
		
		UT_String groupBStr;
		evalString(groupBStr, "groupB", 0, time);
		const GA_PrimitiveGroup *groupB = matchGroup(const_cast<GU_Detail&>(*vdbB), groupBStr.toStdString());
		
		UT_String groupCStr;
		evalString(groupCStr, "groupC", 0, time);
		const GA_PrimitiveGroup *groupC = matchGroup(const_cast<GU_Detail&>(*vdbC), groupCStr.toStdString());
		
		//std::cout << "Required gradient name : " << gradNameStr.toStdString() << std::endl;
		const hvdb::GU_PrimVDB *pgradA = NULL;
		for (hvdb::VdbPrimCIterator it(vdbA, groupA); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
			if ( gridName.compare( valueNameStr.toStdString() ) == 0 )
				pgradA = *it;
		}
		
		const hvdb::GU_PrimVDB *pgradB = NULL;
		for (hvdb::VdbPrimCIterator it(vdbB, groupB); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
			if ( gridName.compare( valueNameStr.toStdString() ) == 0 )
				pgradB = *it;
		}
		
		const hvdb::GU_PrimVDB *pgradC = NULL;
		for (hvdb::VdbPrimCIterator it(vdbC, groupC); it; ++it) {
			const std::string gridName = it.getPrimitiveName().toStdString();
			if ( gridName.compare( valueNameStr.toStdString() ) == 0 )
				pgradC = *it;
		}
		
		openvdb::FloatGrid::Ptr gridA, gridB, gridC;
		if ( pgradA != NULL && pgradB != NULL ) {
			gridA = openvdb::gridPtrCast<openvdb::FloatGrid>( pgradA->getGrid().deepCopyGrid() );	
			gridB = openvdb::gridPtrCast<openvdb::FloatGrid>( pgradB->getGrid().deepCopyGrid() );	
			gridC = openvdb::gridPtrCast<openvdb::FloatGrid>( pgradC->getGrid().deepCopyGrid() );	
		}
		else{
			throw std::runtime_error("Cannot find value grid.");
		}
		
		/*
		 *  Calculate the angle between two gradient vector field
		 */
		hvdb::Interrupter progress("Calculate vector angle based on two gradient vector grid");
				
		UT_String fileNameStr;
		evalString(fileNameStr, "vdbFileName", 0, time);
		const std::string filename = fileNameStr.toStdString();
		if (filename.empty()) {
			addWarning(SOP_MESSAGE, "no name given for the output file");
			return error();
		}
		
		// Create a group for the grid primitives.
        GA_PrimitiveGroup* group = NULL;
        UT_String groupStr;
        evalString(groupStr, "outGroup", 0, time);
        if(groupStr.isstring()) {
            group = gdp->newPrimitiveGroup(groupStr.buffer());
        }
        
		openvdb::FloatGrid::Ptr outGrid = openvdb::FloatGrid::create();

		const openvdb::math::Transform
			&sourceXform = gridA->transform(),
			&targetXform = gridB->transform();
			// Compute a source grid to target grid transform.
			// (For this example, we assume that both grids' transforms are linear,
			// so that they can be represented as 4 x 4 matrices.)
		openvdb::Mat4R xformA =
			sourceXform.baseMap()->getAffineMap()->getMat4();
		std::cout << "grid a' transform is " << xformA << std::endl;
		
		openvdb::Mat4R xformB =
			targetXform.baseMap()->getAffineMap()->getMat4();
		std::cout << "grid b' transform is " << xformB << std::endl;
		
		
		if (gridA->constTransform() != gridB->constTransform()) {
			int resampleWhich = RESAMPLE_OFF;
			int samplingOrder = 1;
			
			const openvdb::Vec3d aVoxSize = gridA->voxelSize(), bVoxSize = gridB->voxelSize();
			const double aVoxVol = aVoxSize[0] * aVoxSize[1] * aVoxSize[2], bVoxVol = bVoxSize[0] * bVoxSize[1] * bVoxSize[2];
			
			resampleWhich = (aVoxVol < bVoxVol) ? RESAMPLE_A : RESAMPLE_B;
			
			openvdb::FloatGrid::Ptr resampledGrid;
			
			if (resampleWhich == RESAMPLE_A) {
				resampledGrid = this->resampleToMatch(*gridA, *gridB, samplingOrder);	
				resampledGrid->setTransform(gridB->transform().copy());
				outGrid->setTransform(resampledGrid->transform().copy());
				outGrid->tree().combine2(resampledGrid->tree(), gridB->tree(), Local::max, false); 
			} else {
				resampledGrid = this->resampleToMatch(*gridB, *gridA, samplingOrder);
				resampledGrid->setTransform(gridA->transform().copy());
				outGrid->setTransform(resampledGrid->transform().copy());
				outGrid->tree().combine2(gridA->tree(), resampledGrid->tree(), Local::max, false); 
			}
		}
		else {
			outGrid->setTransform(gridA->transform().copy());
			outGrid->tree().combine2(gridA->tree(), gridB->tree(), Local::max, false); 
		}
		
		
		/*
		 * read blend grid from file on disk
		 */
		UT_String blendNameStr;
		evalString(blendNameStr, "blendName", 0, time);
		
		openvdb::FloatGrid *  blendPtr = NULL;
		
		openvdb::io::File file(filename);
        file.open();
        // Loop over all grids in the file.		
		for (openvdb::io::File::NameIterator nameIter = file.beginName();
            nameIter != file.endName(); ++nameIter)
        {
			const std::string& gridName = nameIter.gridName();
			if (!UT_String(gridName).multiMatch(blendNameStr.buffer(), 1, " ")) continue;
			
			hvdb::GridPtr blendGrid = file.readGrid(gridName);
	
            if (blendGrid) {
				blendPtr = static_cast<openvdb::FloatGrid *>(blendGrid.get());
			}
			
		}
		if (blendPtr == NULL){
			throw std::runtime_error("Cannot find blend grid file.");
		}	
		openvdb::math::Transform& gridform = blendPtr->transform();
		

		/*
		 * blend two grids based on gradient
		 */
		openvdb::FloatGrid::Accessor accessorA = gridA->getAccessor();
		openvdb::FloatGrid::Accessor accessorB = gridB->getAccessor();
		openvdb::FloatGrid::Accessor accessorC = gridC->getAccessor();
		
		for (openvdb::FloatGrid::ValueOnIter iter = outGrid->beginValueOn(); iter.test(); ++iter) {
			openvdb::Coord xyz = iter.getCoord();
			float a, b, t, value;
			if (iter.isVoxelValue()) { // set a single voxel
				a = accessorA.getValue(xyz);
				b = accessorB.getValue(xyz);
				t = accessorC.getValue(xyz);
				
				openvdb::Vec3f vpos = openvdb::Vec3f(a, b, t);
				vpos = gridform.worldToIndex(vpos);
				openvdb::tools::BoxSampler::sample(blendPtr->tree(), vpos, value);
				
				iter.setValue(value);
				
			}
		}
		
		std::cout << "blend field are : " << std::endl;
		for (openvdb::FloatGrid::ValueOnCIter iter = outGrid->cbeginValueOn(); iter; ++iter) {
			float value = iter.getValue();
			std::cout << "Grid world" << outGrid->constTransform().indexToWorld(iter.getCoord()) << " index" << iter.getCoord() << " = " << value << std::endl;
		}
		
		// use the same name as the output group name 
		GEO_PrimVDB* vdb = hvdb::createVdbPrimitive(*gdp, outGrid, groupStr.toStdString().c_str());
		if (group) group->add(vdb);
		
	} catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
    }
    
    return error();
	
}
