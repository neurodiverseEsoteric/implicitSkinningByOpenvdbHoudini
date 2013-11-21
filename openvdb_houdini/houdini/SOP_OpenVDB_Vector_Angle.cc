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
	enum ResampleMode {
        RESAMPLE_OFF,    // don't auto-resample grids
        RESAMPLE_A,      // resample B to match A
        RESAMPLE_B     	// resample A to match B
    };
	
	SOP_OpenVDB_Vector_Angle(OP_Network *net, const char *name, OP_Operator *op);
    
	virtual ~SOP_OpenVDB_Vector_Angle(){};
	
    //virtual void getDescriptiveParmName(UT_String& s) const { s = "file_name"; }
    
	static OP_Node* factory(OP_Network*, const char* name, OP_Operator*);
	
	template<typename GridT>
	typename GridT::Ptr resampleToMatch(const GridT& src, const GridT& ref, int order);
	
protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);

};


////////////////////////////////////////


struct Local {
	static inline void dot(const openvdb::Vec3f& a, const openvdb::Vec3f& b, openvdb::Vec3f& result) {
		//result.x() =  acos( a.dot(b) ) / PI * 180;
		result.x() =  openvdb::math::angle(a, b) / PI * 180;
		//result.x() =  openvdb::math::angle(a, b);
	}
};


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
	
	
	parms.add(hutil::ParmFactory(PRM_STRING, "groupA",  "Group A")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "groupB",  "Group B")
		.setChoiceList(&hutil::PrimGroupMenu)
		.setHelpText("Choose only a subset of the input vdb grids."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "gradName", "Gradient Grid Name")
		.setDefault(0, "gradGrid")
        .setHelpText("gradient grid name."));
	
	parms.add(hutil::ParmFactory(PRM_STRING, "outGroup",  "output Group")
		.setDefault(0, "angleGrid")
		.setHelpText("Output vector angle grid name."));
	
	// Register this operator.
    hvdb::OpenVDBOpFactory("OpenVDB Vector Angle", SOP_OpenVDB_Vector_Angle::factory, parms, *table)
		.addInput("VDBA")
		.addInput("VDBB");
}


template<typename GridT>
typename GridT::Ptr SOP_OpenVDB_Vector_Angle::resampleToMatch(const GridT& src, const GridT& ref, int order)
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
SOP_OpenVDB_Vector_Angle::cookMySop(OP_Context &context)
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

		UT_String gradNameStr;
		evalString(gradNameStr, "gradName", 0, time);
		
		UT_String groupAStr;
		evalString(groupAStr, "groupA", 0, time);
		const GA_PrimitiveGroup *groupA = matchGroup(const_cast<GU_Detail&>(*vdbA), groupAStr.toStdString());
		
		UT_String groupBStr;
		evalString(groupBStr, "groupB", 0, time);
		const GA_PrimitiveGroup *groupB = matchGroup(const_cast<GU_Detail&>(*vdbB), groupBStr.toStdString());
		
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
		
		openvdb::VectorGrid::Ptr gradA, gradB;
		if ( pgradA != NULL && pgradB != NULL ) {
			gradA = openvdb::gridPtrCast<openvdb::VectorGrid>( pgradA->getGrid().deepCopyGrid() );	
			gradB = openvdb::gridPtrCast<openvdb::VectorGrid>( pgradB->getGrid().deepCopyGrid() );	
		}
		else{
			throw std::runtime_error("Cannot find gradient grid.");
		}
		
		/*
		 *  Calculate the angle between two gradient vector field
		 */
		hvdb::Interrupter progress("Projecting vertex based on vdb grid");
		
		// Create a group for the grid primitives.
        GA_PrimitiveGroup* group = NULL;
        UT_String groupStr;
        evalString(groupStr, "outGroup", 0, time);
        if(groupStr.isstring()) {
            group = gdp->newPrimitiveGroup(groupStr.buffer());
        }
        
		openvdb::VectorGrid::Ptr angleGrid = openvdb::VectorGrid::create();
		openvdb::FloatGrid::Ptr outGrid = openvdb::FloatGrid::create();
		
		const openvdb::math::Transform
			&sourceXform = gradA->transform(),
			&targetXform = gradB->transform();
			// Compute a source grid to target grid transform.
			// (For this example, we assume that both grids' transforms are linear,
			// so that they can be represented as 4 x 4 matrices.)
		openvdb::Mat4R xformA =
			sourceXform.baseMap()->getAffineMap()->getMat4();
		std::cout << "grid a' transform is " << xformA << std::endl;
		
		openvdb::Mat4R xformB =
			targetXform.baseMap()->getAffineMap()->getMat4();
		std::cout << "grid b' transform is " << xformB << std::endl;
		
		
		if (gradA->constTransform() != gradB->constTransform()) {
			std::cout << "grid a' transform != grid b' transform" << std::endl;
			
			int resampleWhich = RESAMPLE_OFF;
			int samplingOrder = 1;
			
			const openvdb::Vec3d aVoxSize = gradA->voxelSize(), bVoxSize = gradA->voxelSize();
			const double aVoxVol = aVoxSize[0] * aVoxSize[1] * aVoxSize[2], bVoxVol = bVoxSize[0] * bVoxSize[1] * bVoxSize[2];
			
			resampleWhich = (aVoxVol < bVoxVol) ? RESAMPLE_A : RESAMPLE_B;
			
			openvdb::VectorGrid::Ptr resampledGrid;
// 			if (gradA->getGridClass() == openvdb::GRID_LEVEL_SET){
// 				std::cout << "grid a' is level set grid" << std::endl;
// 			}
			if (resampleWhich == RESAMPLE_A) {
				std::cout << "Resample grid A" << std::endl;
				resampledGrid = this->resampleToMatch(*gradA, *gradB, samplingOrder);	
				resampledGrid->setTransform(gradB->transform().copy());
				angleGrid->tree().combine2(resampledGrid->tree(), gradB->tree(), Local::dot, false); 
			} else {
				std::cout << "Resample grid B" << std::endl;
				resampledGrid = this->resampleToMatch(*gradB, *gradA, samplingOrder);
				resampledGrid->setTransform(gradA->transform().copy());
				angleGrid->setTransform(resampledGrid->transform().copy());
				outGrid->setTransform(resampledGrid->transform().copy());
				
				const openvdb::math::Transform &Xform = angleGrid->transform();
				openvdb::Mat4R xform = Xform.baseMap()->getAffineMap()->getMat4();
				std::cout << "grid angle' transform is " << xform << std::endl;
		
// 				angleGrid->setTransform(openvdb::math::Transform::createLinearTransform(resampledGrid->voxelSize().x()));
				angleGrid->tree().combine2(gradA->tree(), resampledGrid->tree(), Local::dot, false); 
// 				gradB = resampledGrid.get();
			}
		}
		else {
			std::cout << "grid a' transform == grid b' transform" << std::endl;
			angleGrid->setTransform(gradA->transform().copy());
			angleGrid->tree().combine2(gradA->tree(), gradB->tree(), Local::dot, false); 
		}

		openvdb::FloatGrid::Accessor accessor = outGrid->getAccessor();
		
		for (openvdb::VectorGrid::ValueOnCIter iter = angleGrid->cbeginValueOn(); iter.test(); ++iter) {
			const openvdb::Vec3f& value = *iter;
			if (iter.isVoxelValue()) { // set a single voxel
				accessor.setValue(iter.getCoord(), value.x());
			} else { // fill an entire tile
				openvdb::CoordBBox bbox;
				iter.getBoundingBox(bbox);
				accessor.getTree()->fill(bbox, value.x());
			}
		}
		
		std::cout << "gradient dot angle values are : " << std::endl;
		for (openvdb::FloatGrid::ValueOnCIter iter = outGrid->cbeginValueOn(); iter; ++iter) {
			float value = iter.getValue();
// 			if ( value.x() > 0 ) 
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
