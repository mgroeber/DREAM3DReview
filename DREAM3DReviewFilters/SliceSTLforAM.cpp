/* ============================================================================
* Software developed by US federal government employees (including military personnel) 
* as part of their official duties is not subject to copyright protection and is 
* considered “public domain” (see 17 USC Section 105). Public domain software can be used 
* by anyone for any purpose, and cannot be released under a copyright license 
* (including typical open source software licenses).
* 
* This source code file was originally written by United States DoD employees. The
* original source code files are released into the Public Domain.
* 
* Subsequent changes to the codes by others may elect to add a copyright and license
* for those changes.
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include <memory>

#include "SliceSTLforAM.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Math/GeometryMath.h"
#include "SIMPLib/Math/MatrixMath.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/Core/OrientationTransformation.hpp"
#include "EbsdLib/Core/Quaternion.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SliceSTLforAM::SliceSTLforAM()
: m_CADDataContainerName("TriangleDataContainer")
, m_SliceDataContainerName("SliceDataContainer")
, m_EdgeAttributeMatrixName("EdgeData")
, m_SliceAttributeMatrixName("SliceData")
, m_SliceIdArrayName("SliceIds")
, m_AreasArrayName("SliceAreas")
, m_PerimetersArrayName("SlicePerimeters")
, m_HaveRegionIds(false)
, m_RegionIdArrayPath("", "", "")
, m_SliceResolution(1.0f)
, m_Zstart(0.0)
, m_Zend(0.0)
, m_SliceRange(0)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SliceSTLforAM::~SliceSTLforAM() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceSTLforAM::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  {
    QVector<QString> choices;
    choices.push_back("Full Range");
    choices.push_back("User Defined Range");

    QStringList linkedChoiceProps;
    linkedChoiceProps << "Zstart"
                      << "Zend";

    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Slice Range");
    parameter->setPropertyName("SliceRange");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(SliceSTLforAM, this, SliceRange));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(SliceSTLforAM, this, SliceRange));
    parameter->setChoices(choices);
    parameter->setLinkedProperties(linkedChoiceProps);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slicing Start", Zstart, FilterParameter::Parameter, SliceSTLforAM, 1));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slicing End", Zend, FilterParameter::Parameter, SliceSTLforAM, 1));

  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slice Spacing", SliceResolution, FilterParameter::Parameter, SliceSTLforAM));

  QStringList linkedProps("RegionIdArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Have Region Ids", HaveRegionIds, FilterParameter::Parameter, SliceSTLforAM, linkedProps));
  linkedProps.clear();
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Triangle};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("CAD Geometry", CADDataContainerName, FilterParameter::RequiredArray, SliceSTLforAM, dcsReq));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Face, IGeometry::Type::Triangle);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Region Ids", RegionIdArrayPath, FilterParameter::RequiredArray, SliceSTLforAM, req));
  }
  parameters.push_back(SIMPL_NEW_STRING_FP("Slice Geometry", SliceDataContainerName, FilterParameter::CreatedArray, SliceSTLforAM));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Edge Attribute Matrix", EdgeAttributeMatrixName, SliceDataContainerName, FilterParameter::CreatedArray, SliceSTLforAM));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Slice Ids", SliceIdArrayName, SliceDataContainerName, EdgeAttributeMatrixName, FilterParameter::CreatedArray, SliceSTLforAM));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Slice Attribute Matrix", SliceAttributeMatrixName, SliceDataContainerName, FilterParameter::CreatedArray, SliceSTLforAM));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Areas", AreasArrayName, SliceDataContainerName, SliceAttributeMatrixName, FilterParameter::CreatedArray, SliceSTLforAM));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Perimeters", PerimetersArrayName, SliceDataContainerName, SliceAttributeMatrixName, FilterParameter::CreatedArray, SliceSTLforAM));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceSTLforAM::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setCADDataContainerName(reader->readDataArrayPath("CADDataContainerName", getCADDataContainerName()));
  setSliceDataContainerName(reader->readString("SliceDataContainerName", getSliceDataContainerName()));
  setEdgeAttributeMatrixName(reader->readString("EdgeAttributeMatrixName", getEdgeAttributeMatrixName()));
  setSliceAttributeMatrixName(reader->readString("SliceAttributeMatrixName", getSliceAttributeMatrixName()));
  setSliceIdArrayName(reader->readString("SliceIdArrayName", getSliceIdArrayName()));
  setHaveRegionIds(reader->readValue("HaveRegionIds", getHaveRegionIds()));
  setRegionIdArrayPath(reader->readDataArrayPath("RegionIdArrayPath", getRegionIdArrayPath()));
  setSliceResolution(reader->readValue("SliceResolution", getSliceResolution()));
  setSliceRange(reader->readValue("SliceRange", getSliceRange()));
  setZstart(reader->readValue("Zstart", getZstart()));
  setZend(reader->readValue("Zend", getZend()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceSTLforAM::updateSliceInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_AreaPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Area = m_AreaPtr.lock()->getPointer(0);
  }                                          /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_PerimeterPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Perimeter = m_PerimeterPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceSTLforAM::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  DataArrayPath tempPath;
  QVector<IDataArray::Pointer> dataArrays;

  TriangleGeom::Pointer triangle = getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom>(this, getCADDataContainerName());

  if(m_SliceRange == 1)
  {
    if(m_Zstart >= m_Zend)
    {
      QString message = QObject::tr("Z end must be larger than Z start.");
      setErrorCondition(-62102, message);
    }
  }

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer(this, getSliceDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }
  SharedVertexList::Pointer vertices = EdgeGeom::CreateSharedVertexList(0);
  EdgeGeom::Pointer edge = EdgeGeom::CreateGeometry(0, vertices, SIMPL::Geometry::EdgeGeometry, !getInPreflight());
  m->setGeometry(edge);

  std::vector<size_t> tDims(1, 0);
  m->createNonPrereqAttributeMatrix(this, getEdgeAttributeMatrixName(), tDims, AttributeMatrix::Type::Edge);
  m->createNonPrereqAttributeMatrix(this, getSliceAttributeMatrixName(), tDims, AttributeMatrix::Type::EdgeFeature);
  if(getErrorCode() < 0)
  {
    return;
  }

  std::vector<size_t> cDims(1, 1);
  if(m_HaveRegionIds)
  {
    m_TriRegionIdPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getRegionIdArrayPath(),
                                                                                                           cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_TriRegionIdPtr.lock().get())                                                                   /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_TriRegionId = m_TriRegionIdPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() >= 0)
    {
      dataArrays.push_back(m_TriRegionIdPtr.lock());
    }

    tempPath.update(getSliceDataContainerName(), getEdgeAttributeMatrixName(), getRegionIdArrayPath().getDataArrayName());
    m_RegionIdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(
        this, tempPath, 0, cDims);            /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_RegionIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_RegionId = m_RegionIdPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() < 0)
    {
      return;
    }
  }

  tempPath.update(getSliceDataContainerName(), getEdgeAttributeMatrixName(), getSliceIdArrayName());
  m_SliceIdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0,
                                                                                                                    cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_SliceIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_SliceId = m_SliceIdPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }

  tempPath.update(getSliceDataContainerName(), getSliceAttributeMatrixName(), getAreasArrayName());
  m_AreaPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0,
                                                                                                               cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_AreaPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Area = m_AreaPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }

  tempPath.update(getSliceDataContainerName(), getSliceAttributeMatrixName(), getPerimetersArrayName());
  m_PerimeterPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0,
                                                                                                                    cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_PerimeterPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Perimeter = m_PerimeterPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceSTLforAM::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }


  SliceTriangleGeometry::Pointer sliceTriGeom = SliceTriangleGeometry::New();
  sliceTriGeom->setDataContainerArray(getDataContainerArray());
  sliceTriGeom->setDoRotation(false);
  sliceTriGeom->setCADDataContainerName(getCADDataContainerName());
  sliceTriGeom->setEdgeAttributeMatrixName(getEdgeAttributeMatrixName());
  sliceTriGeom->setHaveRegionIds(getHaveRegionIds());
  sliceTriGeom->setRegionIdArrayPath(getRegionIdArrayPath());
  sliceTriGeom->setSliceAttributeMatrixName(getSliceAttributeMatrixName());
  sliceTriGeom->setSliceDataContainerName(getSliceDataContainerName());

  sliceTriGeom->setSliceIdArrayName(getSliceIdArrayName());
  sliceTriGeom->setSliceRange(getSliceRange());
  sliceTriGeom->setSliceResolution(getSliceResolution());
  sliceTriGeom->setZend(getZend());
  sliceTriGeom->setZstart(getZstart());

  TriangleGeom::Pointer triangle = getDataContainerArray()->getDataContainer(getCADDataContainerName())->getGeometryAs<TriangleGeom>();

  MeshIndexType* tris = triangle->getTriPointer(0);
  float* triVerts = triangle->getVertexPointer(0);
  MeshIndexType numTris = triangle->getNumberOfTris();
  MeshIndexType numTriVerts = triangle->getNumberOfVertices();

  //// calculate slice areas before rotating back so we can use simple area calculation
  //for(size_t i = 0; i < numEdges; i++)
  //{
  //  int32_t sliceId = m_SliceId[i];
  //  float height = 0.5 * (verts[3 * (2 * i) + 1] + verts[3 * (2 * i + 1) + 1]);
  //  float width = (verts[3 * (2 * i + 1)] - verts[3 * (2 * i)]);
  //  float length = ((verts[3 * (2 * i + 1)] - verts[3 * (2 * i)]) * (verts[3 * (2 * i + 1)] - verts[3 * (2 * i)])) +
  //                 ((verts[3 * (2 * i + 1) + 1] - verts[3 * (2 * i) + 1]) * (verts[3 * (2 * i + 1) + 1] - verts[3 * (2 * i) + 1])) +
  //                 ((verts[3 * (2 * i + 1) + 2] - verts[3 * (2 * i) + 2]) * (verts[3 * (2 * i + 1) + 2] - verts[3 * (2 * i) + 2]));
  //  length = sqrt(length);
  //  float area = height * width;
  //  m_Area[sliceId] += area;
  //  m_Perimeter[sliceId] += length;
  //}

  // take absolute value to ensure areas are positive (in case winding is such that areas come out negative)
  for(size_t i = 0; i < m_NumberOfSlices; i++)
  {
    m_Area[i] = fabsf(m_Area[i]);
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer SliceSTLforAM::newFilterInstance(bool copyFilterParameters) const
{
  SliceSTLforAM::Pointer filter = SliceSTLforAM::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceSTLforAM::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceSTLforAM::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceSTLforAM::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceSTLforAM::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceSTLforAM::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceSTLforAM::getHumanLabel() const
{
  return "Slice CAD Geometry";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid SliceSTLforAM::getUuid() const
{
  return QUuid("{222307a4-67fd-5cb5-a12e-d80f9fb970af}");
}

// -----------------------------------------------------------------------------
SliceSTLforAM::Pointer SliceSTLforAM::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<SliceSTLforAM> SliceSTLforAM::New()
{
  struct make_shared_enabler : public SliceSTLforAM
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::getNameOfClass() const
{
  return QString("SliceSTLforAM");
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::ClassName()
{
  return QString("SliceSTLforAM");
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setCADDataContainerName(const DataArrayPath& value)
{
  m_CADDataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath SliceSTLforAM::getCADDataContainerName() const
{
  return m_CADDataContainerName;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setSliceDataContainerName(const QString& value)
{
  m_SliceDataContainerName = value;
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::getSliceDataContainerName() const
{
  return m_SliceDataContainerName;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setEdgeAttributeMatrixName(const QString& value)
{
  m_EdgeAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::getEdgeAttributeMatrixName() const
{
  return m_EdgeAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setSliceAttributeMatrixName(const QString& value)
{
  m_SliceAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::getSliceAttributeMatrixName() const
{
  return m_SliceAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setSliceIdArrayName(const QString& value)
{
  m_SliceIdArrayName = value;
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::getSliceIdArrayName() const
{
  return m_SliceIdArrayName;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setAreasArrayName(const QString& value)
{
  m_AreasArrayName = value;
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::getAreasArrayName() const
{
  return m_AreasArrayName;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setPerimetersArrayName(const QString& value)
{
  m_PerimetersArrayName = value;
}

// -----------------------------------------------------------------------------
QString SliceSTLforAM::getPerimetersArrayName() const
{
  return m_PerimetersArrayName;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setHaveRegionIds(bool value)
{
  m_HaveRegionIds = value;
}

// -----------------------------------------------------------------------------
bool SliceSTLforAM::getHaveRegionIds() const
{
  return m_HaveRegionIds;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setRegionIdArrayPath(const DataArrayPath& value)
{
  m_RegionIdArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath SliceSTLforAM::getRegionIdArrayPath() const
{
  return m_RegionIdArrayPath;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setSliceResolution(float value)
{
  m_SliceResolution = value;
}

// -----------------------------------------------------------------------------
float SliceSTLforAM::getSliceResolution() const
{
  return m_SliceResolution;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setZstart(float value)
{
  m_Zstart = value;
}

// -----------------------------------------------------------------------------
float SliceSTLforAM::getZstart() const
{
  return m_Zstart;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setZend(float value)
{
  m_Zend = value;
}

// -----------------------------------------------------------------------------
float SliceSTLforAM::getZend() const
{
  return m_Zend;
}

// -----------------------------------------------------------------------------
void SliceSTLforAM::setSliceRange(int value)
{
  m_SliceRange = value;
}

// -----------------------------------------------------------------------------
int SliceSTLforAM::getSliceRange() const
{
  return m_SliceRange;
}
