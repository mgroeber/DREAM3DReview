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

#include "LaplacianSmoothPointCloud.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/Math/MatrixMath.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
LaplacianSmoothPointCloud::LaplacianSmoothPointCloud() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
LaplacianSmoothPointCloud::~LaplacianSmoothPointCloud() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  QStringList linkedProps("MaskArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, LaplacianSmoothPointCloud, linkedProps));
  linkedProps.clear();
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Iterations", NumIterations, FilterParameter::Parameter, LaplacianSmoothPointCloud));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Lambda", Lambda, FilterParameter::Parameter, LaplacianSmoothPointCloud));
  {
    DataContainerSelectionFilterParameter::RequirementType req;
    IGeometry::Types geomTypes = {IGeometry::Type::Vertex};
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Data Container to Smooth", DataContainerName, FilterParameter::RequiredArray, LaplacianSmoothPointCloud, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, LaplacianSmoothPointCloud, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setDataContainerName(reader->readDataArrayPath("DataContainerName", getDataContainerName()));
  setNumIterations(reader->readValue("NumIterations", getNumIterations()));
  setLambda(reader->readValue("Lambda", getLambda()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  QVector<IDataArray::Pointer> dataArrays;

  VertexGeom::Pointer vertices = getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, getDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }
  dataArrays.push_back(vertices->getVertices());

  if(getLambda() <= 0)
  {
    QString ss = QObject::tr("Number of Iterations must be greater than 0");
    setErrorCondition(-11000, ss);
  }
  if(getLambda() <= 0 || getLambda() > 1)
  {
    QString ss = QObject::tr("Lambda must be greater than 0 and less than or equal to 1");
    setErrorCondition(-11000, ss);
  }
  if(getErrorCode() < 0)
  {
    return;
  }

  if(getUseMask() == true)
  {
    std::vector<size_t> cDims(1, 1);
    m_MaskPtr =
        getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>>(this, getMaskArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_MaskPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() >= 0)
    {
      dataArrays.push_back(m_MaskPtr.lock());
    }
  }

  getDataContainerArray()->validateNumberOfTuples(this, dataArrays);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::determineNewCoords(float* vertex, int64_t vertId, float* newCoords)
{
  float halfCoord[3];
  for(size_t j = 0; j < 3; j++)
  {
    halfCoord[j] = (vertex[3 * (vertId + 1) + j] + vertex[3 * (vertId - 1) + j]) * 0.5f;
  }
  for(size_t j = 0; j < 3; j++)
  {
    newCoords[3 * vertId + j] = vertex[3 * vertId + j] + (m_Lambda * (halfCoord[j] - vertex[3 * vertId + j]));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  VertexGeom::Pointer vertices = getDataContainerArray()->getDataContainer(getDataContainerName())->getGeometryAs<VertexGeom>();

  int64_t numVerts = vertices->getNumberOfVertices();
  float* vertex = vertices->getVertexPointer(0);

  int64_t progressInt = 0;

  FloatArrayType::Pointer newCoordsPtr = FloatArrayType::CreateArray(3 * numVerts, "newCoords", true);
  newCoordsPtr->initializeWithValue(0.0f);
  float* newCoords = newCoordsPtr->getPointer(0);

  for(int64_t iter = 0; iter < m_NumIterations; iter++)
  {
    progressInt = static_cast<int64_t>((static_cast<float>(iter) / m_NumIterations) * 100.0f);
    QString ss = QObject::tr("Smoothing Point Cloud || %1% Completed").arg(progressInt);
    notifyStatusMessage(ss);
    for(int64_t i = 0; i < numVerts - 1; i++)
    {
      if(!m_UseMask || (m_UseMask && m_Mask[i]))
      {
        determineNewCoords(vertex, i, newCoords);
      }
    }
    for(int64_t i = 0; i < numVerts - 1; i++)
    {
      if(!m_UseMask || (m_UseMask && m_Mask[i]))
      {
        vertex[3 * i + 0] = newCoords[3 * i + 0];
        vertex[3 * i + 1] = newCoords[3 * i + 1];
        vertex[3 * i + 2] = newCoords[3 * i + 2];
      }
    }
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer LaplacianSmoothPointCloud::newFilterInstance(bool copyFilterParameters) const
{
  LaplacianSmoothPointCloud::Pointer filter = LaplacianSmoothPointCloud::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::getGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::PointCloudFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::SmoothingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::getHumanLabel() const
{
  return "Smooth Point Cloud (Laplacian)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid LaplacianSmoothPointCloud::getUuid() const
{
  return QUuid("{1cf52f08-a11a-5870-a38c-ea8958071bd8}");
}

// -----------------------------------------------------------------------------
LaplacianSmoothPointCloud::Pointer LaplacianSmoothPointCloud::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<LaplacianSmoothPointCloud> LaplacianSmoothPointCloud::New()
{
  struct make_shared_enabler : public LaplacianSmoothPointCloud
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::getNameOfClass() const
{
  return QString("LaplacianSmoothPointCloud");
}

// -----------------------------------------------------------------------------
QString LaplacianSmoothPointCloud::ClassName()
{
  return QString("LaplacianSmoothPointCloud");
}

// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::setDataContainerName(const DataArrayPath& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath LaplacianSmoothPointCloud::getDataContainerName() const
{
  return m_DataContainerName;
}

// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::setLambda(float value)
{
  m_Lambda = value;
}

// -----------------------------------------------------------------------------
float LaplacianSmoothPointCloud::getLambda() const
{
  return m_Lambda;
}

// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::setNumIterations(int value)
{
  m_NumIterations = value;
}

// -----------------------------------------------------------------------------
int LaplacianSmoothPointCloud::getNumIterations() const
{
  return m_NumIterations;
}

// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::setUseMask(bool value)
{
  m_UseMask = value;
}

// -----------------------------------------------------------------------------
bool LaplacianSmoothPointCloud::getUseMask() const
{
  return m_UseMask;
}

// -----------------------------------------------------------------------------
void LaplacianSmoothPointCloud::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath LaplacianSmoothPointCloud::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}
