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
#include "MapPointCloudToRegularGrid.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MapPointCloudToRegularGrid::MapPointCloudToRegularGrid()
: m_DataContainerName("")
, m_ImageDataContainerName("ImageDataContainer")
, m_ImageDataContainerPath("")
, m_VoxelIndicesArrayPath("", "", "VoxelIndices")
, m_UseMask(false)
, m_CreateDataContainer(0)
, m_MaskArrayPath("", "", "")
{
  m_GridDimensions[0] = 10;
  m_GridDimensions[1] = 10;
  m_GridDimensions[2] = 10;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MapPointCloudToRegularGrid::~MapPointCloudToRegularGrid() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MapPointCloudToRegularGrid::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Create Data Container");
    parameter->setPropertyName("CreateDataContainer");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(MapPointCloudToRegularGrid, this, CreateDataContainer));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(MapPointCloudToRegularGrid, this, CreateDataContainer));

    QVector<QString> choices;
    choices.push_back("Create New Data Image Data Container");
    choices.push_back("Use Exsting Data Image Data Container");
    parameter->setChoices(choices);
    QStringList linkedProps2;
    linkedProps2 << "GridDimensions"
                 << "ImageDataContainerName"
                 << "ImageDataContainerPath";
    parameter->setLinkedProperties(linkedProps2);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_INT_VEC3_FP("Grid Dimensions", GridDimensions, FilterParameter::Parameter, MapPointCloudToRegularGrid, 0));
  parameters.push_back(SIMPL_NEW_STRING_FP("Image Data Container", ImageDataContainerName, FilterParameter::CreatedArray, MapPointCloudToRegularGrid, 0));
  {
    DataContainerSelectionFilterParameter::RequirementType req;
    IGeometry::Types reqGeom = {IGeometry::Type::Image};
    req.dcGeometryTypes = reqGeom;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("ImageDataContainerPath", ImageDataContainerPath, FilterParameter::RequiredArray, MapPointCloudToRegularGrid, req, 1));
  }
  {
    DataContainerSelectionFilterParameter::RequirementType req;
    IGeometry::Types reqGeom = {IGeometry::Type::Vertex};
    req.dcGeometryTypes = reqGeom;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Data Container to Map", DataContainerName, FilterParameter::RequiredArray, MapPointCloudToRegularGrid, req));
  }
  QStringList linkedProps("MaskArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, MapPointCloudToRegularGrid, linkedProps));
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, MapPointCloudToRegularGrid, req));
  }
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::CreatedArray));
  {
    DataArrayCreationFilterParameter::RequirementType req = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
    parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Voxel Indices", VoxelIndicesArrayPath, FilterParameter::CreatedArray, MapPointCloudToRegularGrid, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MapPointCloudToRegularGrid::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setDataContainerName(reader->readDataArrayPath("DataContainerName", getDataContainerName()));
  setImageDataContainerPath(reader->readDataArrayPath("ImageDataContainerPath", getImageDataContainerPath()));
  setImageDataContainerName(reader->readString("ImageDataContainerName", getImageDataContainerName()));
  setVoxelIndicesArrayPath(reader->readDataArrayPath("VoxelIndicesArrayPath", getVoxelIndicesArrayPath()));
  setGridDimensions(reader->readIntVec3("GridDimensions", getGridDimensions()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setCreateDataContainer(reader->readValue("CreateDataContainer", getCreateDataContainer()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MapPointCloudToRegularGrid::initialize()
{
  m_MeshMinExtents.clear();
  m_MeshMaxExtents.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MapPointCloudToRegularGrid::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

  QVector<IDataArray::Pointer> dataArrays;

  VertexGeom::Pointer vertex = getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom, AbstractFilter>(this, getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  dataArrays.push_back(vertex->getVertices());

  if(m_CreateDataContainer == 0)
  {
    if(getGridDimensions()[0] <= 0 || getGridDimensions()[1] <= 0 || getGridDimensions()[2] <= 0)
    {
      QString ss = QObject::tr("All grid dimensions must be positive.\n "
                               "Current grid dimensions:\n x = %1\n y = %2\n z = %3\n")
                       .arg(getGridDimensions()[0])
                       .arg(getGridDimensions()[1])
                       .arg(getGridDimensions()[2]);
      setErrorCondition(-11000, ss);
    }

    if(getErrorCode() < 0)
    {
      return;
    }

    ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
    size_t dims[3] = {static_cast<size_t>(getGridDimensions()[0]), static_cast<size_t>(getGridDimensions()[1]), static_cast<size_t>(getGridDimensions()[2])};
    image->setDimensions(dims);

    DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getImageDataContainerName());

    if(getErrorCode() < 0)
    {
      return;
    }

    m->setGeometry(image);
  }

  if(m_CreateDataContainer == 1)
  {
    ImageGeom::Pointer vertex = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getImageDataContainerPath());
    if(getErrorCode() < 0)
    {
      return;
    }
  }

  std::vector<size_t> cDims(1, 1);

  m_VoxelIndicesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<size_t>, AbstractFilter, size_t>(
      this, getVoxelIndicesArrayPath(), 0, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_VoxelIndicesPtr.lock().get())    /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_VoxelIndices = m_VoxelIndicesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrays.push_back(m_VoxelIndicesPtr.lock());
  }

  if(getUseMask() == true)
  {
    m_MaskPtr =
        getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getMaskArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_MaskPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() >= 0)
    {
      dataArrays.push_back(m_MaskPtr.lock());
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrays);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MapPointCloudToRegularGrid::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MapPointCloudToRegularGrid::createRegularGrid()
{
  QString ss = QObject::tr("Creating Regular Grid");
  notifyStatusMessage(ss);

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());
  DataContainer::Pointer interpolatedDC = getDataContainerArray()->getDataContainer(getImageDataContainerName());
  VertexGeom::Pointer pointCloud = m->getGeometryAs<VertexGeom>();
  ImageGeom::Pointer image = interpolatedDC->getGeometryAs<ImageGeom>();

  int64_t numVerts = pointCloud->getNumberOfVertices();
  float* vertex = pointCloud->getVertexPointer(0);

  // Find the largest/smallest (x,y,z) dimensions of the incoming data to be used to define the maximum dimensions for the regular grid
  for(size_t i = 0; i < 3; i++)
  {
    m_MeshMaxExtents.push_back(std::numeric_limits<float>::lowest());
    m_MeshMinExtents.push_back(std::numeric_limits<float>::max());
  }

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(!m_UseMask || (m_UseMask && m_Mask[i]))
    {
      if(vertex[3 * i] > m_MeshMaxExtents[0])
      {
        m_MeshMaxExtents[0] = vertex[3 * i];
      }
      if(vertex[3 * i + 1] > m_MeshMaxExtents[1])
      {
        m_MeshMaxExtents[1] = vertex[3 * i + 1];
      }
      if(vertex[3 * i + 2] > m_MeshMaxExtents[2])
      {
        m_MeshMaxExtents[2] = vertex[3 * i + 2];
      }
      if(vertex[3 * i] < m_MeshMinExtents[0])
      {
        m_MeshMinExtents[0] = vertex[3 * i];
      }
      if(vertex[3 * i + 1] < m_MeshMinExtents[1])
      {
        m_MeshMinExtents[1] = vertex[3 * i + 1];
      }
      if(vertex[3 * i + 2] < m_MeshMinExtents[2])
      {
        m_MeshMinExtents[2] = vertex[3 * i + 2];
      }
    }
  }

  SizeVec3Type iDims = image->getDimensions();

  QVector<float> iRes(3, 0.0f);
  QVector<float> iOrigin(3, 0.0f);

  if(iDims[0] > 1)
  {
    iRes[0] = (m_MeshMaxExtents[0] - m_MeshMinExtents[0]) / (static_cast<float>(iDims[0]));
    if(iRes[0] == 0.0)
    {
      iRes[0] = 1.0f;
      iDims[0] = 1;
    }
    else
    {
      iDims[0] += 1;
    }
    iOrigin[0] = m_MeshMinExtents[0] - (iRes[0] / 2.0f);
  }
  else
  {
    iRes[0] = 1.25f * (m_MeshMaxExtents[0] - m_MeshMinExtents[0]) / (static_cast<float>(iDims[0]));
    if(iRes[0] == 0.0)
    {
      iRes[0] = 1.0f;
      iOrigin[0] = -0.5f;
    }
    else
    {
      iOrigin[0] = m_MeshMinExtents[0] - (iRes[0] * 0.1f);
    }
  }

  if(iDims[1] > 1)
  {
    iRes[1] = (m_MeshMaxExtents[1] - m_MeshMinExtents[1]) / (static_cast<float>(iDims[1]));
    if(iRes[1] == 0.0)
    {
      iRes[1] = 1.0f;
      iDims[1] = 1;
    }
    else
    {
      iDims[1] += 1;
    }
    iOrigin[1] = m_MeshMinExtents[1] - (iRes[1] / 2.0f);
  }
  else
  {
    iRes[1] = 1.25f * (m_MeshMaxExtents[1] - m_MeshMinExtents[1]) / (static_cast<float>(iDims[1]));
    if(iRes[1] == 0.0)
    {
      iRes[1] = 1.0f;
      iOrigin[1] = -0.5f;
    }
    else
    {
      iOrigin[1] = m_MeshMinExtents[1] - (iRes[1] * 0.1f);
    }
  }

  if(iDims[2] > 1)
  {
    iRes[2] = (m_MeshMaxExtents[2] - m_MeshMinExtents[2]) / (static_cast<float>(iDims[2]));
    if(iRes[2] == 0.0)
    {
      iRes[2] = 1.0f;
      iDims[2] = 1;
    }
    else
    {
      iDims[2] += 1;
    }
    iOrigin[2] = m_MeshMinExtents[2] - (iRes[2] / 2.0f);
  }
  else
  {
    iRes[2] = 1.25f * (m_MeshMaxExtents[2] - m_MeshMinExtents[2]) / (static_cast<float>(iDims[2]));
    if(iRes[2] == 0.0)
    {
      iRes[2] = 1.0f;
      iOrigin[2] = -0.5f;
    }
    else
    {
      iOrigin[2] = m_MeshMinExtents[1] - (iRes[2] * 0.1f);
    }
  }

  image->setDimensions(iDims);
  image->setSpacing(iRes[0], iRes[1], iRes[2]);
  image->setOrigin(iOrigin[0], iOrigin[1], iOrigin[2]);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MapPointCloudToRegularGrid::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  ImageGeom::Pointer image;
  if(m_CreateDataContainer == 0)
  {
    // Create the regular grid
    createRegularGrid();
    image = getDataContainerArray()->getDataContainer(getImageDataContainerName())->getGeometryAs<ImageGeom>();
  }
  else if(m_CreateDataContainer == 1)
  {
    image = getDataContainerArray()->getDataContainer(getImageDataContainerPath())->getGeometryAs<ImageGeom>();
  }

  VertexGeom::Pointer vertices = getDataContainerArray()->getDataContainer(getDataContainerName())->getGeometryAs<VertexGeom>();

  int64_t numVerts = vertices->getNumberOfVertices();
  SizeVec3Type dims = image->getDimensions();
  FloatVec3Type res = image->getSpacing();
  FloatVec3Type origin = image->getOrigin();
  float coords[3] = {0.0f, 0.0f, 0.0f};
  size_t idxs[3] = {0, 0, 0};
  int64_t progIncrement = numVerts / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  int64_t counter = 0;

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(!m_UseMask || (m_UseMask && m_Mask[i]))
    {
      vertices->getCoords(i, coords);

      for(size_t j = 0; j < 3; j++)
      {
        if((coords[j] - origin[j]) < 0)
        {
          QString ss = QObject::tr("Found negative value for index computation of vertex %1, which may result in unsigned underflow").arg(i);
          setWarningCondition(-1000, ss);
        }
        idxs[j] = int64_t(floor((coords[j] - origin[j]) / res[j]));
      }

      for(size_t j = 0; j < 3; j++)
      {
        if(idxs[j] >= dims[j])
        {
          idxs[j] = (dims[j] - 1);
        }
      }

      m_VoxelIndices[i] = (idxs[2] * dims[1] * dims[0]) + (idxs[1] * dims[0]) + idxs[0];

      if(counter > prog)
      {
        progressInt = static_cast<int64_t>((static_cast<float>(counter) / numVerts) * 100.0f);
        QString ss = QObject::tr("Computing Point Cloud Voxel Indices || %1% Completed").arg(progressInt);
        notifyStatusMessage(ss);
        prog = prog + progIncrement;
      }
      counter++;
    }
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer MapPointCloudToRegularGrid::newFilterInstance(bool copyFilterParameters) const
{
  MapPointCloudToRegularGrid::Pointer filter = MapPointCloudToRegularGrid::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString MapPointCloudToRegularGrid::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString MapPointCloudToRegularGrid::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString MapPointCloudToRegularGrid::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString MapPointCloudToRegularGrid::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString MapPointCloudToRegularGrid::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::MappingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString MapPointCloudToRegularGrid::getHumanLabel() const
{
  return "Map Point Cloud to Regular Grid";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid MapPointCloudToRegularGrid::getUuid()
{
  return QUuid("{9fe34deb-99e1-5f3a-a9cc-e90c655b47ee}");
}
