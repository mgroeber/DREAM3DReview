#include "CutPatchesFromImage.h"

#include <array>
#include <assert.h>
#include <cfenv>
#include <random>
#include <chrono>
#include <unordered_map>

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task_scheduler_init.h>
#endif

#include <QtCore/QDateTime>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QJsonArray>
#include <QtCore/QJsonDocument>

#include <Eigen/Dense>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Utilities/SIMPLibEndian.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputPathFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/RectGridGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Math/SIMPLibRandom.h" 
#include "SIMPLib/Utilities/TimeUtilities.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/**
 * @brief The FindThermalHistoriesImpl class implements a threaded algorithm that computes time vs temperature
 * between two arrays
 */
class CutPatchesImpl
{
public:
  CutPatchesImpl(IntVec3Type patchDims, std::vector<size_t> patchIndices, QVector<IDataArray::WeakPointer> selectedWeakPtrVector, ImageGeom::Pointer image, QString directory, CutPatchesFromImage* filter)
  : m_PatchDims(patchDims)
  , m_PatchIndices(patchIndices)
  , m_Image(image)
  , m_SelectedWeakPtrVector(selectedWeakPtrVector)
  , m_Directory(directory)
  , m_Filter(filter)
  {
  }
  virtual ~CutPatchesImpl()
  {
  }

  void generate(size_t start, size_t end) const
  {
    size_t index = 0;
    int64_t numPoints = 0;
    float delX = 0.0f, delY = 0.0f, delZ = 0.0f;
    double pointXPos = 0.0;
    double pointYPos = 0.0;
    double pointZPos = 0.0;
    int32_t counter = 0;
    bool influenced = false;
    float xPos = 0.0f;
    float yPos = 0.0f;
    float zPos = 0.0f;

    SizeVec3Type udims = m_Image->getDimensions();
    size_t m_ImageDims[3] = {static_cast<size_t>(udims[0]), static_cast<size_t>(udims[1]), static_cast<size_t>(udims[2])};

    int64_t progCounter = 0;
    int64_t totalElements = (end - start);
    int64_t progIncrement = static_cast<int64_t>(totalElements / 100);

	  QVariant var;

    QString filtName = "ITKImageWriter";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    DataContainerArray::Pointer dca = DataContainerArray::New();
    DataContainer::Pointer dc = DataContainer::New("tempDC");
    dca->addOrReplaceDataContainer(dc);

    ImageGeom::Pointer iGeom = ImageGeom::CreateGeometry("testGeom");
    iGeom->setDimensions(m_PatchDims[0], m_PatchDims[1], m_PatchDims[2]);
    iGeom->setSpacing(1.0, 1.0, 1.0);
    iGeom->setOrigin(0.0, 0.0, 0.0);
    dc->setGeometry(iGeom);

    IDataArray::Pointer m_SelectedArray;

    int32_t numArrays = m_SelectedWeakPtrVector.size();
    for(size_t iter = 0; iter < numArrays; iter++)
    {
      m_SelectedArray = m_SelectedWeakPtrVector.at(iter).lock();

      QVector<size_t> tDims(3, 0);
      tDims[0] = m_PatchDims[0];
      tDims[1] = m_PatchDims[1];
      tDims[2] = m_PatchDims[2];
      AttributeMatrix::Pointer am = AttributeMatrix::New(tDims, "tempAM", AttributeMatrix::Type::Cell);
      dc->addOrReplaceAttributeMatrix(am);
      QVector<size_t> cDims = m_SelectedArray->getComponentDimensions();
      IDataArray::Pointer patchPtr = m_SelectedArray->createNewArray(tDims[0] * tDims[1] * tDims[2], cDims, "tempAA", true);
      int err = am->addOrReplaceAttributeArray(patchPtr);

      AbstractFilter::Pointer filter;
      QString filename;
      if(nullptr != filterFactory.get())
      {
        filter = filterFactory->create();
        filter->setDataContainerArray(dca);
        filter->setProperty("Plane", 0);
        DataArrayPath tempPath;
        tempPath.update("tempDC", "tempAM", "tempAA");
        var.setValue(tempPath);
        filter->setProperty("ImageArrayPath", var);
        filename = m_Directory + "\\" + m_SelectedArray->getName() + "_Patch";
      }

      for(size_t i = start; i < end; i++)
      {
        if(m_Filter->getCancel())
        {
          return;
        }

        size_t curIndex = m_PatchIndices[i];

        size_t x = curIndex % m_ImageDims[0];
        size_t y = static_cast<size_t>(curIndex / m_ImageDims[0]) % m_ImageDims[1];
        size_t z = static_cast<size_t>(curIndex / (m_ImageDims[0] * m_ImageDims[1]));

        if(x >= static_cast<size_t>(m_PatchDims[0] / 2) && y >= static_cast<size_t>(m_PatchDims[1] / 2) && z >= static_cast<size_t>(m_PatchDims[2] / 2) &&
           x <= (m_ImageDims[0] - static_cast<size_t>(m_PatchDims[0] / 2)) && y <= (m_ImageDims[1] - static_cast<size_t>(m_PatchDims[1] / 2)) &&
           z <= (m_ImageDims[2] - static_cast<size_t>(m_PatchDims[2] / 2)))
        {
          size_t counter = 0;
          size_t planeChunk, rowChunk;
          size_t iStart = (x - static_cast<size_t>(m_PatchDims[0] / 2));
          size_t jStart = (y - static_cast<size_t>(m_PatchDims[1] / 2));
          size_t kStart = (z - static_cast<size_t>(m_PatchDims[2] / 2));
          for(size_t k = kStart; k < (kStart + m_PatchDims[2]); k++)
          {
            planeChunk = k * m_ImageDims[0] * m_ImageDims[1];
            for(size_t j = jStart; j < (jStart + m_PatchDims[1]); j++)
            {
              rowChunk = j * m_ImageDims[0];
              patchPtr->copyFromArray(counter, m_SelectedArray, (planeChunk + rowChunk + iStart), m_PatchDims[0]);
              counter += m_PatchDims[0];
            }
          }
        }
        else
        {
          continue;
        }

        QString cFilename = filename + '-' + QString::number(curIndex) + ".tif";
        var.setValue(cFilename);
        filter->setProperty("FileName", var);
        filter->execute();

        if(progCounter > progIncrement)
        {
          m_Filter->sendThreadSafeProgressMessage(progCounter);
          progCounter = 0;
        }
        progCounter++;
      }
    }
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    generate(r.begin(), r.end());
  }
#endif

private:
  IntVec3Type m_PatchDims;
  std::vector<size_t> m_PatchIndices;
  ImageGeom::Pointer m_Image;
  QString m_Directory;
  QVector<IDataArray::WeakPointer> m_SelectedWeakPtrVector;
  CutPatchesFromImage* m_Filter;
};

// Include the MOC generated file for this class
#include "moc_CutPatchesFromImage.cpp"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CutPatchesFromImage::CutPatchesFromImage()
: m_OutputDirectory("")
, m_SelectedDataArrayPaths(QVector<DataArrayPath>())
, m_LabelsArrayPath("", "", "")
, m_NumPatches(1)
, m_BalancePatches(false)
{
  m_PatchDims[0] = 1;
  m_PatchDims[1] = 1;
  m_PatchDims[2] = 1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CutPatchesFromImage::~CutPatchesFromImage() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CutPatchesFromImage::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  parameters.push_back(SIMPL_NEW_OUTPUT_PATH_FP("Output Directory", OutputDirectory, FilterParameter::Parameter, CutPatchesFromImage));

  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Patchs", NumPatches, FilterParameter::Parameter, CutPatchesFromImage, 0));

  parameters.push_back(SIMPL_NEW_INT_VEC3_FP("Patch Dimensions", PatchDims, FilterParameter::Parameter, CutPatchesFromImage, 0));
  {
    MultiDataArraySelectionFilterParameter::RequirementType req =
        MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Any, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays to Combine", SelectedDataArrayPaths, FilterParameter::RequiredArray, CutPatchesFromImage, req));
  }
  QStringList linkedProps("LabelsArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Balance Patches", BalancePatches, FilterParameter::Parameter, CutPatchesFromImage, linkedProps));
  linkedProps.clear();
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Labels Array", LabelsArrayPath, FilterParameter::RequiredArray, CutPatchesFromImage, req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CutPatchesFromImage::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setOutputDirectory(reader->readString("OutputDirectory", getOutputDirectory()));
  setSelectedDataArrayPaths(reader->readDataArrayPathVector("SelectedDataArrayPaths", getSelectedDataArrayPaths()));
  setLabelsArrayPath(reader->readDataArrayPath("LabelsArrayPath", getLabelsArrayPath()));
  setNumPatches(reader->readValue("NumPatches", getNumPatches()));
  setBalancePatches(reader->readValue("BalancePatches", getBalancePatches()));
  setPatchDims(reader->readIntVec3("PatchDims", getPatchDims()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CutPatchesFromImage::initialize()
{
  m_ProgressCounter = 0;
  m_TotalElements = 0;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CutPatchesFromImage::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();
  DataArrayPath tempPath;
  QVector<IDataArray::Pointer> dataArrays;

  m_SelectedWeakPtrVector.clear();

  QVector<DataArrayPath> paths = getSelectedDataArrayPaths();

  if(!DataArrayPath::ValidateVector(paths))
  {
    QString ss = QObject::tr("There are Attribute Arrays selected that are not contained in the same Attribute Matrix. All selected Attribute Arrays must belong to the same Attribute Matrix");
    setErrorCondition(-11002, ss);
  }

  for(int32_t i = 0; i < paths.count(); i++)
  {
    DataArrayPath path = paths.at(i);
    IDataArray::WeakPointer ptr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, path);
    if(nullptr != ptr.lock())
    {
      m_SelectedWeakPtrVector.push_back(ptr);
    }
  }

  if(getPatchDims()[0] <= 0 || getPatchDims()[1] <= 0 || getPatchDims()[2] <= 0)
  {
    QString ss = QObject::tr("All patch dimensions must be positive.\n "
                             "Current kernel dimensions:\n x = %1\n y = %2\n z = %3\n")
                     .arg(getPatchDims()[0])
                     .arg(getPatchDims()[1])
                     .arg(getPatchDims()[2]);
    setErrorCondition(-11000, ss);
  }

  if (m_BalancePatches)
  {
    QVector<size_t> cDims(1, 1);
    m_LabelsArrayPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getLabelsArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_LabelsArrayPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_LabelsArray = m_LabelsArrayPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
  }
  if(getErrorCode() >= 0)
  {
    dataArrays.push_back(m_LabelsArrayPtr.lock());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CutPatchesFromImage::preflight()
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
void CutPatchesFromImage::sendThreadSafeProgressMessage(int64_t counter)
{
  m_Mutex.lock();

  m_ProgressCounter += counter;
  int64_t progressInt = static_cast<int64_t>((static_cast<float>(m_ProgressCounter) / m_TotalElements) * 100.0f);

  int64_t progIncrement = m_TotalElements / 100;
  int64_t prog = 1;

  if(m_ProgressCounter > prog && m_LastProgressInt != progressInt)
  {
    QString ss = QObject::tr("Evaluating Thermal Histories || %1% Completed").arg(progressInt);
    notifyStatusMessage(ss);
    prog += progIncrement;
  }

  m_LastProgressInt = progressInt;

  m_Mutex.unlock();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CutPatchesFromImage::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  tbb::task_scheduler_init init;
  bool doParallel = true;
#endif

  ImageGeom::Pointer image = getDataContainerArray()->getDataContainer(getSelectedDataArrayPaths().at(0).getDataContainerName())->getGeometryAs<ImageGeom>();
  int64_t numPossiblePatches = image->getNumberOfElements();
  
  float fraction[2];
  if(!m_BalancePatches)
  {
    fraction[0] = float(m_NumPatches) / float(numPossiblePatches);
    fraction[1] = 0.0f;
  }
  else
  {
    float count1 = 0.0f;
    float count2 = 0.0f;
    for (size_t i = 0; i < numPossiblePatches; i++)
    {
      if(m_LabelsArray[i])
      {
        count1 += 1.0f;
      }
      else
      {
        count2 += 1.0f;
      }
    }
    fraction[0] = float(m_NumPatches) / (2 * count1);
    fraction[1] = float(m_NumPatches) / (2 * count2);
  }

  notifyStatusMessage("Selecting Patches...");

  std::vector<size_t> patchIndices;
  size_t m_Seed = std::chrono::steady_clock::now().time_since_epoch().count();
  SIMPL_RANDOMNG_NEW_SEEDED(m_Seed);
  for(size_t i = 0; i < numPossiblePatches; i++)
  {
    float randNum = rg.genrand_res53();
    if(m_BalancePatches)
    {
      if(m_LabelsArray[i] && randNum > fraction[0])
      {
        continue;
      }
      else if(!m_LabelsArray[i] && randNum > fraction[1])
      {
        continue;
      }
    }
    else
    {
      if(randNum > fraction[0])
      {
        continue;
      }
    }
    patchIndices.push_back(i);
  }

  notifyStatusMessage("Cutting Patches...");
  size_t numPatches = patchIndices.size();

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  if(doParallel == true)
  {

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numPatches), CutPatchesImpl(m_PatchDims, patchIndices, m_SelectedWeakPtrVector, image, m_OutputDirectory, this),
                      tbb::simple_partitioner());
  }
  else
#endif
  {
    CutPatchesImpl serial(m_PatchDims, patchIndices, m_SelectedWeakPtrVector, image, m_OutputDirectory, this);
    serial.generate(0, numPatches);
  }
 
  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer CutPatchesFromImage::newFilterInstance(bool copyFilterParameters) const
{
  CutPatchesFromImage::Pointer filter = CutPatchesFromImage::New();
  if(true == copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CutPatchesFromImage::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CutPatchesFromImage::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CutPatchesFromImage::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CutPatchesFromImage::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CutPatchesFromImage::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CutPatchesFromImage::getHumanLabel() const
{
  return "Cut Patches From Image";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid CutPatchesFromImage::getUuid()
{
  return QUuid("");
}