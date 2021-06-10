/* ============================================================================
 * Copyright (c) 2009-2016 BlueQuartz Software, LLC
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
 * contributors may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The code contained herein was partially funded by the followig contracts:
 *    United States Air Force Prime Contract FA8650-07-D-5800
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *    United States Prime Contract Navy N00173-07-C-2068
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include "SingleDirectionMeanFilter.h"

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#endif

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/**
 * @brief The DirectionalDifferencesImpl class implements a templated threaded algorithm for
 * determining the voxel to voxel differences in a given direction of a 3d image.
 */
template <typename T>
class DirectionalDifferencesImpl
{

public:
  DirectionalDifferencesImpl(T* data, float* filtArray, size_t length, size_t width, size_t depth, int64_t dStride, int64_t lStride, int64_t wStride)
  : m_Data(data)
  , m_FilteredArray(filtArray)
  , m_Length(length)
  , m_Width(width)
  , m_Depth(depth)
  , m_DStride(dStride)
  , m_LStride(lStride)
  , m_WStride(wStride)
  {
  }
  virtual ~DirectionalDifferencesImpl() = default;

  void convert(size_t start, size_t end) const
  {
    int64_t m_DepthStep = fabs(m_DStride);
    for(size_t iter = start; iter < end; iter++)
    {
      for(size_t i = 0; i < m_Length; i++)
      {
        for(size_t j = 0; j < m_Width; j++)
        {
          int64_t kernelIter = 0;
          int64_t curDepth = iter;
          size_t point = (iter * m_DepthStep) + (i * m_LStride) + (j * m_WStride);
          while(kernelIter >= 20 && curDepth >= 0 && curDepth < m_Depth)
          {
            m_FilteredArray[point] += m_Data[point + (kernelIter * m_DStride)];
            kernelIter++;
            if(m_DStride > 0)
            {
              curDepth++;
            }
            else if(m_DStride < 0)
            {
              curDepth--;
            }
          }
          m_FilteredArray[point] /= float(kernelIter);
        }
      }
    }
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    convert(r.begin(), r.end());
  }
#endif
private:
  T* m_Data;
  float* m_FilteredArray;
  size_t m_Length;
  size_t m_Width;
  size_t m_Depth;
  int64_t m_DStride;
  int64_t m_LStride;
  int64_t m_WStride;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SingleDirectionMeanFilter::SingleDirectionMeanFilter() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SingleDirectionMeanFilter::~SingleDirectionMeanFilter() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Direction of Interest");
    parameter->setPropertyName("Direction");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(SingleDirectionMeanFilter, this, Direction));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(SingleDirectionMeanFilter, this, Direction));

    std::vector<QString> choices;
    choices.push_back("+X");
    choices.push_back("-X");
    choices.push_back("+Y");
    choices.push_back("-Y");
    choices.push_back("+Z");
    choices.push_back("-Z");
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    std::vector<QString> daTypes;
    daTypes.push_back(SIMPL::TypeNames::Int8);
    daTypes.push_back(SIMPL::TypeNames::Int16);
    daTypes.push_back(SIMPL::TypeNames::Int32);
    daTypes.push_back(SIMPL::TypeNames::Int64);
    daTypes.push_back(SIMPL::TypeNames::UInt8);
    daTypes.push_back(SIMPL::TypeNames::UInt16);
    daTypes.push_back(SIMPL::TypeNames::UInt32);
    daTypes.push_back(SIMPL::TypeNames::UInt64);
    daTypes.push_back(SIMPL::TypeNames::Float);
    daTypes.push_back(SIMPL::TypeNames::Double);
    req.daTypes = daTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Filter", SelectedArrayPath, FilterParameter::Category::RequiredArray, SingleDirectionMeanFilter, req));
  }
  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Filtered Array", FilteredArrayName, SelectedArrayPath, SelectedArrayPath, FilterParameter::Category::CreatedArray, SingleDirectionMeanFilter));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setFilteredArrayName(reader->readString("FilteredArrayName", getFilteredArrayName()));
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setDirection(reader->readValue("Direction", getDirection()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  DataArrayPath tempPath;

  m_InDataPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, getSelectedArrayPath());
  if(nullptr != m_InDataPtr.lock())
  {
    if(TemplateHelpers::CanDynamicCast<BoolArrayType>()(m_InDataPtr.lock()))
    {
      QString ss = QObject::tr("Selected array cannot be of type bool.  The path is %1").arg(getSelectedArrayPath().serialize());
      setErrorCondition(-11001, ss);
    }
  }

  std::vector<size_t> cDims(1, 1);
  tempPath.update(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getFilteredArrayName());
  m_FilteredArrayPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_FilteredArrayPtr.lock())
  {
    m_FilteredArray = m_FilteredArrayPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  ImageGeom::Pointer image = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom>(this, getSelectedArrayPath().getDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }

  if(image->getXPoints() <= 1 || image->getYPoints() <= 1 || image->getZPoints() <= 1)
  {
    setErrorCondition(-999, "The Image Geometry is not 3D and cannot be run through this Filter");
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedArrayPath().getDataContainerName());

  SizeVec3Type geoDims = m->getGeometryAs<ImageGeom>()->getDimensions();

  int64_t dStride = 0, lStride = 0, wStride;
  size_t count = geoDims[0] * geoDims[1] * geoDims[2];
  size_t length = 0;
  size_t width = 0;
  size_t depth = 0;
  if(m_Direction == 0)
  {
    dStride = 1;
    lStride = geoDims[0];
    wStride = geoDims[0] * geoDims[1];

    depth = geoDims[0];
    length = geoDims[1];
    width = geoDims[2];
  }
  if(m_Direction == 1)
  {
    dStride = -1;
    lStride = geoDims[0];
    wStride = geoDims[0] * geoDims[1];

    depth = geoDims[0];
    length = geoDims[1];
    width = geoDims[2];
  }
  if(m_Direction == 2)
  {
    dStride = geoDims[0];
    lStride = 1;
    wStride = geoDims[0] * geoDims[1];

    depth = geoDims[1];
    length = geoDims[0];
    width = geoDims[2];
  }
  if(m_Direction == 3)
  {
    dStride = -geoDims[0];
    lStride = 1;
    wStride = geoDims[0] * geoDims[1];

    depth = geoDims[1];
    length = geoDims[0];
    width = geoDims[2];
  }
  if(m_Direction == 4)
  {
    dStride = geoDims[0] * geoDims[1];
    lStride = 1;
    wStride = geoDims[0];

    depth = geoDims[2];
    length = geoDims[0];
    width = geoDims[1];
  }
  if(m_Direction == 5)
  {
    dStride = -geoDims[0] * geoDims[1];
    lStride = 1;
    wStride = geoDims[0];

    depth = geoDims[2];
    length = geoDims[0];
    width = geoDims[1];
  }

  IDataArray::Pointer iCellArray = m_InDataPtr.lock();

  if(TemplateHelpers::CanDynamicCast<Int8ArrayType>()(m_InDataPtr.lock()))
  {
    Int8ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int8ArrayType>(m_InDataPtr.lock());
    int8_t* cPtr = cellArray->getPointer(0);
    std::vector<int8_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<int8_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());

#else
    DirectionalDifferencesImpl<int8_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt8ArrayType>()(m_InDataPtr.lock()))
  {
    UInt8ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt8ArrayType>(m_InDataPtr.lock());
    uint8_t* cPtr = cellArray->getPointer(0);
    std::vector<uint8_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<uint8_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());

#else
    DirectionalDifferencesImpl<uint8_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int16ArrayType>()(m_InDataPtr.lock()))
  {
    Int16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int16ArrayType>(m_InDataPtr.lock());
    int16_t* cPtr = cellArray->getPointer(0);
    std::vector<int16_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<int16_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());
#else
    DirectionalDifferencesImpl<int16_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt16ArrayType>()(m_InDataPtr.lock()))
  {
    UInt16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt16ArrayType>(m_InDataPtr.lock());
    uint16_t* cPtr = cellArray->getPointer(0);
    std::vector<uint16_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<uint16_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());
#else
    DirectionalDifferencesImpl<uint16_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int32ArrayType>()(m_InDataPtr.lock()))
  {
    Int32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int32ArrayType>(m_InDataPtr.lock());
    int32_t* cPtr = cellArray->getPointer(0);
    std::vector<int32_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<int32_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());

#else
    DirectionalDifferencesImpl<int32_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt32ArrayType>()(m_InDataPtr.lock()))
  {
    UInt32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt32ArrayType>(m_InDataPtr.lock());
    uint32_t* cPtr = cellArray->getPointer(0);
    std::vector<uint32_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<uint32_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());

#else
    DirectionalDifferencesImpl<uint32_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int64ArrayType>()(m_InDataPtr.lock()))
  {
    Int64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int64ArrayType>(m_InDataPtr.lock());
    int64_t* cPtr = cellArray->getPointer(0);
    std::vector<int64_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<int64_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());

#else
    DirectionalDifferencesImpl<int64_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt64ArrayType>()(m_InDataPtr.lock()))
  {
    UInt64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt64ArrayType>(m_InDataPtr.lock());
    uint64_t* cPtr = cellArray->getPointer(0);
    std::vector<uint64_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<uint64_t>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());
#else
    DirectionalDifferencesImpl<uint64_t> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<FloatArrayType>()(m_InDataPtr.lock()))
  {
    FloatArrayType::Pointer cellArray = std::dynamic_pointer_cast<FloatArrayType>(m_InDataPtr.lock());
    float* cPtr = cellArray->getPointer(0);
    std::vector<float> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<float>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());
#else
    DirectionalDifferencesImpl<float> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<DoubleArrayType>()(m_InDataPtr.lock()))
  {
    DoubleArrayType::Pointer cellArray = std::dynamic_pointer_cast<DoubleArrayType>(m_InDataPtr.lock());
    double* cPtr = cellArray->getPointer(0);
    std::vector<double> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), DirectionalDifferencesImpl<double>(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride), tbb::auto_partitioner());
#else
    DirectionalDifferencesImpl<double> serial(cPtr, m_FilteredArray, length, width, depth, dStride, lStride, wStride);
    serial.convert(0, depth);
#endif
  }
  else
  {
    QString ss = QObject::tr("Selected array is of unsupported type. The type is %1").arg(m_InDataPtr.lock()->getTypeAsString());
    setErrorCondition(-11001, ss);
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer SingleDirectionMeanFilter::newFilterInstance(bool copyFilterParameters) const
{
  SingleDirectionMeanFilter::Pointer filter = SingleDirectionMeanFilter::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getGroupName() const
{
  return SIMPL::FilterGroups::ProcessingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid SingleDirectionMeanFilter::getUuid() const
{
  return QUuid("{577dfdf6-02f8-5284-b45b-e31f5412b166}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::ImageFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getHumanLabel() const
{
  return "Single Direction Mean Filter";
}

// -----------------------------------------------------------------------------
SingleDirectionMeanFilter::Pointer SingleDirectionMeanFilter::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<SingleDirectionMeanFilter> SingleDirectionMeanFilter::New()
{
  struct make_shared_enabler : public SingleDirectionMeanFilter
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getNameOfClass() const
{
  return QString("SingleDirectionMeanFilter");
}

// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::ClassName()
{
  return QString("SingleDirectionMeanFilter");
}

// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath SingleDirectionMeanFilter::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::setDirection(unsigned int value)
{
  m_Direction = value;
}

// -----------------------------------------------------------------------------
unsigned int SingleDirectionMeanFilter::getDirection() const
{
  return m_Direction;
}

// -----------------------------------------------------------------------------
void SingleDirectionMeanFilter::setFilteredArrayName(const QString& value)
{
  m_FilteredArrayName = value;
}

// -----------------------------------------------------------------------------
QString SingleDirectionMeanFilter::getFilteredArrayName() const
{
  return m_FilteredArrayName;
}