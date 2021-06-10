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
#include "FindArrayDifferencesAlongDirection.h"

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
  DirectionalDifferencesImpl(T* data, float* min, size_t length, size_t width, size_t dStride, size_t lStride, size_t wStride, std::vector<T> layerAverages)
  : m_Data(data)
  , m_Min(min)
  , m_Length(length)
  , m_Width(width)
  , m_DStride(dStride)
  , m_LStride(lStride)
  , m_WStride(wStride)
  , m_LayerAverages(layerAverages)
  {
  }
  virtual ~DirectionalDifferencesImpl() = default;

  void convert(size_t start, size_t end) const
  {
    for(size_t iter = start; iter < end; iter++)
    {
      for(size_t i = 0; i < m_Length; i++)
      {
        for(size_t j = 0; j < m_Width; j++)
        {
          size_t point = (iter * m_DStride) + (i * m_LStride) + (j * m_WStride);
          m_Min[3 * point] = (m_Data[point + m_DStride] - m_Data[point]) - (m_LayerAverages[iter + 1] - m_LayerAverages[iter]);
          m_Min[3 * point + 1] = (m_Data[point + 2 * m_DStride] - m_Data[point]) - (m_LayerAverages[iter + 2] - m_LayerAverages[iter]);
          m_Min[3 * point + 2] = (m_Data[point + 3 * m_DStride] - m_Data[point]) - (m_LayerAverages[iter + 3] - m_LayerAverages[iter]);
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
  float* m_Min;
  size_t m_Length;
  size_t m_Width;
  size_t m_DStride;
  size_t m_LStride;
  size_t m_WStride;
  std::vector<T> m_LayerAverages;
};

/**
 * @brief The LayerAveragesImpl class implements a templated threaded algorithm for
 * determining the average values of layers in a 3D image.
 */
template <typename T>
class LayerAveragesImpl
{

public:
  LayerAveragesImpl(T* data, size_t length, size_t width, size_t dStride, size_t lStride, size_t wStride, std::vector<T>& layerAverages)
  : m_Data(data)
  , m_Length(length)
  , m_Width(width)
  , m_DStride(dStride)
  , m_LStride(lStride)
  , m_WStride(wStride)
  , m_LayerAverages(layerAverages)
  {
  }
  virtual ~LayerAveragesImpl() = default;

  void convert(size_t start, size_t end) const
  {
    for(size_t iter = start; iter < end; iter++)
    {
      size_t count = 0;
      for(size_t i = 0; i < m_Length; i ++)
      {
        for(size_t j = 0; j < m_Width; j++)
        {
          size_t point = (iter * m_DStride) + (i * m_LStride) + (j * m_WStride);
          if(m_Data[point] != 0)
          {
            m_LayerAverages[iter] += m_Data[point];
            count++;
          }
        }
      }
      m_LayerAverages[iter] /= count;
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
  size_t m_Length;
  size_t m_Width;
  size_t m_DStride;
  size_t m_LStride;
  size_t m_WStride;
  std::vector<T>& m_LayerAverages;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindArrayDifferencesAlongDirection::FindArrayDifferencesAlongDirection() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindArrayDifferencesAlongDirection::~FindArrayDifferencesAlongDirection() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayDifferencesAlongDirection::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Direction of Interest");
    parameter->setPropertyName("Direction");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(FindArrayDifferencesAlongDirection, this, Direction));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(FindArrayDifferencesAlongDirection, this, Direction));

    std::vector<QString> choices;
    choices.push_back("X");
    choices.push_back("Y");
    choices.push_back("Z");
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
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Quantify", SelectedArrayPath, FilterParameter::Category::RequiredArray, FindArrayDifferencesAlongDirection, req));
  }
  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Projected Image Min", ProjectedImageMinArrayName, SelectedArrayPath, SelectedArrayPath, FilterParameter::Category::CreatedArray, FindArrayDifferencesAlongDirection));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
void FindArrayDifferencesAlongDirection::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setProjectedImageMinArrayName(reader->readString("ProjectedImageMinArrayName", getProjectedImageMinArrayName()));
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setDirection(reader->readValue("Direction", getDirection()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayDifferencesAlongDirection::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayDifferencesAlongDirection::dataCheck()
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

  std::vector<size_t> cDims(1, 3);
  tempPath.update(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getProjectedImageMinArrayName());
  m_ProjectedImageMinPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_ProjectedImageMinPtr.lock())
  {
    m_ProjectedImageMin = m_ProjectedImageMinPtr.lock()->getPointer(0);
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
void FindArrayDifferencesAlongDirection::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedArrayPath().getDataContainerName());

  SizeVec3Type geoDims = m->getGeometryAs<ImageGeom>()->getDimensions();

  size_t dStride = 0, lStride = 0, wStride;
  size_t count = geoDims[0]*geoDims[1]*geoDims[2];
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
    dStride = geoDims[0];
    lStride = 1;
    wStride = geoDims[0] * geoDims[1];
    
    depth = geoDims[1];
    length = geoDims[0];
    width = geoDims[2];
  }
  if(m_Direction == 2)
  {
    dStride = geoDims[0] * geoDims[1];
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
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<int8_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth-1), DirectionalDifferencesImpl<int8_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());

#else
    LayerAveragesImpl<int8_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<int8_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt8ArrayType>()(m_InDataPtr.lock()))
  {
    UInt8ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt8ArrayType>(m_InDataPtr.lock());
    uint8_t* cPtr = cellArray->getPointer(0);
    std::vector<uint8_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<uint8_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<uint8_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());

#else
    LayerAveragesImpl<uint8_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<uint8_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int16ArrayType>()(m_InDataPtr.lock()))
  {
    Int16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int16ArrayType>(m_InDataPtr.lock());
    int16_t* cPtr = cellArray->getPointer(0);
    std::vector<int16_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<int16_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<int16_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());
#else
    LayerAveragesImpl<int16_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<int16_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt16ArrayType>()(m_InDataPtr.lock()))
  {
    UInt16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt16ArrayType>(m_InDataPtr.lock());
    uint16_t* cPtr = cellArray->getPointer(0);
    std::vector<uint16_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<uint16_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<uint16_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());
#else
    LayerAveragesImpl<uint16_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<uint16_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int32ArrayType>()(m_InDataPtr.lock()))
  {
    Int32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int32ArrayType>(m_InDataPtr.lock());
    int32_t* cPtr = cellArray->getPointer(0);
    std::vector<int32_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<int32_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<int32_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());

#else
    LayerAveragesImpl<int32_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<int32_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt32ArrayType>()(m_InDataPtr.lock()))
  {
    UInt32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt32ArrayType>(m_InDataPtr.lock());
    uint32_t* cPtr = cellArray->getPointer(0);
    std::vector<uint32_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<uint32_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<uint32_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());

#else
    LayerAveragesImpl<uint32_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<uint32_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int64ArrayType>()(m_InDataPtr.lock()))
  {
    Int64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int64ArrayType>(m_InDataPtr.lock());
    int64_t* cPtr = cellArray->getPointer(0);
    std::vector<int64_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<int64_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<int64_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());

#else
    LayerAveragesImpl<int64_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<int64_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt64ArrayType>()(m_InDataPtr.lock()))
  {
    UInt64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt64ArrayType>(m_InDataPtr.lock());
    uint64_t* cPtr = cellArray->getPointer(0);
    std::vector<uint64_t> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<uint64_t>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<uint64_t>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());
#else
    LayerAveragesImpl<uint64_t> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<uint64_t> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<FloatArrayType>()(m_InDataPtr.lock()))
  {
    FloatArrayType::Pointer cellArray = std::dynamic_pointer_cast<FloatArrayType>(m_InDataPtr.lock());
    float* cPtr = cellArray->getPointer(0);
    std::vector<float> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<float>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<float>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());
#else
    LayerAveragesImpl<float> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<float> serial(cPtr, m_ProjectedImageMin, dStride, count);
    serial.convert(0, depth - 3);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<DoubleArrayType>()(m_InDataPtr.lock()))
  {
    DoubleArrayType::Pointer cellArray = std::dynamic_pointer_cast<DoubleArrayType>(m_InDataPtr.lock());
    double* cPtr = cellArray->getPointer(0);
    std::vector<double> layerAverages(depth, 0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), LayerAveragesImpl<double>(cPtr, length, width, dStride, lStride, wStride, layerAverages), tbb::auto_partitioner());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, depth - 3), DirectionalDifferencesImpl<double>(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages),
                      tbb::auto_partitioner());
#else
    LayerAveragesImpl<double> serial(cPtr, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth);
    DirectionalDifferencesImpl<double> serial(cPtr, m_ProjectedImageMin, length, width, dStride, lStride, wStride, layerAverages);
    serial.convert(0, depth - 3);
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
AbstractFilter::Pointer FindArrayDifferencesAlongDirection::newFilterInstance(bool copyFilterParameters) const
{
  FindArrayDifferencesAlongDirection::Pointer filter = FindArrayDifferencesAlongDirection::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getGroupName() const
{
  return SIMPL::FilterGroups::ProcessingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FindArrayDifferencesAlongDirection::getUuid() const
{
  return QUuid("{577dfdf6-02f8-5284-b45b-e31f5392a196}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::ImageFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getHumanLabel() const
{
  return "Find Array Differences Along Image Direction";
}

// -----------------------------------------------------------------------------
FindArrayDifferencesAlongDirection::Pointer FindArrayDifferencesAlongDirection::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<FindArrayDifferencesAlongDirection> FindArrayDifferencesAlongDirection::New()
{
  struct make_shared_enabler : public FindArrayDifferencesAlongDirection
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getNameOfClass() const
{
  return QString("FindArrayDifferencesAlongDirection");
}

// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::ClassName()
{
  return QString("FindArrayDifferencesAlongDirection");
}

// -----------------------------------------------------------------------------
void FindArrayDifferencesAlongDirection::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindArrayDifferencesAlongDirection::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void FindArrayDifferencesAlongDirection::setDirection(unsigned int value)
{
  m_Direction = value;
}

// -----------------------------------------------------------------------------
unsigned int FindArrayDifferencesAlongDirection::getDirection() const
{
  return m_Direction;
}

// -----------------------------------------------------------------------------
void FindArrayDifferencesAlongDirection::setProjectedImageMinArrayName(const QString& value)
{
  m_ProjectedImageMinArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayDifferencesAlongDirection::getProjectedImageMinArrayName() const
{
  return m_ProjectedImageMinArrayName;
}