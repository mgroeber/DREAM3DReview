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
#include "FindRidgeLines.h"

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
class FindRidgesImpl
{

public:
  FindRidgesImpl(T *data, DataContainer::Pointer geom, bool* ridgeFlag)
  : m_Data(data)
  , m_Geom(geom)
  , m_RidgeFlag(ridgeFlag)
  {
  }
  virtual ~FindRidgesImpl() = default;

  void convert(size_t start, size_t end) const
  {

    SizeVec3Type geoDims = m_Geom->getGeometryAs<ImageGeom>()->getDimensions();

    std::vector<size_t> neighbors(6);
    neighbors[0] = -(geoDims[0] * geoDims[1]);
    neighbors[1] = (geoDims[0] * geoDims[1]);
    neighbors[2] = -(geoDims[0]);
    neighbors[3] = (geoDims[0]);
    neighbors[4] = -1;
    neighbors[5] = 1;

    for(size_t iter = start; iter < end; iter++)
    {
      size_t col = iter % geoDims[0];
      size_t row = (iter / geoDims[0]) % geoDims[1];
      size_t plane = iter / (geoDims[0] * geoDims[1]);

      bool xflag = true;
      bool yflag = true;
      bool zflag = true;

      for(size_t i = 0; i < 2; i++)
      {
        if(i == 0 && plane == 0)
        {
          continue;
        }
        if(i == 1 && plane == (geoDims[2] - 1))
        {
          continue;
        }
        if(m_Data[iter] < m_Data[iter + neighbors[i]])
        {
          zflag = false;
        }
      }
      for(size_t i = 2; i < 4; i++)
      {
        if(i == 2 && row == 0)
        {
          continue;
        }
        if(i == 3 && row == (geoDims[1] - 1))
        {
          continue;
        }
        if(m_Data[iter] < m_Data[iter + neighbors[i]])
        {
          yflag = false;
        }
      }
      for(size_t i = 4; i < 6; i++)
      {
        if(i == 4 && col == 0)
        {
          continue;
        }
        if(i == 5 && col == (geoDims[0] - 1))
        {
          continue;
        }
        if(m_Data[iter] < m_Data[iter + neighbors[i]])
        {
          xflag = false;
        }
      }
      if(xflag && yflag && zflag)
      {
        m_RidgeFlag[iter] = true;
      }
      else
      {
        m_RidgeFlag[iter] = false;
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
  DataContainer::Pointer m_Geom;
  bool* m_RidgeFlag;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindRidgeLines::FindRidgeLines() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindRidgeLines::~FindRidgeLines() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindRidgeLines::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Direction of Interest");
    parameter->setPropertyName("Direction");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(FindRidgeLines, this, Direction));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(FindRidgeLines, this, Direction));

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
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Quantify", SelectedArrayPath, FilterParameter::Category::RequiredArray, FindRidgeLines, req));
  }
  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Ridge Flag Array", RidgeFlagArrayName, SelectedArrayPath, SelectedArrayPath, FilterParameter::Category::CreatedArray, FindRidgeLines));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
void FindRidgeLines::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setRidgeFlagArrayName(reader->readString("RidgeFlagArrayName", getRidgeFlagArrayName()));
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setDirection(reader->readValue("Direction", getDirection()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindRidgeLines::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindRidgeLines::dataCheck()
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
  tempPath.update(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getRidgeFlagArrayName());
  m_RidgeFlagPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, tempPath, false, cDims);
  if(nullptr != m_RidgeFlagPtr.lock())
  {
    m_RidgeFlag = m_RidgeFlagPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  
  ImageGeom::Pointer image = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom>(this, getSelectedArrayPath().getDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindRidgeLines::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedArrayPath().getDataContainerName());

  SizeVec3Type geoDims = m->getGeometryAs<ImageGeom>()->getDimensions();
  size_t totalPoints = geoDims[0] * geoDims[1] * geoDims[2];

  IDataArray::Pointer iCellArray = m_InDataPtr.lock();


  if(TemplateHelpers::CanDynamicCast<Int8ArrayType>()(m_InDataPtr.lock()))
  {
    Int8ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int8ArrayType>(m_InDataPtr.lock());
    int8_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<int8_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());

#else
    FindRidgesImpl<int8_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt8ArrayType>()(m_InDataPtr.lock()))
  {
    UInt8ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt8ArrayType>(m_InDataPtr.lock());
    uint8_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<uint8_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());

#else
    FindRidgesImpl<uint8_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int16ArrayType>()(m_InDataPtr.lock()))
  {
    Int16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int16ArrayType>(m_InDataPtr.lock());
    int16_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS 
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<int16_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());
#else
DirectionalDifferencesImpl<int16_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt16ArrayType>()(m_InDataPtr.lock()))
  {
    UInt16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt16ArrayType>(m_InDataPtr.lock());
    uint16_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<uint16_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());
#else
    FindRidgesImpl<uint16_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int32ArrayType>()(m_InDataPtr.lock()))
  {
    Int32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int32ArrayType>(m_InDataPtr.lock());
    int32_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<int32_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());

#else
    FindRidgesImpl<int32_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt32ArrayType>()(m_InDataPtr.lock()))
  {
    UInt32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt32ArrayType>(m_InDataPtr.lock());
    uint32_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS 
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<uint32_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());

#else
    FindRidgesImpl<uint32_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<Int64ArrayType>()(m_InDataPtr.lock()))
  {
    Int64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int64ArrayType>(m_InDataPtr.lock());
    int64_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<int64_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());

#else
    FindRidgesImpl<int64_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<UInt64ArrayType>()(m_InDataPtr.lock()))
  {
    UInt64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt64ArrayType>(m_InDataPtr.lock());
    uint64_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<uint64_t>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());
#else
    FindRidgesImpl<uint64_t> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<FloatArrayType>()(m_InDataPtr.lock()))
  {
    FloatArrayType::Pointer cellArray = std::dynamic_pointer_cast<FloatArrayType>(m_InDataPtr.lock());
    float* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<float>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());
#else
    FindRidgesImpl<float> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
#endif
  }
  else if(TemplateHelpers::CanDynamicCast<DoubleArrayType>()(m_InDataPtr.lock()))
  {
    DoubleArrayType::Pointer cellArray = std::dynamic_pointer_cast<DoubleArrayType>(m_InDataPtr.lock());
    double* cPtr = cellArray->getPointer(0);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPoints), FindRidgesImpl<double>(cPtr, m, m_RidgeFlag),
                      tbb::auto_partitioner());
#else
    FindRidgesImpl<double> serial(cPtr, m, m_RidgeFlag);
    serial.convert(0, totalPoints);
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
AbstractFilter::Pointer FindRidgeLines::newFilterInstance(bool copyFilterParameters) const
{
  FindRidgeLines::Pointer filter = FindRidgeLines::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindRidgeLines::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindRidgeLines::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindRidgeLines::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindRidgeLines::getGroupName() const
{
  return SIMPL::FilterGroups::ProcessingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FindRidgeLines::getUuid() const
{
  return QUuid("{577dfdf6-02f8-5284-b45b-e31f5392a178}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindRidgeLines::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::ImageFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindRidgeLines::getHumanLabel() const
{
  return "Find Ridge Lines";
}

// -----------------------------------------------------------------------------
FindRidgeLines::Pointer FindRidgeLines::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<FindRidgeLines> FindRidgeLines::New()
{
  struct make_shared_enabler : public FindRidgeLines
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString FindRidgeLines::getNameOfClass() const
{
  return QString("FindRidgeLines");
}

// -----------------------------------------------------------------------------
QString FindRidgeLines::ClassName()
{
  return QString("FindRidgeLines");
}

// -----------------------------------------------------------------------------
void FindRidgeLines::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindRidgeLines::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void FindRidgeLines::setDirection(unsigned int value)
{
  m_Direction = value;
}

// -----------------------------------------------------------------------------
unsigned int FindRidgeLines::getDirection() const
{
  return m_Direction;
}

// -----------------------------------------------------------------------------
void FindRidgeLines::setRidgeFlagArrayName(const QString& value)
{
  m_RidgeFlagArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindRidgeLines::getRidgeFlagArrayName() const
{
  return m_RidgeFlagArrayName;
}