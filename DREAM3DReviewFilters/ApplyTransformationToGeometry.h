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

#pragma once

#include <memory>

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/FilterParameters/DynamicTableData.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

/**
 * @brief The ApplyTransformationToGeometry class. See [Filter documentation](@ref applytransformationtogeometry) for details.
 */
class DREAM3DReview_EXPORT ApplyTransformationToGeometry : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(ApplyTransformationToGeometry SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(ApplyTransformationToGeometry)
  PYB11_FILTER_NEW_MACRO(ApplyTransformationToGeometry)
  PYB11_FILTER_PARAMETER(DynamicTableData, ManualTransformationMatrix)
  PYB11_FILTER_PARAMETER(DataArrayPath, ComputedTransformationMatrix)
  PYB11_FILTER_PARAMETER(DataArrayPath, GeometryToTransform)
  PYB11_FILTER_PARAMETER(int, TransformationMatrixType)
  PYB11_FILTER_PARAMETER(FloatVec3Type, RotationAxis)
  PYB11_FILTER_PARAMETER(float, RotationAngle)
  PYB11_FILTER_PARAMETER(FloatVec3Type, Translation)
  PYB11_FILTER_PARAMETER(FloatVec3Type, Scale)
  PYB11_PROPERTY(DynamicTableData ManualTransformationMatrix READ getManualTransformationMatrix WRITE setManualTransformationMatrix)
  PYB11_PROPERTY(DataArrayPath ComputedTransformationMatrix READ getComputedTransformationMatrix WRITE setComputedTransformationMatrix)
  PYB11_PROPERTY(DataArrayPath GeometryToTransform READ getGeometryToTransform WRITE setGeometryToTransform)
  PYB11_PROPERTY(int TransformationMatrixType READ getTransformationMatrixType WRITE setTransformationMatrixType)
  PYB11_PROPERTY(FloatVec3Type RotationAxis READ getRotationAxis WRITE setRotationAxis)
  PYB11_PROPERTY(float RotationAngle READ getRotationAngle WRITE setRotationAngle)
  PYB11_PROPERTY(FloatVec3Type Translation READ getTranslation WRITE setTranslation)
  PYB11_PROPERTY(FloatVec3Type Scale READ getScale WRITE setScale)
#endif

public:
  using Self = ApplyTransformationToGeometry;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ApplyTransformationToGeometry> New();

  /**
   * @brief Returns the name of the class for ApplyTransformationToGeometry
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ApplyTransformationToGeometry
   */
  static QString ClassName();

  ~ApplyTransformationToGeometry() override;

  /**
   * @brief Setter property for ManualTransformationMatrix
   */
  void setManualTransformationMatrix(const DynamicTableData& value);
  /**
   * @brief Getter property for ManualTransformationMatrix
   * @return Value of ManualTransformationMatrix
   */
  DynamicTableData getManualTransformationMatrix() const;

  Q_PROPERTY(DynamicTableData ManualTransformationMatrix READ getManualTransformationMatrix WRITE setManualTransformationMatrix)

  /**
   * @brief Setter property for ComputedTransformationMatrix
   */
  void setComputedTransformationMatrix(const DataArrayPath& value);
  /**
   * @brief Getter property for ComputedTransformationMatrix
   * @return Value of ComputedTransformationMatrix
   */
  DataArrayPath getComputedTransformationMatrix() const;

  Q_PROPERTY(DataArrayPath ComputedTransformationMatrix READ getComputedTransformationMatrix WRITE setComputedTransformationMatrix)

  /**
   * @brief Setter property for GeometryToTransform
   */
  void setGeometryToTransform(const DataArrayPath& value);
  /**
   * @brief Getter property for GeometryToTransform
   * @return Value of GeometryToTransform
   */
  DataArrayPath getGeometryToTransform() const;

  Q_PROPERTY(DataArrayPath GeometryToTransform READ getGeometryToTransform WRITE setGeometryToTransform)

  /**
   * @brief Setter property for TransformationMatrixType
   */
  void setTransformationMatrixType(int value);
  /**
   * @brief Getter property for TransformationMatrixType
   * @return Value of TransformationMatrixType
   */
  int getTransformationMatrixType() const;

  Q_PROPERTY(int TransformationMatrixType READ getTransformationMatrixType WRITE setTransformationMatrixType)

  /**
   * @brief Setter property for RotationAxis
   */
  void setRotationAxis(const FloatVec3Type& value);
  /**
   * @brief Getter property for RotationAxis
   * @return Value of RotationAxis
   */
  FloatVec3Type getRotationAxis() const;

  Q_PROPERTY(FloatVec3Type RotationAxis READ getRotationAxis WRITE setRotationAxis)

  /**
   * @brief Setter property for RotationAngle
   */
  void setRotationAngle(float value);
  /**
   * @brief Getter property for RotationAngle
   * @return Value of RotationAngle
   */
  float getRotationAngle() const;

  Q_PROPERTY(float RotationAngle READ getRotationAngle WRITE setRotationAngle)

  /**
   * @brief Setter property for Translation
   */
  void setTranslation(const FloatVec3Type& value);
  /**
   * @brief Getter property for Translation
   * @return Value of Translation
   */
  FloatVec3Type getTranslation() const;

  Q_PROPERTY(FloatVec3Type Translation READ getTranslation WRITE setTranslation)

  /**
   * @brief Setter property for Scale
   */
  void setScale(const FloatVec3Type& value);
  /**
   * @brief Getter property for Scale
   * @return Value of Scale
   */
  FloatVec3Type getScale() const;

  Q_PROPERTY(FloatVec3Type Scale READ getScale WRITE setScale)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
  */
  QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

protected:
  ApplyTransformationToGeometry();

  /**
   * @brief applyTransformation
   */
  void applyTransformation();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

private:
  std::weak_ptr<DataArray<float>> m_TransformationMatrixPtr;
  float* m_TransformationMatrix = nullptr;

  DynamicTableData m_ManualTransformationMatrix = {};
  DataArrayPath m_ComputedTransformationMatrix = {};
  DataArrayPath m_GeometryToTransform = {};
  int m_TransformationMatrixType = {};
  FloatVec3Type m_RotationAxis = {};
  float m_RotationAngle = {};
  FloatVec3Type m_Translation = {};
  FloatVec3Type m_Scale = {};

  FloatArrayType::Pointer m_TransformationReference;

public:
  ApplyTransformationToGeometry(const ApplyTransformationToGeometry&) = delete; // Copy Constructor Not Implemented
  ApplyTransformationToGeometry(ApplyTransformationToGeometry&&) = delete;      // Move Constructor Not Implemented
  ApplyTransformationToGeometry& operator=(const ApplyTransformationToGeometry&) = delete; // Copy Assignment Not Implemented
  ApplyTransformationToGeometry& operator=(ApplyTransformationToGeometry&&) = delete;      // Move Assignment Not Implemented
};

