/* ============================================================================
 * Copyright (c) 2016 Czech Academy of Sciences, Institute of Physics,
 * Group of Bulk Nanomaterials and Interfaces, http://ams.fzu.cz
 *
 * Redistribution and use in source and binary forms, with or without
 * modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of the Czech Academy of Sciences, nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
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
 * The code contained herein was partially funded by the followig grants:
 *    Czech Science Foundation (GA CR), project no. GBP108/12/G043
 *    Czech Ministry of Education, Youth and Sports (MSMT), project no. LM2015087
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"

#include "DREAM3DReview/DREAM3DReviewFilters/AdaptiveAlignment.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The AdaptiveAlignmentFeature class. See [Filter documentation](@ref AdaptiveAlignmentfeature) for details.
 */
class DREAM3DReview_EXPORT AdaptiveAlignmentFeature : public AdaptiveAlignment
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(AdaptiveAlignmentFeature SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(AdaptiveAlignmentFeature)
  PYB11_FILTER_NEW_MACRO(AdaptiveAlignmentFeature)
  PYB11_FILTER_PARAMETER(DataArrayPath, GoodVoxelsArrayPath)
  PYB11_PROPERTY(DataArrayPath GoodVoxelsArrayPath READ getGoodVoxelsArrayPath WRITE setGoodVoxelsArrayPath)
#endif

public:
  using Self = AdaptiveAlignmentFeature;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<AdaptiveAlignmentFeature> New();

  /**
   * @brief Returns the name of the class for AdaptiveAlignmentFeature
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for AdaptiveAlignmentFeature
   */
  static QString ClassName();

  virtual ~AdaptiveAlignmentFeature();

  /**
   * @brief Setter property for GoodVoxelsArrayPath
   */
  void setGoodVoxelsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for GoodVoxelsArrayPath
   * @return Value of GoodVoxelsArrayPath
   */
  DataArrayPath getGoodVoxelsArrayPath() const;

  Q_PROPERTY(DataArrayPath GoodVoxelsArrayPath READ getGoodVoxelsArrayPath WRITE setGoodVoxelsArrayPath)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is
   * a tag
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
  AdaptiveAlignmentFeature();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and
   * availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief find_shifts Reimplemented from @see AdaptiveAlignment class
   * @param xshifts
   * @param yshifts
   * @param xneedshifts
   * @param yneedshifts
   */
  virtual void find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts);

  /**
   * @brief compute_error1 Determines error between the current and desired shifts
   * (based on discrepancy of slopes)
   * @param iter
   * @param index
   * @param xneedtrend
   * @param yneedtrend
   * @param newxshift
   * @param newyshift
   * @param curindex
   * @return
   */
  virtual float compute_error1(uint64_t iter, uint64_t index, float xneedtrend, float yneedtrend, std::vector<std::vector<int64_t>>& newxshift, std::vector<std::vector<int64_t>>& newyshift,
                               std::vector<uint64_t>& curindex);

  /**
   * @brief compute_error1 Determines error between the current and desired shifts
   * (based on discrepancy from the estimates obtained from SEM images)
   * @param iter
   * @param index
   * @param xshiftsest
   * @param yshiftsest
   * @param newxshift
   * @param newyshift
   * @param curindex
   * @return
   */
  virtual float compute_error2(uint64_t iter, uint64_t index, std::vector<float>& xshiftsest, std::vector<float>& yshiftsest, std::vector<std::vector<int64_t>>& newxshift,
                               std::vector<std::vector<int64_t>>& newyshift, std::vector<uint64_t>& curindex);

private:
  std::weak_ptr<DataArray<bool>> m_GoodVoxelsPtr;
  bool* m_GoodVoxels = nullptr;

  DataArrayPath m_GoodVoxelsArrayPath = {};

public:
  AdaptiveAlignmentFeature(const AdaptiveAlignmentFeature&) = delete;            // Copy Constructor Not Implemented
  AdaptiveAlignmentFeature(AdaptiveAlignmentFeature&&) = delete;                 // Move Constructor Not Implemented
  AdaptiveAlignmentFeature& operator=(const AdaptiveAlignmentFeature&) = delete; // Copy Assignment Not Implemented
  AdaptiveAlignmentFeature& operator=(AdaptiveAlignmentFeature&&) = delete;      // Move Assignment Not Implemented
};
