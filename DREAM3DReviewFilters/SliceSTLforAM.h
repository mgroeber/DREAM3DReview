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
#ifndef _SliceSTLforAM_h_
#define _SliceSTLforAM_h_

#include <memory>

#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Geometry/IGeometry.h"

#include "DREAM3DReview\DREAM3DReviewFilters\SliceTriangleGeometry.h"

/**
 * @brief The SliceSTLforAM class. See [Filter documentation](@ref SliceSTLforAM) for details.
 */
class SliceSTLforAM : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = SliceSTLforAM;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<SliceSTLforAM> New();

  /**
   * @brief Returns the name of the class for SliceSTLforAM
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for SliceSTLforAM
   */
  static QString ClassName();

  ~SliceSTLforAM() override;

  /**
   * @brief Setter property for CADDataContainerName
   */
  void setCADDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for CADDataContainerName
   * @return Value of CADDataContainerName
   */
  DataArrayPath getCADDataContainerName() const;

  Q_PROPERTY(DataArrayPath CADDataContainerName READ getCADDataContainerName WRITE setCADDataContainerName)

  /**
   * @brief Setter property for SliceDataContainerName
   */
  void setSliceDataContainerName(const QString& value);
  /**
   * @brief Getter property for SliceDataContainerName
   * @return Value of SliceDataContainerName
   */
  QString getSliceDataContainerName() const;

  Q_PROPERTY(QString SliceDataContainerName READ getSliceDataContainerName WRITE setSliceDataContainerName)

  /**
   * @brief Setter property for EdgeAttributeMatrixName
   */
  void setEdgeAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for EdgeAttributeMatrixName
   * @return Value of EdgeAttributeMatrixName
   */
  QString getEdgeAttributeMatrixName() const;

  Q_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)

  /**
   * @brief Setter property for SliceAttributeMatrixName
   */
  void setSliceAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for SliceAttributeMatrixName
   * @return Value of SliceAttributeMatrixName
   */
  QString getSliceAttributeMatrixName() const;

  Q_PROPERTY(QString SliceAttributeMatrixName READ getSliceAttributeMatrixName WRITE setSliceAttributeMatrixName)

  /**
   * @brief Setter property for SliceIdArrayName
   */
  void setSliceIdArrayName(const QString& value);
  /**
   * @brief Getter property for SliceIdArrayName
   * @return Value of SliceIdArrayName
   */
  QString getSliceIdArrayName() const;

  Q_PROPERTY(QString SliceIdArrayName READ getSliceIdArrayName WRITE setSliceIdArrayName)

  /**
   * @brief Setter property for AreasArrayName
   */
  void setAreasArrayName(const QString& value);
  /**
   * @brief Getter property for AreasArrayName
   * @return Value of AreasArrayName
   */
  QString getAreasArrayName() const;

  Q_PROPERTY(QString AreasArrayName READ getAreasArrayName WRITE setAreasArrayName)

  /**
   * @brief Setter property for PerimetersArrayName
   */
  void setPerimetersArrayName(const QString& value);
  /**
   * @brief Getter property for PerimetersArrayName
   * @return Value of PerimetersArrayName
   */
  QString getPerimetersArrayName() const;

  Q_PROPERTY(QString PerimetersArrayName READ getPerimetersArrayName WRITE setPerimetersArrayName)

  /**
   * @brief Setter property for HaveRegionIds
   */
  void setHaveRegionIds(bool value);
  /**
   * @brief Getter property for HaveRegionIds
   * @return Value of HaveRegionIds
   */
  bool getHaveRegionIds() const;

  Q_PROPERTY(bool HaveRegionIds READ getHaveRegionIds WRITE setHaveRegionIds)

  /**
   * @brief Setter property for RegionIdArrayPath
   */
  void setRegionIdArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for RegionIdArrayPath
   * @return Value of RegionIdArrayPath
   */
  DataArrayPath getRegionIdArrayPath() const;

  Q_PROPERTY(DataArrayPath RegionIdArrayPath READ getRegionIdArrayPath WRITE setRegionIdArrayPath)

  /**
   * @brief Setter property for SliceResolution
   */
  void setSliceResolution(float value);
  /**
   * @brief Getter property for SliceResolution
   * @return Value of SliceResolution
   */
  float getSliceResolution() const;

  Q_PROPERTY(float SliceResolution READ getSliceResolution WRITE setSliceResolution)

  /**
   * @brief Setter property for Zstart
   */
  void setZstart(float value);
  /**
   * @brief Getter property for Zstart
   * @return Value of Zstart
   */
  float getZstart() const;

  Q_PROPERTY(float Zstart READ getZstart WRITE setZstart)

  /**
   * @brief Setter property for Zend
   */
  void setZend(float value);
  /**
   * @brief Getter property for Zend
   * @return Value of Zend
   */
  float getZend() const;

  Q_PROPERTY(float Zend READ getZend WRITE setZend)

  /**
   * @brief Setter property for SliceRange
   */
  void setSliceRange(int value);
  /**
   * @brief Getter property for SliceRange
   * @return Value of SliceRange
   */
  int getSliceRange() const;

  Q_PROPERTY(int SliceRange READ getSliceRange WRITE setSliceRange)

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

  /**
  * @brief getUuid Return the unique identifier for this filter.
  * @return A QUuid object.
  */
  QUuid getUuid() const override;

protected:
  SliceSTLforAM();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
  * @brief updateSliceInstancePointers
  */
  void updateSliceInstancePointers();

private:
  std::weak_ptr<DataArray<int32_t>> m_SliceIdPtr;
  int32_t* m_SliceId = nullptr;
  std::weak_ptr<DataArray<float>> m_AreaPtr;
  float* m_Area = nullptr;
  std::weak_ptr<DataArray<float>> m_PerimeterPtr;
  float* m_Perimeter = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_RegionIdPtr;
  int32_t* m_RegionId = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_TriRegionIdPtr;
  int32_t* m_TriRegionId = nullptr;

  DataArrayPath m_CADDataContainerName = {};
  QString m_SliceDataContainerName = {};
  QString m_EdgeAttributeMatrixName = {};
  QString m_SliceAttributeMatrixName = {};
  QString m_SliceIdArrayName = {};
  QString m_AreasArrayName = {};
  QString m_PerimetersArrayName = {};
  bool m_HaveRegionIds = {};
  DataArrayPath m_RegionIdArrayPath = {};
  float m_SliceResolution = {};
  float m_Zstart = {};
  float m_Zend = {};
  int m_SliceRange = {};

  enum RotationDirection
  {
    rotForward,
    rotBackward
  };

  int32_t m_NumberOfSlices = 0;

  SliceSTLforAM(const SliceSTLforAM&) = delete;   // Copy Constructor Not Implemented
  SliceSTLforAM(SliceSTLforAM&&) = delete;        // Move Constructor Not Implemented
  void operator=(const SliceSTLforAM&) = delete; // Operator '=' Not Implemented
};

#endif /* SliceSTLforAM */
