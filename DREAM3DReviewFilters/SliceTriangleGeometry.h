/* ============================================================================
 * Software developed by US federal government employees (including military personnel)
 * as part of their official duties is not subject to copyright protection and is
 * considered "public domain" (see 17 USC Section 105). Public domain software can be used
 * by anyone for any purpose, and cannot be released under a copyright license
 * (including typical open source software licenses).
 *
 * This source code file was originally written by United States DoD employees. The
 * original source code files are released into the Public Domain.
 *
 * Subsequent changes to the codes by others may elect to add a copyright and license
 * for those changes.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/IGeometry.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The SliceTriangleGeometry class. See [Filter documentation](@ref SliceTriangleGeometry) for details.
 */
class DREAM3DReview_EXPORT SliceTriangleGeometry : public AbstractFilter
{
  Q_OBJECT

  // clang-format off
  PYB11_BEGIN_BINDINGS(SliceTriangleGeometry SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(SliceTriangleGeometry)
  PYB11_FILTER_NEW_MACRO(SliceTriangleGeometry)

  PYB11_PROPERTY(DataArrayPath CADDataContainerName READ getCADDataContainerName WRITE setCADDataContainerName)
  PYB11_PROPERTY(QString SliceDataContainerName READ getSliceDataContainerName WRITE setSliceDataContainerName)
  PYB11_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)
  PYB11_PROPERTY(QString SliceAttributeMatrixName READ getSliceAttributeMatrixName WRITE setSliceAttributeMatrixName)
  PYB11_PROPERTY(QString SliceIdArrayName READ getSliceIdArrayName WRITE setSliceIdArrayName)
  PYB11_PROPERTY(bool HaveRegionIds READ getHaveRegionIds WRITE setHaveRegionIds)
  PYB11_PROPERTY(bool DoRotation READ getDoRotation WRITE setDoRotation)
  PYB11_PROPERTY(DataArrayPath RegionIdArrayPath READ getRegionIdArrayPath WRITE setRegionIdArrayPath)
  PYB11_PROPERTY(FloatVec3Type SliceDirection READ getSliceDirection WRITE setSliceDirection)
  PYB11_PROPERTY(float SliceResolution READ getSliceResolution WRITE setSliceResolution)
  PYB11_PROPERTY(float Zstart READ getZstart WRITE setZstart)
  PYB11_PROPERTY(float Zend READ getZend WRITE setZend)
  PYB11_PROPERTY(int SliceRange READ getSliceRange WRITE setSliceRange)
  PYB11_END_BINDINGS()
  // clang-format on

public:
  using Self = SliceTriangleGeometry;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<SliceTriangleGeometry> New();

  /**
   * @brief Returns the name of the class for SliceTriangleGeometry
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for SliceTriangleGeometry
   */
  static QString ClassName();

  ~SliceTriangleGeometry() override;

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
   * @brief Setter property for HaveRegionIds
   */
  void setDoRotation(bool value);
  /**
   * @brief Getter property for HaveRegionIds
   * @return Value of HaveRegionIds
   */

  bool getDoRotation() const;
  //don't actually want to display this to user, so no need for a Q_PROPERTY
  //Q_PROPERTY(bool DoRotation READ getDoRotation WRITE setDoRotation)

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
   * @brief Setter property for SliceDirection
   */
  void setSliceDirection(const FloatVec3Type& value);
  /**
   * @brief Getter property for SliceDirection
   * @return Value of SliceDirection
   */
  FloatVec3Type getSliceDirection() const;
  Q_PROPERTY(FloatVec3Type SliceDirection READ getSliceDirection WRITE setSliceDirection)

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
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

protected:
  SliceTriangleGeometry();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief rotateVertices either rotates CAD triangles and/or edges into or out of the sectioning orientation so the sectioning is always along 001
   */
  void rotateVertices(unsigned int direction, float* n, int64_t numVerts, float* verts);

  /**
   * @brief Determines if a segment between two points intersects a plane defined by a normal and distance
   * @param n
   * @param d
   * @param q
   * @param r
   * @param p
   * @return
   */
  static char rayIntersectsPlane(const float d, const float* q, const float* r, float* p);

  /**
   * @brief updateEdgeInstancePointers
   */
  void determineBoundsAndNumSlices(float& minDim, float& maxDim, MeshIndexType numTris, MeshIndexType* tris, float* triVerts);

private:
  std::weak_ptr<Int32ArrayType> m_SliceIdPtr;
  std::weak_ptr<Int32ArrayType> m_RegionIdPtr;
  std::weak_ptr<Int32ArrayType> m_TriRegionIdPtr;

  DataArrayPath m_CADDataContainerName = {"TriangleDataContainer", "", ""};
  QString m_SliceDataContainerName = {"SliceDataContainer"};
  QString m_EdgeAttributeMatrixName = {"EdgeData"};
  QString m_SliceAttributeMatrixName = {"SliceData"};
  QString m_SliceIdArrayName = {"SliceIds"};

  bool m_HaveRegionIds = {false};
  bool m_DoRotation = {false};
  DataArrayPath m_RegionIdArrayPath = {"", "", ""};
  FloatVec3Type m_SliceDirection = {0.0F, 0.0F, 1.0F};
  float m_SliceResolution = {1.0F};
  float m_Zstart = {0.0F};
  float m_Zend = {0.0F};
  int m_SliceRange = {0};

  enum RotationDirection
  {
    rotForward,
    rotBackward
  };

  int32_t m_NumberOfSlices = 0;

  SliceTriangleGeometry(const SliceTriangleGeometry&) = delete; // Copy Constructor Not Implemented
  SliceTriangleGeometry(SliceTriangleGeometry&&) = delete;      // Move Constructor Not Implemented
  void operator=(const SliceTriangleGeometry&) = delete;        // Operator '=' Not Implemented
};
