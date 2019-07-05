#ifndef _CutPatchesFromImage_h_
#define _CutPatchesFromImage_h_

#include <array>
#include <unordered_map>

#include <QMutex>
#include "H5Support/QH5Utilities.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/FilterParameters/IntVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

#include "SIMPLib/Geometry/ImageGeom.h"

/**
 * @brief The CutPatchesFromImage class. See [Filter documentation](@ref CutPatchesFromImage) for details.
 */
class CutPatchesFromImage : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(CutPatchesFromImage)
  SIMPL_FILTER_NEW_MACRO(CutPatchesFromImage)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(CutPatchesFromImage, AbstractFilter)

  ~CutPatchesFromImage() override;

  SIMPL_FILTER_PARAMETER(QString, OutputDirectory)
  Q_PROPERTY(QString OutputDirectory READ getOutputDirectory WRITE setOutputDirectory)

  SIMPL_FILTER_PARAMETER(QVector<DataArrayPath>, SelectedDataArrayPaths)
  Q_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)

  SIMPL_FILTER_PARAMETER(IntVec3Type, PatchDims)
  Q_PROPERTY(IntVec3Type PatchDims READ getPatchDims WRITE setPatchDims)

  SIMPL_FILTER_PARAMETER(int, NumPatches)
  Q_PROPERTY(int NumPatches READ getNumPatches WRITE setNumPatches)

  SIMPL_FILTER_PARAMETER(bool, BalancePatches)
  Q_PROPERTY(bool BalancePatches READ getBalancePatches WRITE setBalancePatches)

  SIMPL_FILTER_PARAMETER(DataArrayPath, LabelsArrayPath)
  Q_PROPERTY(DataArrayPath LabelsArrayPath READ getLabelsArrayPath WRITE setLabelsArrayPath)

    /**
     * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
     */
  const QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
  */
  const QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  const QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getSubGroupName() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  const QString getHumanLabel() const override;

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
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  void preflight() override;

  /**
  * @brief getUuid Return the unique identifier for this filter.
  * @return A QUuid object.
  */
  const QUuid getUuid() override;

signals:
  /**
   * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
   * be pushed from a user-facing control (such as a widget)
   * @param filter Filter instance pointer
   */
  void updateFilterParameters(AbstractFilter* filter);

  /**
   * @brief parametersChanged Emitted when any Filter parameter is changed internally
   */
  void parametersChanged();

  /**
   * @brief preflightAboutToExecute Emitted just before calling dataCheck()
   */
  void preflightAboutToExecute();

  /**
   * @brief preflightExecuted Emitted just after calling dataCheck()
   */
  void preflightExecuted();

protected:
  CutPatchesFromImage();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  void sendThreadSafeProgressMessage(int64_t counter);

private:
  DEFINE_DATAARRAY_VARIABLE(bool, LabelsArray)

  QVector<IDataArray::WeakPointer> m_SelectedWeakPtrVector;

  QMutex m_Mutex;
  int64_t m_ProgressCounter;
  int64_t m_TotalElements;
  int64_t m_LastProgressInt;

  friend class CutPatchesImpl;

  CutPatchesFromImage(const CutPatchesFromImage&) = delete; // Copy Constructor Not Implemented
  CutPatchesFromImage(CutPatchesFromImage&&) = delete;      // Move Constructor Not Implemented
  void operator=(const CutPatchesFromImage&) = delete;                           // Operator '=' Not Implemented
};

#endif /* _CutPatchesFromImage_h_ */
