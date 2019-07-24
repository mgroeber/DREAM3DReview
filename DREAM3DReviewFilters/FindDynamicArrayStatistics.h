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
#ifndef _finddynamicarraystatistics_h_
#define _finddynamicarraystatistics_h_

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The FindDynamicArrayStatistics class. See [Filter documentation](@ref finddynamicarraystatistics) for details.
 */
class FindDynamicArrayStatistics : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(FindDynamicArrayStatistics)
  SIMPL_FILTER_NEW_MACRO(FindDynamicArrayStatistics)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(FindDynamicArrayStatistics, AbstractFilter)

  ~FindDynamicArrayStatistics() override;

  SIMPL_FILTER_PARAMETER(bool, FindLength)
  Q_PROPERTY(bool FindLength READ getFindLength WRITE setFindLength)

  SIMPL_FILTER_PARAMETER(bool, FindMin)
  Q_PROPERTY(bool FindMin READ getFindMin WRITE setFindMin)

  SIMPL_FILTER_PARAMETER(bool, FindMax)
  Q_PROPERTY(bool FindMax READ getFindMax WRITE setFindMax)

  SIMPL_FILTER_PARAMETER(bool, FindMean)
  Q_PROPERTY(bool FindMean READ getFindMean WRITE setFindMean)

  SIMPL_FILTER_PARAMETER(bool, FindMedian)
  Q_PROPERTY(bool FindMedian READ getFindMedian WRITE setFindMedian)

  SIMPL_FILTER_PARAMETER(bool, FindStdDeviation)
  Q_PROPERTY(bool FindStdDeviation READ getFindStdDeviation WRITE setFindStdDeviation)

  SIMPL_FILTER_PARAMETER(bool, FindSummation)
  Q_PROPERTY(bool FindSummation READ getFindSummation WRITE setFindSummation)

  SIMPL_FILTER_PARAMETER(DataArrayPath, LengthArrayPath)
  Q_PROPERTY(DataArrayPath LengthArrayPath READ getLengthArrayPath WRITE setLengthArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MinimumArrayPath)
  Q_PROPERTY(DataArrayPath MinimumArrayPath READ getMinimumArrayPath WRITE setMinimumArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MaximumArrayPath)
  Q_PROPERTY(DataArrayPath MaximumArrayPath READ getMaximumArrayPath WRITE setMaximumArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MeanArrayPath)
  Q_PROPERTY(DataArrayPath MeanArrayPath READ getMeanArrayPath WRITE setMeanArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MedianArrayPath)
  Q_PROPERTY(DataArrayPath MedianArrayPath READ getMedianArrayPath WRITE setMedianArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, StdDeviationArrayPath)
  Q_PROPERTY(DataArrayPath StdDeviationArrayPath READ getStdDeviationArrayPath WRITE setStdDeviationArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SummationArrayPath)
  Q_PROPERTY(DataArrayPath SummationArrayPath READ getSummationArrayPath WRITE setSummationArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SelectedArrayPath)
  Q_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)

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
  /**
   * @brief createCompatibleArrays Creates the output statistics arrays with compatible types to the
   * incoming NeighborList
   * @param dataArrayPaths Vector used to store paths for later tuple checking
   */
  template <typename T> void createCompatibleArrays(QVector<DataArrayPath>& dataArrayPaths)
  {
    std::vector<size_t> cDims(1, 1);

    if(m_FindMin)
    {
      m_MinimumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMinimumArrayPath(), 0,
                                                                                                            cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMinimumArrayPath());
      }
    }

    if(m_FindMax)
    {
      m_MaximumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMaximumArrayPath(), 0,
                                                                                                            cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMaximumArrayPath());
      }
    }

    if(m_FindMean)
    {
      m_MeanPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMeanArrayPath(), 0,
                                                                                                         cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMeanArrayPath());
      }
    }

    if(m_FindMedian)
    {
      m_MedianPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMedianArrayPath(), 0,
                                                                                                           cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMedianArrayPath());
      }
    }

    if(m_FindStdDeviation)
    {
      m_StandardDeviationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(
          this, getStdDeviationArrayPath(), 0, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getStdDeviationArrayPath());
      }
    }

    if(m_FindSummation)
    {
      m_SummationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getSummationArrayPath(), 0,
                                                                                                              cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getSummationArrayPath());
      }
    }
  }

  FindDynamicArrayStatistics();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DEFINE_DATAARRAY_VARIABLE(int32_t, Length)
  DEFINE_IDATAARRAY_VARIABLE(Minimum)
  DEFINE_IDATAARRAY_VARIABLE(Maximum)
  DEFINE_IDATAARRAY_VARIABLE(Mean)
  DEFINE_IDATAARRAY_VARIABLE(Median)
  DEFINE_IDATAARRAY_VARIABLE(StandardDeviation)
  DEFINE_IDATAARRAY_VARIABLE(Summation)
  DEFINE_IDATAARRAY_VARIABLE(InputArray)

  FindDynamicArrayStatistics(const FindDynamicArrayStatistics&) = delete; // Copy Constructor Not Implemented
  FindDynamicArrayStatistics(FindDynamicArrayStatistics&&) = delete;      // Move Constructor Not Implemented
  void operator=(const FindDynamicArrayStatistics&) = delete;             // Operator '=' Not Implemented
};

#endif /* _FindDynamicArrayStatistics_H_ */
