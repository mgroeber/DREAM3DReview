#ifndef _ThermoElasticFFT_h_
#define _ThermoElasticFFT_h_

#include <array>
#include <unordered_map>

#include <QMutex>
#include "H5Support/QH5Utilities.h"
#include "SIMPLib/CoreFilters/FileReader.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/FilterParameters/IntVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

#include "SIMPLib/Geometry/ImageGeom.h"

/**
 * @brief The ThermoElasticFFT class. See [Filter documentation](@ref ThermoElasticFFT) for details.
 */
class ThermoElasticFFT : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(ThermoElasticFFT)
  SIMPL_FILTER_NEW_MACRO(ThermoElasticFFT)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(ThermoElasticFFT, AbstractFilter)

  ~ThermoElasticFFT() override;

  SIMPL_FILTER_PARAMETER(QString, SetupFile)
  Q_PROPERTY(QString SetupFile READ getSetupFile WRITE setSetupFile)

  SIMPL_FILTER_PARAMETER(QString, ModuliFile)
  Q_PROPERTY(QString ModuliFile READ getModuliFile WRITE setModuliFile)

  SIMPL_FILTER_PARAMETER(DataArrayPath, FeatureIdsArrayPath)
  Q_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, CellPhasesArrayPath)
  Q_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, EulerAnglesArrayPath)
  Q_PROPERTY(DataArrayPath EulerAnglesArrayPath READ getEulerAnglesArrayPath WRITE setEulerAnglesArrayPath)

  SIMPL_FILTER_PARAMETER(QString, StressArrayName)
  Q_PROPERTY(QString StressArrayName READ getStressArrayName WRITE setStressArrayName)

  SIMPL_FILTER_PARAMETER(QString, StrainArrayName)
  Q_PROPERTY(QString StrainArrayName READ getStrainArrayName WRITE setStrainArrayName)

  SIMPL_FILTER_PARAMETER(bool, SaveEigenStresses)
  Q_PROPERTY(bool SaveEigenStresses READ getSaveEigenStresses WRITE setSaveEigenStresses)

  SIMPL_FILTER_PARAMETER(QString, EigenStressArrayName)
  Q_PROPERTY(QString EigenStressArrayName READ getEigenStressArrayName WRITE setEigenStressArrayName)

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
  ThermoElasticFFT();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DEFINE_DATAARRAY_VARIABLE(int32_t, FeatureIds)
  DEFINE_DATAARRAY_VARIABLE(int32_t, CellPhases)
  DEFINE_DATAARRAY_VARIABLE(float, EulerAngles)
  DEFINE_DATAARRAY_VARIABLE(float, EigenStress)
  DEFINE_DATAARRAY_VARIABLE(float, Stress)
  DEFINE_DATAARRAY_VARIABLE(float, Strain)

  FILE* m_InStream1;
  FILE* m_InStream2;

  int32_t npts1;
  int32_t npts2;
  int32_t npts3;

  float dbar[6], sbar[6], delt[3], xlsec[6][6], xlsec33[3][3][3][3], cmat[6][6];
  std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> cmat33;
  float sref, dref, errs, errd;
  float ctexlsec[3][3], cteaux33[3][3], ctebar[6], dtilbar[6];
  std::vector<std::vector<float>> xlsecinv;

  float wgt; 
  std::vector<std::vector<float>> sg, dtilde, cteloc6;
  std::vector < std::vector<std::vector<float>>> velgrad, velgradim, cloc, fsloc;
  

  float udot[3][3], dsim[3][3], stens[3][3], tomtot[3][3], iudot[3][3], idsim[6], dtens[3][3];
  float error, deltat;
  int32_t ictrl, ictrl1, ictrl2, nsteps, itmax, thermctrl;
  std::vector<std::vector<std::vector<float>>> ctemat;

  void fourn(std::vector<float>& data, int32_t* nn, int32_t ndim, int32_t isign);
  void voigt(float c2[6][6], float c4[3][3][3][3], int iopt);
  void input();
  void data_crystal();
  void data_grain();
  void chg_basis(float ce2[6], float c2[3][3], float ce4[6][6], float c4[3][3][3][3], int iopt);
  void augm_lagr();
  float tnorm(float v[36], int32_t nrows, int32_t ncols);

  ThermoElasticFFT(const ThermoElasticFFT&) = delete; // Copy Constructor Not Implemented
  ThermoElasticFFT(ThermoElasticFFT&&) = delete;      // Move Constructor Not Implemented
  void operator=(const ThermoElasticFFT&) = delete;                           // Operator '=' Not Implemented
};

#endif /* _ThermoElasticFFT_h_ */
