#include "ThermoElasticFFT.h"

#include <array>
#include <assert.h>
#include <cfenv>
#include <random>
#include <unordered_map>

#include <QtCore/QDateTime>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QJsonArray>
#include <QtCore/QJsonDocument>

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Utilities/SIMPLibEndian.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/RectGridGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Math/SIMPLibRandom.h"
#include "SIMPLib/Math/MatrixMath.h"

#include "OrientationLib/OrientationMath/OrientationTransforms.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#define BUF_SIZE 1024

// Include the MOC generated file for this class
#include "moc_ThermoElasticFFT.cpp"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ThermoElasticFFT::ThermoElasticFFT()
: m_ModuliFile("")
, m_SetupFile("")
, m_FeatureIdsArrayPath("", "", "")
, m_CellPhasesArrayPath("", "", "")
, m_EulerAnglesArrayPath("", "", "")
, m_StressArrayName("")
, m_StrainArrayName("")
, m_FeatureIds(nullptr)
, m_CellPhases(nullptr)
, m_EulerAngles(nullptr)
, m_Stress(nullptr)
, m_Strain(nullptr)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ThermoElasticFFT::~ThermoElasticFFT() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Setup File", SetupFile, FilterParameter::Parameter, ThermoElasticFFT, "*.txt", "Text"));
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Elastic Moduli File", ModuliFile, FilterParameter::Parameter, ThermoElasticFFT, "*.txt", "Text"));

  QStringList linkedProps("EigenStressArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Save Eigen Stresses", SaveEigenStresses, FilterParameter::Parameter, ThermoElasticFFT, linkedProps));
  linkedProps.clear();

  parameters.push_back(SeparatorFilterParameter::New("Element Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Category::Element);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, ThermoElasticFFT, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Category::Element);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phases", CellPhasesArrayPath, FilterParameter::RequiredArray, ThermoElasticFFT, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::TypeNames::Float, 3, AttributeMatrix::Category::Element);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Euler Angles", EulerAnglesArrayPath, FilterParameter::RequiredArray, ThermoElasticFFT, req));
  }

  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Stress Data", StressArrayName, FilterParameter::CreatedArray, ThermoElasticFFT));
  parameters.push_back(SIMPL_NEW_STRING_FP("Strain Data", StrainArrayName, FilterParameter::CreatedArray, ThermoElasticFFT));
  parameters.push_back(SIMPL_NEW_STRING_FP("Eigen Stress Data", EigenStressArrayName, FilterParameter::CreatedArray, ThermoElasticFFT));

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSetupFile(reader->readString("SetupFile", getSetupFile()));
  setModuliFile(reader->readString("ModuliFile", getModuliFile()));
  setEulerAnglesArrayPath(reader->readDataArrayPath("EulerAnglesArrayPath", getEulerAnglesArrayPath()));
  setCellPhasesArrayPath(reader->readDataArrayPath("CellPhasesArrayPath", getCellPhasesArrayPath()));
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setSaveEigenStresses(reader->readValue("SaveEigenStresses", getSaveEigenStresses()));
  setEigenStressArrayName(reader->readString("EigenStressArrayName", getEigenStressArrayName()));
  setStressArrayName(reader->readString("StressArrayName", getStressArrayName()));
  setStrainArrayName(reader->readString("StrainArrayName", getStrainArrayName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();
  DataArrayPath tempPath;
  QVector<IDataArray::Pointer> dataArrays;
  QVector<DataArrayPath> dataArrayPaths;

  QFileInfo fi(getModuliFile());
  if(getModuliFile().isEmpty() == true)
  {
    QString ss = QObject::tr("The moduli file must be set");
    setErrorCondition(-387, ss);
  }
  else if(fi.exists() == false)
  {
    QString ss = QObject::tr("The moduli file does not exist");
    setErrorCondition(-388, ss);
  }

  QFileInfo fi2(getSetupFile());
  if(getSetupFile().isEmpty() == true)
  {
    QString ss = QObject::tr("The setup file must be set");
    setErrorCondition(-387, ss);
  }
  else if(fi2.exists() == false)
  {
    QString ss = QObject::tr("The setup file does not exist");
    setErrorCondition(-388, ss);
  }

  QVector<size_t> cDims(1, 1);
  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getFeatureIdsArrayPath(),
                                                                                                        cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_FeatureIdsPtr.lock())                                                                         /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getFeatureIdsArrayPath());
  }

  m_CellPhasesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getCellPhasesArrayPath(),
                                                                                                        cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_CellPhasesPtr.lock())                                                                         /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getCellPhasesArrayPath());
  }

  cDims[0] = 3;
  m_EulerAnglesPtr =
      getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>, AbstractFilter>(this, getEulerAnglesArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_EulerAnglesPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_EulerAngles = m_EulerAnglesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getEulerAnglesArrayPath());
  }
  tempPath.update(getEulerAnglesArrayPath().getDataContainerName(), getEulerAnglesArrayPath().getAttributeMatrixName(), getStressArrayName());
  cDims[0] = 6;
  m_StressPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0, cDims);             /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_StressPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Stress = m_StressPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }
  tempPath.update(getEulerAnglesArrayPath().getDataContainerName(), getEulerAnglesArrayPath().getAttributeMatrixName(), getStrainArrayName());
  cDims[0] = 6;
  m_StrainPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0,
                                                                                                                     cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_StrainPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Strain = m_StrainPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }

  if (m_SaveEigenStresses)
  {
    tempPath.update(getEulerAnglesArrayPath().getDataContainerName(), getEulerAnglesArrayPath().getAttributeMatrixName(), getEigenStressArrayName());
    cDims[0] = 9;
    m_EigenStressPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0,
                                                                                                                 cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_EigenStressPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_EigenStress = m_EigenStressPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() < 0)
    {
      return;
    }
  }
  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::fourn(std::vector<float>& data, int32_t* nn, int32_t ndim, int32_t isign)
{
  int32_t i1, i2, i2rev, i3, i3rev, idim, ibit, ifp1, ifp2, ip1, ip2, ip3;
  int32_t k1, k2, n, nprev, nrem, ntot;
  float tempi, tempr;
  double theta, wi, wpi, wpr, wr, wtemp;
  ntot = 1;
  for(idim = 0; idim < ndim; idim++)
  {
    ntot *= nn[idim];
  }
  nprev = 1;
  for(idim = 0; idim < ndim; idim++)
  {
    n = nn[idim];
    nrem = ntot / (n * nprev);
    ip1 = 2 * nprev;
    ip2 = ip1 * n;
    ip3 = ip2 * nrem;
    i2rev = 1; // j in four1
    for(i2 = 1; i2 < ip2; i2 += ip1) // i loop in four1
    {
      if(i2 < i2rev)
      {
        for(i1 = i2; i1 <= (i2 + ip1 - 2); i1 += 2)
        {
          for(i3 = i1; i3 < ip3; i3 += ip2)
          {
            i3rev = i2rev + i3 - i2;
            std::swap(data[i3 - 1], data[i3rev - 1]);
            std::swap(data[i3], data[i3rev]);
          }
        }
      }
      ibit = ip2 / 2; // m in four1
      while((ibit >= ip1) && (i2rev > ibit))
      {
        i2rev = i2rev - ibit;
        ibit = ibit / 2;
      }
      i2rev = i2rev + ibit;
    }
    ifp1 = ip1;
    while(ifp1 < ip2)
    {
      ifp2 = 2 * ifp1;
      theta = isign * SIMPLib::Constants::k_2Pi / (ifp2 / ip1);
      wpr = -2.0f * sinf(0.5f * theta) * sinf(0.5f * theta);
      wpi = sinf(theta);
      wr = 1.0f;
      wi = 0.0f;
      for(i3 = 1; i3 < ifp1; i3 += ip1)
      {
        for(i1 = i3; i1 <= (i3 + ip1 - 2); i1 += 2)
        {
          for(i2 = i1; i2 < ip3; i2 += ifp2)
          {
            k1 = i2;
            k2 = k1 + ifp1;
            tempr = static_cast<float>(wr) * data[k2 - 1] - static_cast<float>(wi) * data[k2];
            tempi = static_cast<float>(wr) * data[k2] + static_cast<float>(wi) * data[k2 - 1];
            data[k2 - 1] = data[k1 - 1] - tempr;
            data[k2] = data[k1] - tempi;
            data[k1 - 1] += tempr;
            data[k1] += tempi;
          }
        }
        wtemp = wr;
        wr = wr * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
      }
      ifp1 = ifp2;
    }
    nprev = n * nprev;
  }
  return;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::voigt(float c2[6][6], float c4[3][3][3][3], int iopt)
{
  // *** TRANSFORMS SECOND ORDER MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF
  // *** IOPT=1 AND VICEVERSA IF IOPT=2. IF IOPT=3,TRANSFORMS WITH INV.FACT.
  // *** IOPT=4 FOR GO FROM 6x6 TO 3x3x3x3 WITH Aijkl ANTISYMMETRY

  float F[6][6];
  int32_t IJV[6][2];
  IJV[0][0] = 0;
  IJV[1][0] = 1;
  IJV[2][0] = 2;
  IJV[3][0] = 1;
  IJV[4][0] = 0;
  IJV[5][0] = 0;
  IJV[0][1] = 0;
  IJV[1][1] = 1;
  IJV[2][1] = 2;
  IJV[3][1] = 2;
  IJV[4][1] = 2;
  IJV[5][1] = 1;

  if(iopt == 1)
  {
    for(int32_t i = 0; i < 6; i++)
    {
      int32_t I1 = IJV[i][0];
      int32_t I2 = IJV[i][1];
      for(int32_t j = 0; j < 6; j++)
      {
        int32_t J1 = IJV[j][0];
        int32_t J2 = IJV[j][1];
        c4[I1][I2][J1][J2] = c2[i][j];
        c4[I2][I1][J1][J2] = c2[i][j];
        c4[I1][I2][J2][J1] = c2[i][j];
        c4[I2][I1][J2][J1] = c2[i][j];
      }
    }
  }

  if(iopt == 2)
  {
    for(int32_t i = 0; i < 6; i++)
    {
      int32_t I1 = IJV[i][0];
      int32_t I2 = IJV[i][1];
      for(int32_t j = 0; j < 6; j++)
      {
        int32_t J1 = IJV[j][0];
        int32_t J2 = IJV[j][1];
        c2[i][j] = c4[i][I2][J1][J2];
      }
    }
  }

  if(iopt == 3)
  {
    for(int32_t i = 0; i < 6; i++)
    {
      for(int32_t j = 0; j < 6; j++)
      {
        F[i][j] = 1.0f;
        if(i > 3)
        {
          F[i][j] = 0.5;
        }
        if(j > 3)
        {
          F[i][j] = 0.5 * F[i][j];
        }
      }
    }

    for(int32_t i = 0; i < 6; i++)
    {
      int32_t I1 = IJV[i][0];
      int32_t I2 = IJV[i][1];
      for(int32_t j = 0; j < 6; j++)
      {
        int32_t J1 = IJV[j][0];
        int32_t J2 = IJV[j][1];
        c4[I1][I2][J1][J2] = F[i][j] * c2[i][j];
        c4[I2][I1][J1][J2] = F[i][j] * c2[i][j];
        c4[I1][I2][J2][J1] = F[i][j] * c2[i][j];
        c4[I2][I1][J2][J1] = F[i][j] * c2[i][j];
      }
    }
  }

  if(iopt == 4)
  {
    for(int32_t i = 0; i < 6; i++)
    {
      int32_t I1 = IJV[i][0];
      int32_t I2 = IJV[i][1];
      for(int32_t j = 0; j < 6; j++)
      {
        int32_t J1 = IJV[j][0];
        int32_t J2 = IJV[j][1];
        if(i < 3)
        {
          c4[I1][I2][J1][J2] = c2[i][j];
          c4[I2][I1][J1][J2] = c2[i][j];
          c4[I1][I2][J2][J1] = c2[i][j];
          c4[I2][I1][J2][J1] = c2[i][j];
        }
        else
        {
          c4[I1][I2][J1][J2] = c2[i][j];
          c4[I2][I1][J1][J2] = -c2[i][j];
          c4[I1][I2][J2][J1] = c2[i][j];
          c4[I2][I1][J2][J1] = -c2[i][j];
        }
      }
    }
  }

  return;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::input()
{
  float aux66[6][6], aux3333[3][3][3][3];

  m_InStream1 = fopen(getSetupFile().toLatin1().data(), "r");
  if(m_InStream1 == nullptr)
  {
    QString ss = QObject::tr("Error opening setup file '%1'").arg(getSetupFile());
    setErrorCondition(-48030, ss);
    return;
  }

  //*********INITIALIZATION BLOCK * **************************
  //     CALCULATES TENSORS OF THE SYMMETRIC BASIS 'B(3,3,6)'

  data_crystal();

  //   note that it is essential to call DATA_GRAIN after DATA_CRYSTAL
  data_grain();

  // ****************************************************************************
  // Read RVE DIMENSIONS
  // ****************************************************************************
  char buf[BUF_SIZE];
  // Read Line #2 and dump it
  ::memset(buf, 0, BUF_SIZE);
  fgets(buf, BUF_SIZE, m_InStream1);
  fscanf(m_InStream1, "%f %f %f\n", &delt[0], &delt[1], &delt[2]);

  // ****************************************************************************
  //     READ BOUNDARY CONDITIONS ON OVERALL STRAIN
  // ****************************************************************************
  ::memset(buf, 0, BUF_SIZE);
  fgets(buf, BUF_SIZE, m_InStream1);
  for(int32_t i = 0; i < 3; i++)
  {
    fscanf(m_InStream1, "%f %f %f\n", &udot[i][0], &udot[i][1], &udot[i][2]);
  }

  //     SYMMETRIC STRAIN, ANTISYMMETRIC ROTATION TENSORS
  //     AND INDICES OF IMPOSED COMPONENTS

  for(int32_t i = 0; i < 6; i++)
  {
    for(int32_t j = 0; j < 6; j++)
    {
      dsim[i][j] = (udot[i][j] + udot[j][i]) / 2.0f;
    }
  }

  //     WRITES STRAIN DSIM(I,J) IN b-BASIS AS A 5-DIM VECTOR DBAR(K)

  chg_basis(dbar, dsim, aux66, aux3333, 2);

  // ****************************************************************************
  //     READ ICTRL, THERMCTRL, NSTEPS, ERROR, ITMAX, DELTAT
  // ****************************************************************************
  ::memset(buf, 0, BUF_SIZE);
  fgets(buf, BUF_SIZE, m_InStream1);
  fscanf(m_InStream1, "%i\n", &ictrl);
  fscanf(m_InStream1, "%i\n", &thermctrl);
  fscanf(m_InStream1, "%i\n", &nsteps);
  fscanf(m_InStream1, "%f\n", &error);
  fscanf(m_InStream1, "%i\n", &itmax);
  fscanf(m_InStream1, "%f\n", &deltat);

  if(ictrl <= 3)
  {
    ictrl1 = ictrl-1;
    ictrl2 = ictrl-1;
  }
  else if(ictrl == 4)
  {
    ictrl1 = 1;
    ictrl2 = 2;
  }
  else if(ictrl == 5)
  {
    ictrl1 = 0;
    ictrl2 = 2;
  }
  else if(ictrl == 6)
  {
    ictrl1 = 0;
    ictrl2 = 1;
  }

  chg_basis(ctebar, ctexlsec, aux66, aux3333, 2);

  for(int32_t i = 0; i < 6; i++)
  {
    ctebar[i] = ctebar[i] * deltat;
  }
  int64_t count = 0;
  for(int32_t k = 0; k < npts3; k++)
  {
    for(int32_t j = 0; j < npts2; j++)
    {
      for(int32_t i = 0; i < npts1; i++)
      {
        for(int32_t ii = 0; ii < 6; ii++)
        {
          cteloc6[ii][count] *= deltat;
          //  make CTELOC6 the fluctuation CTE strain only - Sep 09
          //   as of Mar 10, keep the whole CTE eigenstrain (but subtract out ave from DBAR)
        }
        count++;
      }
    }
  }

  dref = dsim[ictrl1][ictrl2] + (ctemat[1][thermctrl-1][thermctrl-1] * deltat);

  //  we'll let dbar change with CTE

  return;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::data_crystal()
{
  float dde[3][3], xid4[3][3][3][3];
  float cmat33local[3][3][3][3];

  for(int32_t i = 0; i < 3; i++)
  {
    for(int32_t j = 0; j < 3; j++)
    {
      dde[i][j] = 0.0f;
      if(i == j)
      {
        dde[i][j] = 1.0f;
      }
    }
  }

  for(int32_t i = 0; i < 3; i++)
  {
    for(int32_t j = 0; j < 3; j++)
    {
      for(int32_t k = 0; k < 3; k++)
      {
        for(int32_t l = 0; l < 3; l++)
        {
          xid4[i][j][k][l] = (dde[i][k] * dde[j][l] + dde[i][l] * dde[j][k]) / 2.0f;
        }
      }
    }
  }

  m_InStream2 = fopen(getModuliFile().toLatin1().data(), "r");
  if(m_InStream2 == nullptr)
  {
    QString ss = QObject::tr("Error opening moduli file '%1'").arg(getModuliFile());
    setErrorCondition(-48030, ss);
    return;
  }

  int32_t nPhases = 0;
  fscanf(m_InStream2, "%i\n", &nPhases);
  // size arrays associted with phases
  nPhases += 1;
  ctemat.resize(nPhases);
  cmat33.resize(nPhases);
  for (int32_t i = 1; i < nPhases; i++)
  {
    ctemat[i].resize(3);
    cmat33[i].resize(3);
    for (int32_t j = 0; j < 3; j++)
    {
      ctemat[i][j].resize(3, 0.0f);
      cmat33[i][j].resize(3);
      for(int32_t k = 0; k < 3; k++)
      {
        cmat33[i][j][k].resize(3);
        for(int32_t l = 0; l < 3; l++)
        {
          cmat33[i][j][k][l].resize(3, 0.0f);
        }
      }
    }
  }

  for(int32_t iter = 1; iter < nPhases; iter++)
  {
    float fact = 0.5f;
    int32_t iso;

    fscanf(m_InStream2, "%i\n", &iso);

    if(iso == 0)
    {
      for(int32_t i = 0; i < 6; i++)
      {
        fscanf(m_InStream2, "%f %f %f %f %f %f\n", &cmat[i][0], &cmat[i][1], &cmat[i][2], &cmat[i][3], &cmat[i][4], &cmat[i][5]);
      }
      voigt(cmat, cmat33local, 1);
      for(int32_t i = 0; i < 3; i++)
      {
        for(int32_t j = 0; j < 3; j++)
        {
          for(int32_t k = 0; k < 3; k++)
          {
            for(int32_t l = 0; l < 3; l++)
            {
              cmat33[iter][i][j][k][l] = cmat33local[i][j][k][l];
            }
          }
        }
      }

      for(int32_t i = 0; i < 3; i++)
      {
        fscanf(m_InStream2, "%f %f %f\n", &ctemat[iter][i][0], &ctemat[iter][i][1], &ctemat[iter][i][2]);
      }

      for(int32_t i = 0; i < 3; i++)
      {
        for(int32_t j = 0; j < 3; j++)
        {
          ctemat[iter][i][j] *= 0.000001;
          // CTE coefficients assumed to be in units of 10 ^ -6
        }
      }
    }
    else
    {
      float young, tmu, tnu, tla;
      fscanf(m_InStream2, "%f %f\n", &young, &tnu);
      // read(UR1,*) young,tnu;
      tmu = young / (2.0f * (1. + tnu));
      tla = 2.0f * tmu * tnu / (1.0f - 2.0f * tnu);

      for(int32_t i = 0; i < 3; i++)
      {
        for(int32_t j = 0; j < 3; j++)
        {
          for(int32_t k = 0; k < 3; k++)
          {
            for(int32_t l = 0; l < 3; l++)
            {
              cmat33[iter][i][j][k][l] = tla * dde[i][j] * dde[k][l] + 2.0f * tmu * xid4[i][j][k][l];
            }
          }
        }
      }
    }
  }

  fclose(m_InStream2);

  return;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::data_grain()
{
  float aa[3][3];
  float caux33[3][3][3][3];
  float caux[6][6];
  std::vector<std::vector<float>> saux;
  std::vector<std::vector<float>> taux;
  float aux6[6];
  float aux33[3][3];
  float aux66[6][6], aux3333[3][3][3][3];

  Eigen::Matrix<float, 6, 6> xlSecTmp, xlSecInvTmp;
  Eigen::Matrix<float, 6, 6> sauxTmp, sauxInvTmp;
  Eigen::Matrix<float, 6, 6> tauxTmp, tauxInvTmp;
  FOrientArrayType om(9, 0.0);

  saux.resize(6);
  taux.resize(6);
  xlsecinv.resize(6);
  for(int32_t iter = 0; iter < 6; iter++)
  {
    saux[iter].resize(6, 0.0f);
    taux[iter].resize(6, 0.0f);
    xlsecinv[iter].resize(6, 0.0f);
  }

  for(int32_t i = 0; i < 6; i++)
  {
    for(int32_t j = 0; j < 6; j++)
    {
      xlsec[i][j] = 0.0f;
    }
  }

  for(int32_t i = 0; i < 3; i++)
  {
    for(int32_t j = 0; j < 3; j++)
    {
      ctexlsec[i][j] = 0.0f;
    }
  }

  int64_t totalpts = (npts1 * npts2 * npts3);
  for(int32_t iter = 0; iter < totalpts; iter++)
  {
    int32_t feature = m_FeatureIds[iter];
    int32_t phase = m_CellPhases[iter];

    //     CALCULATES THE TRANSFORMATION
    //     MATRIX AA, WHICH TRANSFORMS FROM SAMPLE TO CRYSTAL.
    //     AG TRANSFORMS FROM CRYSTAL TO SAMPLE.

    FOrientTransformsType::eu2om(FOrientArrayType(&(m_EulerAngles[3 * iter]), 3), om);
    om.toGMatrix(aa);

    float dum = 0.0f;
    if(phase != 0)
    {
      for(int32_t i1 = 0; i1 < 3; i1++)
      {
        for(int32_t j1 = 0; j1 < 3; j1++)
        {
          for(int32_t k1 = 0; k1 < 3; k1++)
          {
            for(int32_t l1 = 0; l1 < 3; l1++)
            {
              dum = 0.0f;
              for(int32_t i2 = 0; i2 < 3; i2++)
              {
                for(int32_t j2 = 0; j2 < 3; j2++)
                {
                  for(int32_t k2 = 0; k2 < 3; k2++)
                  {
                    for(int32_t l2 = 0; l2 < 3; l2++)
                    {
                      dum += (aa[i2][i1] * aa[j2][j1] * aa[k2][k1] * aa[l2][l1] * cmat33[phase][i2][j2][k2][l2]);
                    }
                  }
                }
              }
              caux33[i1][j1][k1][l1] = dum;
            }
          }
        }
      }
    }
    else
    {
      for(int32_t i1 = 0; i1 < 3; i1++)
      {
        for(int32_t j1 = 0; j1 < 3; j1++)
        {
          for(int32_t k1 = 0; k1 < 3; k1++)
          {
            for(int32_t l1 = 0; l1 < 3; l1++)
            {
              caux33[i1][j1][k1][l1] = 0.0f;
            }
          }
        }
      }
    }
    //  CAUX33 = 4th rank C0
    chg_basis(aux6, aux33, caux, caux33, 4);
    //  converts C0 to 3x3
    //xlsec = Co - average over all points
    //cloc = pointwise stiffness
    for(int32_t i = 0; i < 6; i++)
    {
      for(int32_t j = 0; j < 6; j++)
      {
        cloc[i][j][iter] = caux[i][j];
        xlsec[i][j] += (caux[i][j] * wgt);
      }
    }

    //     Calculations for CTEs
    float dum2 = 0.0f;
    if(phase != 0)
    {
      for(int32_t i1 = 0; i1 < 3; i1++)
      {
        for(int32_t j1 = 0; j1 < 3; j1++)
        {
          dum2 = 0.0f;
          for(int32_t i2 = 0; i2 < 3; i2++)
          {
            for(int32_t j2 = 0; j2 < 3; j2++)
            {
              dum2 += (aa[i2][i1] * aa[j2][j1] * ctemat[phase][i2][j2]);
            }
          }
          cteaux33[i1][j1] = dum2;
        }
      }
    }
    else
    {
      for(int32_t i1 = 0; i1 < 3; i1++)
      {
        for(int32_t j1 = 0; j1 < 3; j1++)
        {
          cteaux33[i1][j1] = 0.0f;
        }
      }
    }

    chg_basis(aux6, cteaux33, aux66, aux3333, 2);
    //  converts to 1x6 vector

    //cteloc = pointwise cte coefficients - alphaIJ
    for(int32_t i = 0; i < 6; i++)
    {
      cteloc6[i][iter] = aux6[i];
    }

    // at this point, the CTE strain is the full eigenstrain at each point
    //  CTEaux = 2nd rank strain tensor
    //ctexlsec = cte coefficients averaged over all points
    for(int32_t i = 0; i < 3; i++)
    {
      for(int32_t j = 0; j < 3; j++)
      {
        ctexlsec[i][j] += (cteaux33[i][j] * wgt);
      }
    }
    //  CTExlsec is weighted version of CTE strain
  } // Ending the kkk = npts1 * npts2 * npts3 loop

  for(int32_t p = 0; p < 6; p++)
  {
    for(int32_t q = 0; q < 6; q++)
    {
      xlSecTmp(p, q) = xlsec[p][q];
    }
  }
  xlSecInvTmp = xlSecTmp.inverse();
  for(int32_t p = 0; p < 3; p++)
  {
    for(int32_t q = 0; q < 3; q++)
    {
      xlsecinv[p][q] = xlSecInvTmp(p, q);
    }
  }

  chg_basis(aux6, aux33, xlsec, xlsec33, 3);
  // convert to 4th rank
  int64_t count = 0;
  for(int32_t kk = 0; kk < npts3; kk++)
  {
    for(int32_t jj = 0; jj < npts2; jj++)
    {
      for(int32_t ii = 0; ii < npts1; ii++)
      {
        for(int32_t i = 0; i < 6; i++)
        {
          for(int32_t j = 0; j < 6; j++)
          {
            saux[i][j] = cloc[i][j][count];
            sauxTmp(i, j) = saux[i][j];
          }
        }
        sauxInvTmp = sauxTmp.inverse();
        for(int32_t i = 0; i < 6; i++)
        {
          for(int32_t j = 0; j < 6; j++)
          {
            saux[i][j] = sauxInvTmp(i, j);
          }
        }
        // invert CLOC(like XLSEC, unweighted however) !so that SAUX is now C ^-1
        float dum = 0.0f;
        for(int32_t i = 0; i < 6; i++)
        {
          for(int32_t j = 0; j < 6; j++)
          {
            dum = 0.0f;
            for(int32_t k = 0; k < 6; k++)
            {
              dum += (xlsec[i][k] * saux[k][j]);
              // product of C0 and C ^ -1
            }
            if(i == j)
            {
              dum += 1.0f;
              //              !add the identity
            }
            taux[i][j] = dum;
            tauxTmp(i, j) = taux[i][j];
          }
        }
        tauxInvTmp = tauxTmp.inverse();
        for(int32_t i = 0; i < 6; i++)
        {
          for(int32_t j = 0; j < 6; j++)
          {
            taux[i][j] = tauxInvTmp(i, j);
          }
        }

        // and get its inverse
        for(int32_t i = 0; i < 6; i++)
        {
          for(int32_t j = 0; j < 6; j++)
          {
            dum = 0.0f;
            for(int32_t k = 0; k < 6; k++)
            {
              dum += (saux[i][k] * taux[k][j]);
            }
            fsloc[i][j][count] = dum;
            // finally get the product of C ^ -1 with I + (C0: C ^ -1)
          }
        }
        count++;
      }
    }
  }

  return;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::chg_basis(float ce2[6], float c2[3][3], float ce4[6][6], float c4[3][3][3][3], int iopt)
{
  int32_t kDim = 6;
  float SQR2 = SIMPLib::Constants::k_Sqrt2;
  float RSQ2 = SIMPLib::Constants::k_1OverRoot2;
  float RSQ3 = SIMPLib::Constants::k_1OverRoot3;
  float RSQ6 = 1.0f / sqrtf(6.0f);

  float b[3][3][6];

  // CALCULATES BASIS TENSORS B(N)
  for(int32_t k = 0; k < 6; k++)
  {
    for(int32_t j = 0; j < 3; j++)
    {
      for(int32_t i = 0; i < 3; i++)
      {
        b[i][j][k] = 0.0f;
      }
    }
  }

  b[0][0][1] = -RSQ6;
  b[1][1][1] = -RSQ6;
  b[2][2][1] = 2.0f * RSQ6;
  // N = 0 = axisymmetric along 3
  b[0][0][0] = -RSQ2;
  b[1][1][0] = RSQ2;
  //  N = 1  simple shear (110) + [-1 1 0 ]
  b[1][2][2] = RSQ2;
  b[2][1][2] = RSQ2;
  //  N = 2,4 = 001 - 010 shear modes
  b[0][2][3] = RSQ2;
  b[2][0][3] = RSQ2;
  b[0][1][4] = RSQ2;
  b[1][0][4] = RSQ2;

  b[0][0][5] = RSQ3;
  b[1][1][5] = RSQ3;
  b[2][2][5] = RSQ3;
  //  N=5  is the hydrostatic component
  
  if(iopt == 1)
  {
    // CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
    for(int32_t i = 0; i < 3; i++)
    {
      for(int32_t j = 0; j < 3; j++)
      {
        c2[i][j] = 0.0;
        for(int32_t n = 0; n < 6; n++)
        {
          c2[i][j] += (ce2[n] * b[i][j][n]);
        }
      }
    }
  }
  else if(iopt == 2)
  {
    // CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
    for(int32_t n = 0; n < 6; n++)
    {
      ce2[n] = 0.0;
      for(int32_t i = 0; i < 3; i++)
      {
        for(int32_t j = 0; j < 3; j++)
        {
          ce2[n] += (c2[i][j] * b[i][j][n]);
        }
      }
    }
  }
  else if(iopt == 3)
  {
    // CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
    for(int32_t i = 0; i < 3; i++)
    {
      for(int32_t j = 0; j < 3; j++)
      {
        for(int32_t k = 0; k < 3; k++)
        {
          for(int32_t l = 0; l < 3; l++)
          {
            c4[i][j][k][l] = 0.0;
            for(int32_t n = 0; n < 6; n++)
            {
              for(int32_t m = 0; m < 6; m++)
              {
                c4[i][j][k][l] += (ce4[n][m] * b[i][j][n] * b[k][l][m]);
              }
            }
          }
        }
      }
    }
  }
  else if(iopt == 4)
  {
    // CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
    for(int32_t n = 0; n < 6; n++)
    {
      for(int32_t m = 0; m < 6; m++)
      {
        ce4[n][m] = 0.0;
        for(int32_t i = 0; i < 3; i++)
        {
          for(int32_t j = 0; j < 3; j++)
          {
            for(int32_t k = 0; k < 3; k++)
            {
              for(int32_t l = 0; l < 3; l++)
              {
                ce4[n][m] += (c4[i][j][k][l] * b[i][j][n] * b[k][l][m]);
              }
            }
          }
        }
      }
    }
  }
  return;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::augm_lagr()
{
  //      INCLUDE 'elfft.dim'

  float x[6], dg[6], edot[6], dsg[6], ddg[6];

  errd = 0.0f;
  errs = 0.0f;

  int64_t count = 0;
  for(int32_t kk = 0; kk < npts3; kk++)
  {
    for(int32_t jj = 0; jj < npts2; jj++)
    {
      for(int32_t ii = 0; ii < npts1; ii++)
      {
        if(m_CellPhases[count] >= 1)
        {
          //  Sep 09 - need to apply CTE strain to all valid phases
          //     DTILDE comes from FFT^-1 of Gamma:tau
          //      tau is the fluctuation part of the stress field (not the total elastic stress)
          //  thus, DTILDE should be the fluctuation in the (complete) elastic strain,
          //    but NOT the total strain (DBAR)
          for(int32_t i = 0; i < 6; i++)
          {
            dg[i] = dbar[i] + dtilde[i][count] - cteloc6[i][count];
          }
          //     get the local strain, "epsilon", as average, "E",
          //     - the fluctuation term  (changed sign as of 8 viii 09)
          //     Eq 4  -  DG is the actual strain at each point
          //     now we subtract the thermoelastic source term
          //     Sep 09: note that the CTE strain is fixed by the input,
          //          whereas DTILDE changes
          //  Mar 10: should write this as, elastic strain = total strain - eigenstrain

          for(int32_t i = 0; i < 6; i++)
          {
            x[i] = sg[i][count];
            for(int32_t j = 0; j < 6; j++)
            {
              x[i] += (xlsec[i][j] * dg[j]);
            }
          }
          //     get fluctuation stress as previous value + C0 * strain
          //     Eq 5
          //     contribution from thermoelastic term does NOT enter here
          //        because it was subtracted above

          for(int32_t i = 0; i < 6; i++)
          {
            //     EDOT is "e" in the eqs
            edot[i] = 0;
            //     changed from "-" to "+" as of 8 viii 09
            //     also re-positioned the TE strain contribution to avoid
            //     doing it 6 times over inside the J loop (Thanks to Ben!)
            //  Sep09: removed this addition of the CTE strain, since
            //    the local calculated fluctuation strain should include the
            //    CTE contribution

            for(int32_t j = 0; j < 6; j++)
            {
              edot[i] += (fsloc[i][j][count] * x[j]);
            }
          }
          //     get strain, "e", as C^-1:(I+C0:C^-1) * {new stress}
          //     Eq 6
          //     and add on the thermoelastic term

          //     UPDATE SG (LAGRANGE MULTIPLIER)

          for(int32_t i = 0; i < 6; i++)
          {
            ddg[i] = dg[i] - edot[i];
            // change in strain is "epsilon" minus "e"
            dsg[i] = 0.0f;
            for(int32_t j = 0; j < 6; j++)
            {
              dsg[i] += (xlsec[i][j] * (dg[j] - edot[j]));
              //     change in stress = C0 * { change in strain }
              //     Eq 7
              //       Sep 09 - but, only the elastic strain part,
              //           not the entire eigenstrain (CTE)
            }
            sg[i][count] += dsg[i];
          }

          errd += tnorm(ddg, 6, 1) * wgt;
          errs += tnorm(dsg, 6, 1) * wgt;
        }
        count++;
      }
    }
  }

  return;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float ThermoElasticFFT::tnorm(float v[36], int32_t nrows, int32_t ncols)
{
  float tnorm = 0.0f;
  for(int32_t i = 0; i < (nrows * ncols); i++)
  {
    tnorm = tnorm + v[i] * v[i];
  }
  tnorm = sqrtf(tnorm);

  return tnorm;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ThermoElasticFFT::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  ImageGeom::Pointer image = getDataContainerArray()->getDataContainer(m_FeatureIdsArrayPath.getDataContainerName())->getGeometryAs<ImageGeom>();
  SizeVec3Type udims = image->getDimensions();
  size_t npts1 = static_cast<size_t>(udims[0]);
  size_t npts2 = static_cast<size_t>(udims[1]);
  size_t npts3 = static_cast<size_t>(udims[2]);

  int32_t nn[3], nn2[2];
  int32_t kx, ky, kz;
  float xk[3];
  float g1[3][3][3][3];
  std::vector<std::vector<float>> a;
  a.resize(3);
  for (int32_t iter = 0; iter < 3; iter++)
  {
    a[iter].resize(3, 0.0f);
  }

  Eigen::Matrix<float, 3, 3> aMat, aInv;

  int64_t totalpts = npts1 * npts2 * npts3;
  std::vector<float> data(2 * totalpts);
  std::vector<std::vector<float>> delta(6);
  std::vector<std::vector<float>> deltaim(6);
  std::vector<std::vector<std::vector<float>>> cnew(6);
  sg.resize(6);
  dtilde.resize(6);
  cteloc6.resize(6);
  velgrad.resize(6);
  velgradim.resize(6);
  cloc.resize(6);
  fsloc.resize(6);
  for(int32_t iter = 0; iter < 6; iter++)
  {
    delta[iter].resize(totalpts, 0.0f);
    deltaim[iter].resize(totalpts, 0.0f);
    cnew[iter].resize(6);
    sg[iter].resize(totalpts, 0.0f);
    dtilde[iter].resize(totalpts, 0.0f);
    cteloc6[iter].resize(totalpts, 0.0f);
    velgrad[iter].resize(6);
    velgradim[iter].resize(6);
    cloc[iter].resize(6);
    fsloc[iter].resize(6);
    for(int32_t iter2 = 0; iter2 < 6; iter2++)
    {
      cnew[iter][iter2].resize(totalpts, 0.0f);
      velgrad[iter][iter2].resize(totalpts, 0.0f);
      velgradim[iter][iter2].resize(totalpts, 0.0f);
      cloc[iter][iter2].resize(totalpts, 0.0f);
      fsloc[iter][iter2].resize(totalpts, 0.0f);
    }
  }
  float aux6[6], aux33[3][3];
  float aux66[6][6], aux3333[3][3][3][3];

  float d6[6], d6im[6], d33[3][3], d33im[3][3];

  float daux[6], dloc[3][3];
  float saux[6], sloc[3][3];

  float sdevstrain[6], sdevstress[6];
  float avestrain[6]; // to include thermoelastic strain
  float avetotstrain[6];
  float aveCTEstrain[6]; // estimated thermoelastic strain over polyxtal !Sep 09 - but, this is same as CTEBAR

  float sbarMax, dbarMax, dum;

  nn[0] = npts1;
  nn[1] = npts2;
  nn[2] = npts3;

  nn2[0] = npts1;
  nn2[1] = npts2;

  int32_t prodnn = float(totalpts);
  wgt = 1.0f / float(prodnn);

  //setup structure, moduli, boundary conditions, etc
  input();

  //going to use this throughout as a counter in triple nested loops for indexing into flattened out arrays
  int64_t count = 0;

  for(int32_t ii = 0; ii < 6; ii++)
  {
    sbar[ii] = 0.0f;
    dtilbar[ii] = 0.0f;
    sdevstrain[ii] = 0.0f;
    sdevstress[ii] = 0.0f;
    aveCTEstrain[ii] = 0.0f;
  }

  // inserted extra loop to get initial CTE eigenstrain,2 iii 10
  for(int32_t k = 0; k < npts3; k++)
  {
    for(int32_t j = 0; j < npts2; j++)
    {
      for(int32_t i = 0; i < npts1; i++)
      {
        // INITIALIZE D ~, SG
        for(int32_t ii = 0; ii < 6; ii++)
        {
          dtilde[ii][count] = 0.0f;
        }

        for(int32_t ii = 0; ii < 6; ii++)
        {
          aveCTEstrain[ii] += (cteloc6[ii][count] * wgt);
        }
        count++;
      }
    }
  }

  for(int32_t i = 0; i < 6; i++)
  {
    dbar[i] += aveCTEstrain[i];
  }
  // and here we simply add the ave CTE eigenstrain off the overall strain

  count = 0;
  for(int32_t k = 0; k < npts3; k++)
  {
    for(int32_t j = 0; j < npts2; j++)
    {
      for(int32_t i = 0; i < npts1; i++)
      {
        for(int32_t ii = 0; ii < 6; ii++)
        {
          if(m_CellPhases[count] == 0)
          {
            // viii 09, changed this BACK to phase = 0, ADR
            sg[ii][count] = 0.0f;
            // clearly, this had to do with composites with void space(carbon foam)
          }
          else
          {
            for(int32_t jj = 0; jj < 6; jj++)
            {
              sg[ii][count] += (cloc[ii][jj][count] * (dbar[jj] + dtilde[jj][count] - cteloc6[jj][count]));
              if(m_SaveEigenStresses)
              {
                m_EigenStress[9 * count + ii] += (cloc[ii][jj][count] * fabsf(dbar[jj] + dtilde[jj][count] - cteloc6[jj][count]));
              }
            }
            // Mar 10 - initial estimate of elastic stress includes !CTE eigenstrain(as local fluctuation) !i.e.elastic stress = total strain - eigenstrain
          }

          sbar[ii] += sg[ii][count] * wgt;
        }
        count++;
      }
    }
  }

  chg_basis(sbar, stens, aux66, aux3333, 1);

  sref = stens[ictrl1][ictrl2];

  for(int32_t imicro = 0; imicro < nsteps; imicro++)
  {
    int32_t iter = 0;
    float err2mod = 2 * error;
    while(iter < itmax && err2mod > error)
    {
      iter++;
      for(int32_t ii = 0; ii < 6; ii++)
      {
        int32_t k1 = 0;
        count = 0;
        for(int32_t k = 0; k < npts3; k++)
        {
          for(int32_t j = 0; j < npts2; j++)
          {
            for(int32_t i = 0; i < npts1; i++)
            {
              delta[ii][count] = sg[ii][count];
              for(int32_t jj = 0; jj < 6; jj++)
              {
                delta[ii][count] -= (xlsec[ii][jj] * dtilde[jj][count]);
                // as in the ordinary Elastic FFT !$ - xlsec(ii, jj) * (dtilde(jj, i, j, k) + cteloc6(jj, i, j, k)) !Sep 09 - DTILDE varies,
                // but CTELOC6 is(local) fixed CTE strain !Sep 09 - DTILDE at this point is ONLY the fluctuation strain but !NOT including the CTE strain, obviously !so,
                // we have to add it back in to estimate the stress !Mar 10 : we do not do anything with CTE strain at this point

                // Added to see if this would change the stress error - BSA 4 / 2 / 09 !delta(ii, i, j, k) = delta(ii, i, j, k) - xlsec(ii, jj) *
                //(dtilde(jj, i, j, k) - cteloc6(jj, i, j, k) * deltat)

                // Necessary according to Vinogradov algorithm !cnew(ii, jj, i, j, k) = cloc(ii, jj, i, j, k) - xlsec(ii, jj) !delta(ii, i, j, k) = delta(ii, i, j, k)
                // 1 + cnew(ii, jj, i, j, k) * dtilde(jj, i, j, k) !removed, as of 19 vii 09, ADR
              }

              data[k1] = delta[ii][count];
              k1++;
              data[k1] = 0.0f;
              k1++;
              count++;
            }
          }
        }
        if(npts3 > 1)
        {
          fourn(data, nn, 3, 1);
        }
        else
        {
          fourn(data, nn2, 2, 1);
        }
        k1 = 0;
        count = 0;
        for(int32_t kzz = 0; kzz < npts3; kzz++)
        {
          for(int32_t kyy = 0; kyy < npts2; kyy++)
          {
            for(int32_t kxx = 0; kxx < npts1; kxx++)
            {
              delta[ii][count] = data[k1];
              k1++;
              deltaim[ii][count] = data[k1];
              k1++;
              count++;
            }
          }
        }
      }
      count = 0;
      for(int32_t kzz = 0; kzz < npts3; kzz++)
      {
        for(int32_t kyy = 0; kyy < npts2; kyy++)
        {
          for(int32_t kxx = 0; kxx < npts1; kxx++)
          {
            for(int32_t i = 0; i < 6; i++)
            {
              d6[i] = delta[i][count];
              d6im[i] = deltaim[i][count];
            }
            chg_basis(d6, d33, aux66, aux3333, 1);
            chg_basis(d6im, d33im, aux66, aux3333, 1);
            if(kxx <= npts1 / 2)
            {
              kx = kxx;
            }
            if(kxx > npts1 / 2)
            {
              kx = kxx - npts1;
            }
            if(kyy <= npts2 / 2)
            {
              ky = kyy;
            }
            if(kyy > npts2 / 2)
            {
              ky = kyy - npts2;
            }
            if(kzz <= npts3 / 2)
            {
              kz = kzz;
            }
            if(kzz > npts3 / 2)
            {
              kz = kzz - npts3;
            }

            xk[0] = kx / (delt[0] * nn[0]);
            xk[1] = ky / (delt[1] * nn[1]);

            if(npts3 > 1)
            {
              xk[2] = kz / (delt[2] * nn[2]);
            }
            else
            {
              xk[2] = 0.0f;
            }

            float xknorm = sqrtf((xk[0] * xk[0]) + (xk[1] * xk[1]) + (xk[2] * xk[2]));
            if(xknorm != 0.0f)
            {
              xk[0] /= xknorm;
              xk[1] /= xknorm;
              xk[2] /= xknorm;
            }

            for(int32_t i = 0; i < 3; i++)
            {
              for(int32_t k = 0; k < 3; k++)
              {
                a[i][k] = 0.0f;
                for(int32_t j = 0; j < 3; j++)
                {
                  for(int32_t l = 0; l < 3; l++)
                  {
                    a[i][k] += (xlsec33[i][j][k][l] * xk[j] * xk[l]);
                  }
                }
                aMat(i, k) = a[i][k];
              }
            }

            aInv = aMat.inverse();
            for(int32_t p = 0; p < 3; p++)
            {
              for(int32_t q = 0; q < 3; q++)
              {
                a[p][q] = aInv(p, q);
              }
            }


            for(int32_t p = 0; p < 3; p++)
            {
              for(int32_t q = 0; q < 3; q++)
              {
                for(int32_t i = 0; i < 3; i++)
                {
                  for(int32_t j = 0; j < 3; j++)
                  {
                    g1[p][q][i][j] = -a[p][i] * xk[q] * xk[j];
                  }
                }
              }
            }

            //                   !Apply De Graef filter if(iq_filter /= 0)
            //                       then do p = 1,
            // 3 do q = 1, 3 do i = 1, 3 do j = 1,
            // 3 if(iq_filter == 1.AND.xknorm /= 0) then g1(p, q, i, j) =
            //    g1(p, q, i, j) &
            //    *(3 / ((xknorm) * *2) & *((sin(xknorm) / (xknorm)) * *2 & -sin(2 * xknorm) / (2 * xknorm)) & +alpha & !ALPHA value * 2 / (xknorm * *2) &
            //      *(3 * (sin(2 * xknorm) / (2 * xknorm)) * *2 & -2 * sin(2 * xknorm) / (2 * xknorm) & -sin(4 * xknorm) / (4 * xknorm))) endif

            //    if(iq_filter == 2.OR.iq_filter == 3 &.OR.iq_filter == 6.OR. iq_filter == 7) then if(xk(1) /= 0) then g1(p, q, i, j) =
            //        g1(p, q, i, j) &
            //        *(3 / ((xk(1)) * *2) & *((sin(xk(1)) / (xk(1))) * *2 & -sin(2 * xk(1)) / (2 * xk(1))) & +alpha & !ALPHA value * 2 / ((xk(1)) * *2) &
            //          *(3 * (sin(2 * xk(1)) / (2 * xk(1))) * *2 & -2 * sin(2 * xk(1)) / (2 * xk(1)) & -sin(4 * xk(1)) / (4 * xk(1)))) endif endif

            //        if(iq_filter == 2.OR.iq_filter == 4 &.OR.iq_filter == 6.OR.iq_filter == 8) then if(xk(2) /= 0) then g1(p, q, i, j) =
            //            g1(p, q, i, j) &
            //            *(3 / ((xk(2)) * *2) & *((sin(xk(2)) / (xk(2))) * *2 & -sin(2 * xk(2)) / (2 * xk(2))) & +alpha & !ALPHA VALUE * 2 / ((xk(2)) * *2) &
            //              *(3 * (sin(2 * xk(2)) / (2 * xk(2))) * *2 & -2 * sin(2 * xk(2)) / (2 * xk(2)) & -sin(4 * xk(1)) / (4 * xk(2)))) endif endif

            //            if(iq_filter == 2.OR.iq_filter == 5 &.OR.iq_filter == 7.OR.iq_filter == 8) then if(xk(3) /= 0) then g1(p, q, i, j) =
            //                g1(p, q, i, j) &
            //                *(3 / ((xk(3)) * *2) & *((sin(xk(3)) / (xk(3))) * *2 & -sin(2 * xk(3)) / (2 * xk(3))) & +alpha & !ALPHA VALUE * 2 / ((xk(3)) * *2) &
            //                  *(3 * (sin(2 * xk(3)) / (2 * xk(3))) * *2 & -2 * sin(2 * xk(3)) / (2 * xk(3)) & -sin(4 * xk(3)) / (4 * xk(3)))) endif endif

            //                if(iq_filter == 9.AND.xknorm /= 0) then g1(p, q, i, j) =
            //                    g1(p, q, i, j) * sin(xknorm) /
            //                    (xknorm)endif

            // enddo enddo enddo enddo endif

            for(int32_t i = 0; i < 3; i++)
            {
              for(int32_t j = 0; j < 3; j++)
              {
                velgrad[i][j][count] = 0.0f;
                velgradim[i][j][count] = 0.0f;

                if(kx != 0 || ky != 0 || kz != 0)
                {
                  for(int32_t k = 0; k < 3; k++)
                  {
                    for(int32_t l = 0; l < 3; l++)
                    {
                      velgrad[i][j][count] += (g1[i][j][k][l] * d33[k][l]);
                      velgradim[i][j][count] += (g1[i][j][k][l] * d33im[k][l]);
                    }
                  }
                }
              }
            }
            count++;
          }
        }
      }
      // which at this point includes both CTE and regular strain

      for(int32_t m = 0; m < 3; m++)
      {
        for(int32_t n = 0; n < 3; n++)
        {
          int32_t k1 = 0;
          count = 0;
          for(int32_t k = 0; k < npts3; k++)
          {
            for(int32_t j = 0; j < npts2; j++)
            {
              for(int32_t i = 0; i < npts1; i++)
              {
                data[k1] = velgrad[m][n][count];
                k1 = k1 + 1;
                data[k1] = velgradim[m][n][count];
                k1 = k1 + 1;
                count++;
              }
            }
          }

          if(npts3 > 1)
          {
            fourn(data, nn, 3, -1);
          }
          else
          {
            fourn(data, nn2, 2, -1);
          }

          for(int32_t i = 0; i < (2 * totalpts); i++)
          {
            data[i] /= prodnn;
          }

          k1 = 0;
          count = 0;
          for(int32_t kzz = 0; kzz < npts3; kzz++)
          {
            for(int32_t kyy = 0; kyy < npts2; kyy++)
            {
              for(int32_t kxx = 0; kxx < npts1; kxx++)
              {
                velgrad[m][n][count] = data[k1];
                k1++;
                k1++;
                count++;
              }
            }
          }
        }
      }
      count = 0;
      for(int32_t kzz = 0; kzz < npts3; kzz++)
      {
        for(int32_t kyy = 0; kyy < npts2; kyy++)
        {
          for(int32_t kxx = 0; kxx < npts1; kxx++)
          {
            for(int32_t ii = 0; ii < 3; ii++)
            {
              for(int32_t jj = 0; jj < 3; jj++)
              {
                aux33[ii][jj] = (velgrad[ii][jj][count] + velgrad[jj][ii][count]) / 2.0f;
              }
            }

            chg_basis(aux6, aux33, aux66, aux3333, 2);

            for(int32_t m = 0; m < 6; m++)
            {
              dtilde[m][count] = aux6[m];
            }
            count++;
          }
        }
      }

      //********Sep 09 - CAREFUL !!!DTILDE at this point has full strain in it
      augm_lagr();
      // !********Sep 09 - however, DTILDE emerging from AUGM_LAGR !is only the regular strain, not including CTE ! % % % % % % % % Mar 10 - DTILDE now is the elastic strain

      for(int32_t ii = 0; ii < 6; ii++)
      {
        sbar[ii] = 0.0f;
        dtilbar[ii] = 0.0f;
        avestrain[ii] = 0.0f;
        avetotstrain[ii] = 0.0f;
        sdevstrain[ii] = 0.0f;
        sdevstress[ii] = 0.0f;
        dbarMax = 0.0f;
      }
      count = 0;
      for(int32_t k = 0; k < npts3; k++)
      {
        for(int32_t j = 0; j < npts2; j++)
        {
          for(int32_t i = 0; i < npts1; i++)
          {
            for(int32_t ii = 0; ii < 6; ii++)
            {
              sbar[ii] += sg[ii][count] * wgt;
              dtilbar[ii] += (dtilde[ii][count] * wgt);
              avestrain[ii] += (dbar[ii] + dtilde[ii][count] - cteloc6[ii][count]) * wgt;
              avetotstrain[ii] += (dbar[ii] + dtilde[ii][count]) * wgt;
              // Sep 09; Mar 10 : CTE strain already included, avestrain is average ELASTIC strain
            }
            count++;
          }
        }
      }

      if(iter < 3)
      {
        sbarMax = 0.0f;
        for(int32_t ii = 0; ii < 6; ii++)
        {
          if(abs(sbar[ii]) > sbarMax)
          {
            sbarMax = abs(sbar[ii]);
            // logic : capture an early value of max stress, and then scale changes to that
          }
        }
      }

      for(int32_t ii = 0; ii < 6; ii++)
      {
        if(abs(dbar[ii]) > dbarMax)
        {
          dbarMax = abs(dbar[ii]);
        }
      }
      count = 0;
      for(int32_t k = 0; k < npts3; k++)
      {
        for(int32_t j = 0; j < npts2; j++)
        {
          for(int32_t i = 0; i < npts1; i++)
          {
            for(int32_t ii = 0; ii < 6; ii++)
            {
              sdevstrain[ii] += ((dtilde[ii][count] - dtilbar[ii]) * (dtilde[ii][count] - dtilbar[ii])) * wgt;
              sdevstress[ii] += ((sg[ii][count] - sbar[ii]) * (sg[ii][count] - sbar[ii])) * wgt;
            }
            count++;
          }
        }
      }

      for(int32_t i = 0; i < 6; i++)
      {
        dum = 0.0f;
        for(int32_t j = 0; j < 6; j++)
        {
          dum -= (xlsecinv[i][j] * sbar[j]);
        }
        dbar[i] = avetotstrain[i] + dum;
      }

      chg_basis(sbar, stens, aux66, aux3333, 1);
      chg_basis(dbar, dtens, aux66, aux3333, 1);

      dref = dtens[ictrl1][ictrl2];
      sref = stens[ictrl1][ictrl2];

      err2mod = fabsf(errs / sref);
//      err2mod = (errs / sref) > (errd / dref) ? (errs / sref) : (errd / dref);

		  QString ss = QObject::tr("Iteration: %1 || Error: %2").arg(iter).arg(err2mod);
      notifyStatusMessage(ss);

    } // ends while iterations
  } // ends microSteps

  count = 0;
  for(int32_t k = 0; k < npts3; k++)
  {
    for(int32_t j = 0; j < npts2; j++)
    {
      for(int32_t i = 0; i < npts1; i++)
      {
        for(int32_t ii = 0; ii < 6; ii++)
        {
          daux[ii] = dbar[ii] + dtilde[ii][count] - cteloc6[ii][count];
          saux[ii] = sg[ii][count];
        }

        chg_basis(daux, dloc, aux66, aux3333, 1);
        daux[0] = dloc[0][0];
        daux[1] = dloc[1][1];
        daux[2] = dloc[2][2];
        daux[3] = dloc[1][2];
        daux[4] = dloc[0][2];
        daux[5] = dloc[0][1];
        chg_basis(saux, sloc, aux66, aux3333, 1);
        saux[0] = sloc[0][0];
        saux[1] = sloc[1][1];
        saux[2] = sloc[2][2];
        saux[3] = sloc[1][2];
        saux[4] = sloc[0][2];
        saux[5] = sloc[0][1];
        
        m_Stress[6 * count + 0] = saux[0];
        m_Stress[6 * count + 1] = saux[1];
        m_Stress[6 * count + 2] = saux[2];
        m_Stress[6 * count + 3] = saux[3];
        m_Stress[6 * count + 4] = saux[4];
        m_Stress[6 * count + 5] = saux[5];
        m_Strain[6 * count + 0] = daux[0];
        m_Strain[6 * count + 1] = daux[1];
        m_Strain[6 * count + 2] = daux[2];
        m_Strain[6 * count + 3] = daux[3];
        m_Strain[6 * count + 4] = daux[4];
        m_Strain[6 * count + 5] = daux[5];
        count++;
      }
    }
  }

  notifyStatusMessage("Complete");
}

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  AbstractFilter::Pointer ThermoElasticFFT::newFilterInstance(bool copyFilterParameters) const
  {
    ThermoElasticFFT::Pointer filter = ThermoElasticFFT::New();
    if(true == copyFilterParameters)
    {
      copyFilterParameterInstanceVariables(filter.get());
    }
    return filter;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  const QString ThermoElasticFFT::getCompiledLibraryName() const
  {
    return DREAM3DReviewConstants::DREAM3DReviewBaseName;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  const QString ThermoElasticFFT::getBrandingString() const
  {
    return "DREAM3DReview";
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  const QString ThermoElasticFFT::getFilterVersion() const
  {
    QString version;
    QTextStream vStream(&version);
    vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
    return version;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  const QString ThermoElasticFFT::getGroupName() const
  {
    return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  const QString ThermoElasticFFT::getSubGroupName() const
  {
    return "SimulationFilters";
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  const QString ThermoElasticFFT::getHumanLabel() const
  {
    return "Thermal Elastic FFT";
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  const QUuid ThermoElasticFFT::getUuid()
  {
    return QUuid("");
  }