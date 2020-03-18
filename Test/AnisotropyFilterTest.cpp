/* ============================================================================
 * Copyright (c) 2016 Czech Academy of Sciences, Institute of Physics,
 * Group of Bulk Nanomaterials and Interfaces, http://ams.fzu.cz
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

#include <QtCore/QDir>
#include <QtCore/QJsonDocument>

#include <QtCore/QDebug>

#include "SIMPLib/DataArrays/DataArray.hpp"

#include "SIMPLib/FilterParameters/FileListInfoFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/ComparisonInputs.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"

#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

class AnisotropyFilterTest
{
public:
  AnisotropyFilterTest() = default;
  ~AnisotropyFilterTest() = default;

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void RemoveTestFiles()
  {
#if REMOVE_TEST_FILES
    // remove output files
    QFile::remove(UnitTest::AnisotropyTest::TestOutput1);
    QFile::remove(UnitTest::AnisotropyTest::TestOutput2);
    QFile::remove(UnitTest::AnisotropyTest::TestOutput3);
    QFile::remove(UnitTest::AnisotropyTest::TestOutput4);
    QFile::remove(UnitTest::AnisotropyTest::TestOutput5);
#endif
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestFilterAvailability()
  {
    {
      // Now instantiate the AdaptiveAlignmentFeature Filter from the FilterManager
      QString filtName = "AdaptiveAlignmentFeature";
      FilterManager* fm = FilterManager::Instance();
      IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
      if(nullptr == filterFactory.get())
      {
        std::stringstream ss;
        ss << "Unable to initialize the AdaptiveAlignmentFeature while executing the AnisotropyTest.";
        DREAM3D_TEST_THROW_EXCEPTION(ss.str())
      }
    }

    {
      // Now instantiate the AdaptiveAlignmentMisorientation Filter from the FilterManager
      QString filtName = "AdaptiveAlignmentMisorientation";
      FilterManager* fm = FilterManager::Instance();
      IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
      if(nullptr == filterFactory.get())
      {
        std::stringstream ss;
        ss << "Unable to initialize the AdaptiveAlignmentMisorientation while executing the AnisotropyTest.";
        DREAM3D_TEST_THROW_EXCEPTION(ss.str())
      }
    }

    {
      // Now instantiate the AdaptiveAlignmentMutualInformation Filter from the FilterManager
      QString filtName = "AdaptiveAlignmentMutualInformation";
      FilterManager* fm = FilterManager::Instance();
      IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
      if(nullptr == filterFactory.get())
      {
        std::stringstream ss;
        ss << "Unable to initialize the AdaptiveAlignmentMutualInformation while executing the AnisotropyTest.";
        DREAM3D_TEST_THROW_EXCEPTION(ss.str())
      }
    }

    {
      // Now instantiate the SteinerCompact Filter from the FilterManager
      QString filtName = "SteinerCompact";
      FilterManager* fm = FilterManager::Instance();
      IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
      if(nullptr == filterFactory.get())
      {
        std::stringstream ss;
        ss << "Unable to initialize the SteinerCompact while executing the AnisotropyTest.";
        DREAM3D_TEST_THROW_EXCEPTION(ss.str())
      }
    }

    return 0;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addReadH5EBSDFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "ReadH5Ebsd";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    DataContainerArray::Pointer dca = DataContainerArray::New();

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      bool propWasSet;
      QVariant var;

      propWasSet = filter->setProperty("InputFile", UnitTest::AnisotropyTest::TestInput);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("ZStartIndex", UnitTest::AnisotropyTest::TestTifStartIndex);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("ZEndIndex", UnitTest::AnisotropyTest::TestTifEndIndex);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      QSet<QString> SelectedArrays = {"Phases", "EulerAngles", "Image Quality"};
      var.setValue(SelectedArrays);
      propWasSet = filter->setProperty("SelectedArrayNames", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("UseTransformations", true);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      filter->setDataContainerArray(dca);
      filter->preflight();
      int err = filter->getErrorCode();
      DREAM3D_REQUIRE_EQUAL(err, 0);

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addImportImageStackFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "ITKImportImageStack";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      bool propWasSet;
      QVariant var;

      StackFileListInfo input;
      input.InputPath = UnitTest::AnisotropyTest::InputDir + "/" + UnitTest::AnisotropyTest::TestTifExtension;
      input.StartIndex = UnitTest::AnisotropyTest::TestTifStartIndex;
      input.EndIndex = UnitTest::AnisotropyTest::TestTifEndIndex;
      input.FilePrefix = UnitTest::AnisotropyTest::TestTifPrefix;
      input.FileSuffix = UnitTest::AnisotropyTest::TestTifSuffix;
      input.IncrementIndex = 1;
      input.FileExtension = UnitTest::AnisotropyTest::TestTifExtension;
      input.PaddingDigits = UnitTest::AnisotropyTest::TestTifPaddingDigits;
      input.Ordering = 0;
      var.setValue(input);
      propWasSet = filter->setProperty("InputFileListInfo", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      FloatVec3Type origin = {0, 0, 0};
      var.setValue(origin);
      propWasSet = filter->setProperty("Origin", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      FloatVec3Type resolution = {1, 1, 1};
      var.setValue(resolution);
      propWasSet = filter->setProperty("Spacing", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      var.setValue(DataArrayPath("SEMImageDataContainer"));
      propWasSet = filter->setProperty("DataContainerName", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addMultiThresholdObjectsFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "MultiThresholdObjects";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      QVariant var;
      bool propWasSet;

      ComparisonInputs Thresholds;
      Thresholds.addInput("ImageDataContainer", "CellData", "Image Quality", 1, 800.0);

      var.setValue(Thresholds);
      propWasSet = filter->setProperty("SelectedThresholds", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("DestinationArrayName", SIMPL::GeneralData::Mask);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addConvertOrientationsFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "ConvertOrientations";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      bool propWasSet;
      QVariant var;

      propWasSet = filter->setProperty("InputType", 0);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("OutputType", 2);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      DataArrayPath path("ImageDataContainer", "CellData", "EulerAngles");
      var.setValue(path);
      propWasSet = filter->setProperty("InputOrientationArrayPath", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("OutputOrientationArrayName", "Quats");
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addEBSDSegmentFeatures(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "EBSDSegmentFeatures";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      bool propWasSet;

      propWasSet = filter->setProperty("MisorientationTolerance", 5.0);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("UseGoodVoxels", false);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addAdaptiveAlignmentFeatureFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "AdaptiveAlignmentFeature";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      bool propWasSet;

      propWasSet = filter->setProperty("GlobalCorrection", 2);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("ShiftX", 0);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("ShiftY", 0);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("WriteAlignmentShifts", true);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("AlignmentShiftFileName", UnitTest::AnisotropyTest::TestOutput1);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addAdaptiveAlignmentMisorientationFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "AdaptiveAlignmentMisorientation";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();

      QVariant var;
      bool propWasSet;

      propWasSet = filter->setProperty("GlobalCorrection", 1);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("MisorientationTolerance", 5.0f);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("UseGoodVoxels", false);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      DataArrayPath path("SEMImageDataContainer", "CellData", "ImageData");
      var.setValue(path);
      propWasSet = filter->setProperty("ImageDataArrayPath", var);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("WriteAlignmentShifts", true);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("AlignmentShiftFileName", UnitTest::AnisotropyTest::TestOutput2);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addAdaptiveAlignmentMutualInformationFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "AdaptiveAlignmentMutualInformation";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      bool propWasSet;

      propWasSet = filter->setProperty("GlobalCorrection", 0);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("UseGoodVoxels", false);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("MisorientationTolerance", 5.0f);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("WriteAlignmentShifts", true);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("AlignmentShiftFileName", UnitTest::AnisotropyTest::TestOutput3);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void addSteinerCompactFilter(const FilterPipeline::Pointer& pipeline)
  {
    QString filtName = "SteinerCompact";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

    if(nullptr != filterFactory.get())
    {
      AbstractFilter::Pointer filter = filterFactory->create();
      bool propWasSet;

      propWasSet = filter->setProperty("Plane", 0);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("Sites", 1);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("VtkOutput", true);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("VtkFileName", UnitTest::AnisotropyTest::TestOutput4);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("TxtOutput", true);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)
      propWasSet = filter->setProperty("TxtFileName", UnitTest::AnisotropyTest::TestOutput5);
      DREAM3D_REQUIRE_EQUAL(propWasSet, true)

      pipeline->pushBack(filter);
    }
    else
    {
      QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
      qDebug() << ss;
      DREAM3D_REQUIRE_EQUAL(0, 1)
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestAnisotropy()
  {
    Observer obs;
    FilterPipeline::Pointer pipeline = FilterPipeline::New();
    pipeline->addMessageReceiver(&obs);

    {
      addReadH5EBSDFilter(pipeline);
      int err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addImportImageStackFilter(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addConvertOrientationsFilter(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addAdaptiveAlignmentMisorientationFilter(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addMultiThresholdObjectsFilter(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addAdaptiveAlignmentFeatureFilter(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addAdaptiveAlignmentMutualInformationFilter(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addEBSDSegmentFeatures(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)

      addSteinerCompactFilter(pipeline);
      err = pipeline->preflightPipeline();
      DREAM3D_REQUIRE_EQUAL(err, 0)
      pipeline->execute();
      DREAM3D_REQUIRE_EQUAL(pipeline->getErrorCode(), 0)
    }

    if(!UnitTest::AnisotropyTest::PipelineJsonFile.isEmpty())
    {
      QFile outputFile(UnitTest::AnisotropyTest::PipelineJsonFile);
      QFileInfo info(outputFile);
      QString parentPath = info.absolutePath();
      QDir parentDir(parentPath);

      if(!parentDir.exists())
      {
        parentDir.mkpath(parentPath);
      }

      QJsonDocument doc(pipeline->toJson());

      if(outputFile.exists())
      {
        outputFile.remove();
      }
      if(outputFile.open(QIODevice::WriteOnly))
      {
        outputFile.write(doc.toJson());
        outputFile.close();
      }
    }
    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    std::cout << "###### AnisotropyFilterTest ######" << std::endl;
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestFilterAvailability());
    DREAM3D_REGISTER_TEST(TestAnisotropy())
    DREAM3D_REGISTER_TEST(RemoveTestFiles())
  }

public:
  AnisotropyFilterTest(const AnisotropyFilterTest&) = delete;            // Copy Constructor Not Implemented
  AnisotropyFilterTest(AnisotropyFilterTest&&) = delete;                 // Move Constructor Not Implemented
  AnisotropyFilterTest& operator=(const AnisotropyFilterTest&) = delete; // Copy Assignment Not Implemented
  AnisotropyFilterTest& operator=(AnisotropyFilterTest&&) = delete;      // Move Assignment Not Implemented
};
