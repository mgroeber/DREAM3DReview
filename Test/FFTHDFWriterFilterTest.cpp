// -----------------------------------------------------------------------------
// Insert your license & copyright information here
// -----------------------------------------------------------------------------

#include <QtCore/QFile>

#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"
#include "SIMPLib/SIMPLib.h"
#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

class FFTHDFWriterFilterTest
{

  public:
    FFTHDFWriterFilterTest() = default;
    virtual ~FFTHDFWriterFilterTest() = default;

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    void RemoveTestFiles()
    {
#if REMOVE_TEST_FILES
      QFile::remove(UnitTest::FFTHDFWriterFilterTest::TestFile1);
      QFile::remove(UnitTest::FFTHDFWriterFilterTest::TestFile2);
#endif
    }

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    int TestFilterAvailability()
    {
      // Now instantiate the FFTHDFWriterFilterTest Filter from the FilterManager
      QString filtName = "FFTHDFWriterFilter";
      FilterManager* fm = FilterManager::Instance();
      IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
      if(nullptr == filterFactory.get())
      {
        std::stringstream ss;
        ss << "The FFTHDFWriterFilterTest Requires the use of the " << filtName.toStdString() << " filter which is found in the MASSIFUtilities Plugin";
        DREAM3D_TEST_THROW_EXCEPTION(ss.str())
      }
      return 0;
    }

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    int TestFFTHDFWriterFilterTest()
    {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      /* Please write FFTHDFWriterFilterTest test code here.
       *
       * Your IO test files are:
       * UnitTest::FFTHDFWriterFilterTest::TestFile1
       * UnitTest::FFTHDFWriterFilterTest::TestFile2
       *
       * SIMPLib provides some macros that will throw exceptions when a test fails
       * and thus report that during testing. These macros are located in the
       * SIMPLib/Common/UnitTestSupport.hpp file. Some examples are:
       *
       * SIMPLib_REQUIRE_EQUAL(foo, 0)
       * This means that if the variable foo is NOT equal to Zero then test will fail
       * and the current test will exit immediately. If there are more tests registered
       * with the SIMPLib_REGISTER_TEST() macro, the next test will execute. There are
       * lots of examples in the SIMPLib/Test folder to look at.
       */
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      int foo = 0;
      DREAM3D_REQUIRE_EQUAL(foo, 0)

      return EXIT_SUCCESS;
    }

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    void operator()()
    {
      std::cout << "###### FFTHDFWriterFilterTest ######" << std::endl;

      int err = EXIT_SUCCESS;

      DREAM3D_REGISTER_TEST(TestFilterAvailability());

      DREAM3D_REGISTER_TEST(TestFFTHDFWriterFilterTest())

      DREAM3D_REGISTER_TEST(RemoveTestFiles())
    }

  public:
    FFTHDFWriterFilterTest(const FFTHDFWriterFilterTest&) = delete;            // Copy Constructor Not Implemented
    FFTHDFWriterFilterTest(FFTHDFWriterFilterTest&&) = delete;                 // Move Constructor Not Implemented
    FFTHDFWriterFilterTest& operator=(const FFTHDFWriterFilterTest&) = delete; // Copy Assignment Not Implemented
    FFTHDFWriterFilterTest& operator=(FFTHDFWriterFilterTest&&) = delete;      // Move Assignment Not Implemented
};
