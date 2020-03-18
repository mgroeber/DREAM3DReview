/*
 * Your License or Copyright can go here
 */

#include "DREAM3DReviewPlugin.h"

#include <QtCore/QFile>
#include <QtCore/QFileInfo>

#include <QtCore/QTextStream>

#include "SIMPLib/Filtering/FilterFactory.hpp"

#include "SIMPLib/Filtering/FilterManager.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DREAM3DReviewPlugin::DREAM3DReviewPlugin()
: m_Version(DREAM3DReview::Version::Package())
, m_CompatibilityVersion(DREAM3DReview::Version::Package())
, m_Vendor("Open-Source")
, // Initialize DREAM3DReview's Vendor Name Here
    m_URL("http://www.github.com/dream3d/DREAM3DReview")
, // Initialize Company URL Here
    m_Location("")
, // Initialize DREAM3DReview library Location Here
    m_Description("Staging plugin for filters that are ready to be incorporated into official DREAM.3D plugin repositories")
, // Initialize DREAM3DReview's Description Here
    m_Copyright("Copyright is assigned by Filter")
, // Initialize DREAM3DReview's Copyright Here
    m_Filters(QList<QString>())
, // Initialize DREAM3DReview's List of Dependencies Here
    m_DidLoad(false)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DREAM3DReviewPlugin::~DREAM3DReviewPlugin() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getPluginFileName()
{
  return DREAM3DReviewConstants::DREAM3DReviewPluginFile;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getPluginDisplayName()
{
  return DREAM3DReviewConstants::DREAM3DReviewPluginDisplayName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getPluginBaseName()
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getVersion()
{
  return m_Version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getCompatibilityVersion()
{
  return m_CompatibilityVersion;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getVendor()
{
  return m_Vendor;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getURL()
{
  return m_URL;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getLocation()
{
  return m_Location;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getDescription()
{
  /* PLEASE UPDATE YOUR PLUGIN'S DESCRIPTION FILE.
  It is located at DREAM3DReview/Resources/DREAM3DReview/DREAM3DReviewDescription.txt */

  QFile licenseFile(":/DREAM3DReview/DREAM3DReviewDescription.txt");
  QFileInfo licenseFileInfo(licenseFile);
  QString text = "<<--Description was not read-->>";

  if(licenseFileInfo.exists())
  {
    if(licenseFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      QTextStream in(&licenseFile);
      text = in.readAll();
    }
  }
  return text;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getCopyright()
{
  return m_Copyright;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DREAM3DReviewPlugin::getLicense()
{
  /* PLEASE UPDATE YOUR PLUGIN'S LICENSE FILE.
  It is located at DREAM3DReview/Resources/DREAM3DReview/DREAM3DReviewLicense.txt */

  QFile licenseFile(":/DREAM3DReview/DREAM3DReviewLicense.txt");
  QFileInfo licenseFileInfo(licenseFile);
  QString text = "<<--License was not read-->>";

  if(licenseFileInfo.exists())
  {
    if(licenseFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      QTextStream in(&licenseFile);
      text = in.readAll();
    }
  }
  return text;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QMap<QString, QString> DREAM3DReviewPlugin::getThirdPartyLicenses()
{
  QMap<QString, QString> licenseMap;
  QList<QString> fileStrList;
  fileStrList.push_back(":/ThirdParty/HDF5.txt");

  fileStrList.push_back(":/ThirdParty/Qt.txt");
  fileStrList.push_back(":/ThirdParty/Qwt.txt");

  for(auto& iter : fileStrList)
  {
    QFile file(iter);
    QFileInfo licenseFileInfo(file);

    if(licenseFileInfo.exists())
    {
      if(file.open(QIODevice::ReadOnly | QIODevice::Text))
      {
        QTextStream in(&file);
        licenseMap.insert(licenseFileInfo.baseName(), in.readAll());
      }
    }
  }

  return licenseMap;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool DREAM3DReviewPlugin::getDidLoad()
{
  return m_DidLoad;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DREAM3DReviewPlugin::setDidLoad(bool didLoad)
{
  m_DidLoad = didLoad;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DREAM3DReviewPlugin::setLocation(QString filePath)
{
  m_Location = filePath;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DREAM3DReviewPlugin::writeSettings(QSettings& /*prefs*/)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DREAM3DReviewPlugin::readSettings(QSettings& /*prefs*/)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DREAM3DReviewPlugin::registerFilterWidgets(FilterWidgetManager* fwm)
{
}

#include "DREAM3DReviewFilters/RegisterKnownFilters.cpp"
