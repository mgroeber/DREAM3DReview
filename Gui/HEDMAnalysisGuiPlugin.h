#pragma once

#include "HEDMAnalysis/HEDMAnalysisPlugin.h"

class HEDMAnalysisGuiPlugin : public HEDMAnalysisPlugin
{
  Q_OBJECT
  Q_INTERFACES(ISIMPLibPlugin)
  Q_PLUGIN_METADATA(IID "net.bluequartz.dream3d.HEDMAnalysisGuiPlugin")

public:
  HEDMAnalysisGuiPlugin();
  ~HEDMAnalysisGuiPlugin() override;

public:
  HEDMAnalysisGuiPlugin(const HEDMAnalysisGuiPlugin&) = delete;            // Copy Constructor Not Implemented
  HEDMAnalysisGuiPlugin(HEDMAnalysisGuiPlugin&&) = delete;                 // Move Constructor
  HEDMAnalysisGuiPlugin& operator=(const HEDMAnalysisGuiPlugin&) = delete; // Copy Assignment Not Implemented
  HEDMAnalysisGuiPlugin& operator=(HEDMAnalysisGuiPlugin&&) = delete;      // Move Assignment Not Implemented
};
