#ifndef _tdmsfileproxy_h
#define _tdmsfileproxy_h

#include <memory>

#include "TDMSObject.h"
#include "TDMSSegment.h"

class TDMSFileProxy
{
public:
  virtual ~TDMSFileProxy();
  TDMSFileProxy(const TDMSFileProxy&) = delete;
  TDMSFileProxy& operator=(const TDMSFileProxy&) = delete;

  typedef std::shared_ptr<TDMSFileProxy> Pointer;
  static Pointer New(const std::string& file);

  void readMetaData();

  void readRawData();

  void allocateObjects();

  std::unordered_map<std::string, TDMSObject::Pointer> objects()
  {
    return m_Objects;
  }

  TDMSObject::Pointer rootObject();

  std::unordered_map<std::string, TDMSObject::Pointer> groupObjects();

  std::unordered_map<std::string, TDMSObject::Pointer> channelObjects();

private:
  TDMSFileProxy(const std::string& file);

  std::unordered_map<std::string, TDMSObject::Pointer> extractObjectsOfType(TDMSObject::Type type);

  std::string m_File;
  std::ifstream m_FileStream;
  std::list<TDMSSegment::Pointer> m_Segments;
  std::unordered_map<std::string, TDMSObject::Pointer> m_Objects;
  std::vector<std::string> m_ObjectOrder;
  bool m_ObjectsAllocated;
  bool m_MetaDataRead;
};

#endif
