#ifndef _tdmsmetadata_h
#define _tdmsmetadata_h

#include <memory>

#include <unordered_map>

#include "TDMSDataType.hpp"
#include "TDMSProperty.h"

class TDMSObject;
class TDMSSegment;

class TDMSMetaData
{
public:
  virtual ~TDMSMetaData();
  TDMSMetaData(const TDMSMetaData&) = delete;
  TDMSMetaData& operator=(const TDMSMetaData&) = delete;

  struct MetaData
  {
    uint32_t RawDataIndex = 0;
    TDMSDataType::Pointer DataType;
    uint32_t ArrayDimension = 0;
    uint64_t NumberOfValues = 0;
    uint64_t TotalSegmentSize = 0;
    bool HasData = false;
  };

  typedef std::shared_ptr<TDMSMetaData> Pointer;

private:
  friend class TDMSObject;
  friend class TDMSSegment;

  TDMSMetaData();
  static Pointer New();

  void extendSegmentMetaData();

  void extendSegmentMetaData(size_t count);

  void extendSegmentMetaData(bool hasData);

  void resizeSegmentMetaData(size_t size);

  void readSegmentMetaData(std::ifstream& filestream, size_t index);

  std::vector<MetaData> m_SegmentMetaData;
  std::unordered_map<std::string, TDMSProperty::Pointer> m_Properties;
};

#endif
