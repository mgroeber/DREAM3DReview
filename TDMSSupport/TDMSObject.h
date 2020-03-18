#ifndef _tdmsobject_h

class IDataArray;
using IDataArrayShPtrType = std::shared_ptr<IDataArray>;

#define _tdmsobject_h

#include <memory>

#include "TDMSDataType.hpp"
#include "TDMSMetaData.h"

class TDMSFileProxy;
class TDMSSegment;

class TDMSObject
{
public:
  virtual ~TDMSObject();
  TDMSObject(const TDMSObject&) = delete;
  TDMSObject& operator=(const TDMSObject&) = delete;

  typedef std::shared_ptr<TDMSObject> Pointer;

  using EnumType = uint8_t;

  enum class Type : EnumType
  {
    Root,
    Group,
    Channel
  };

  std::string path()
  {
    return m_Path;
  }

  std::string baseName()
  {
    return m_BaseName;
  }

  TDMSDataType::Pointer dataType()
  {
    return m_DataType;
  }

  IDataArrayShPtrType data()
  {
    return m_Data;
  }

  Type objectType()
  {
    return m_ObjectType;
  }

  std::unordered_map<std::string, TDMSProperty::Pointer> properties()
  {
    return m_MetaData->m_Properties;
  }

private:
  friend class TDMSFileProxy;
  friend class TDMSSegment;

  TDMSObject(const std::string& path);
  static Pointer New(const std::string& path);

  void populateMetaData();

  void updateSizeInformation(uint64_t numChunks, size_t index);

  void generateDataArray();

  void allocate();

  void readRawData(std::ifstream& filestream, uint64_t index);

  std::string parseChannelName();

  std::string parseGroupName();

  void determineObjectType();

  std::string m_Path;
  std::string m_BaseName;
  TDMSMetaData::Pointer m_MetaData;
  TDMSDataType::Pointer m_DataType;
  uint32_t m_ArrayDimension;
  uint64_t m_NumberOfValues;
  uint64_t m_TotalSize;
  size_t m_CurrentDataPosition;
  bool m_HasData;
  bool m_HasInitializedMetaData;
  IDataArrayShPtrType m_Data;
  Type m_ObjectType;
};

#endif
