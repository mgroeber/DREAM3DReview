#ifndef _tdmsproperty_h
#define _tdmsproperty_h

#include <memory>

#include "TDMSDataType.hpp"

#include "SIMPLib/DataArrays/DataArray.hpp"

class IDataArray;
using IDataArrayShPtrType = std::shared_ptr<IDataArray>;

class TDMSMetaData;

class TDMSProperty
{
public:
  virtual ~TDMSProperty();
  TDMSProperty(const TDMSProperty&) = delete;
  TDMSProperty& operator=(const TDMSProperty&) = delete;

  typedef std::shared_ptr<TDMSProperty> Pointer;

  TDMSDataType::Pointer dataType()
  {
    return m_Type;
  }

  IDataArrayShPtrType value()
  {
    return m_Value;
  }

private:
  friend class TDMSMetaData;

  TDMSProperty(TDMSDataType::Pointer type, IDataArrayShPtrType value);
  static Pointer New(TDMSDataType::Pointer type, IDataArrayShPtrType value);

  TDMSDataType::Pointer m_Type;
  IDataArrayShPtrType m_Value;
};

#endif
