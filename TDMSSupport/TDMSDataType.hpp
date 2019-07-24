#ifndef _tdmsdatatype_h
#define _tdmsdatatype_h

#include <cmath>
#include <fstream>

#include <QtCore/QDateTime>

#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/StringDataArray.h"

namespace TDMSDataTypeHelpers
{
struct TDMSTimeStamp
{
  uint64_t FractionsSecond = 0;
  int64_t SecondsSinceEpoch = 0;
};

inline QDateTime TDMSTimeStampToQDateTime(uint8_t* timestamp)
{
  TDMSTimeStamp time = *(reinterpret_cast<TDMSTimeStamp*>(timestamp));
  QDateTime dateTime(QDate(1904, 1, 1), QTime(0, 0, 0, 0), Qt::UTC);
  dateTime = dateTime.addSecs(time.SecondsSinceEpoch);
  double msecs = std::floor((static_cast<double>(time.FractionsSecond) * std::pow(2.0, -64.0)) * 1000.0);
  dateTime = dateTime.addMSecs(static_cast<qint64>(msecs));
  return dateTime;
}

inline QList<QDateTime> TDMSTimeStampsToQDateTimes(DataArray<uint8_t>::Pointer ptr)
{
  if(ptr->getSize() % 16 != 0)
  {
    return QList<QDateTime>();
  }
  QList<QDateTime> dateTimeList;
  for(size_t i = 0; i < ptr->getSize(); i += 16)
  {
    dateTimeList.push_back(TDMSTimeStampToQDateTime(ptr->getPointer(i)));
  }
  return dateTimeList;
}

template <typename T> inline typename DataArray<T>::Pointer ReadValueFromFile(std::ifstream& filestream, std::string name)
{
  typename DataArray<T>::Pointer data = DataArray<T>::CreateArray(1, QString::fromStdString(name), true);
  T* p = data->getPointer(0);
  filestream.read(reinterpret_cast<char*>(p), sizeof(T));
  return data;
}

inline DataArray<uint8_t>::Pointer ReadTimeStampFromFile(std::ifstream& filestream, std::string name)
{
  DataArray<uint8_t>::Pointer data = DataArray<uint8_t>::CreateArray(16, QString::fromStdString(name), true);
  uint8_t* p = data->getPointer(0);
  filestream.read(reinterpret_cast<char*>(p), 16);
  return data;
}

inline StringDataArray::Pointer ReadStringFromFile(std::ifstream& filestream, std::string name)
{
  StringDataArray::Pointer data = StringDataArray::CreateArray(1, QString::fromStdString(name), true);
  char buffer[4];
  filestream.read(buffer, 4);
  uint32_t length = *(reinterpret_cast<uint32_t*>(buffer));
  if(length != 0)
  {
    std::vector<char> string(length);
    filestream.read(string.data(), length);
    data->setValue(0, QString(string.data()));
  }
  return data;
}

template <typename T> inline void ReadArrayFromFile(std::ifstream& filestream, IDataArray::Pointer ptr, uint64_t pos, uint64_t bytes)
{
  typename DataArray<T>::Pointer data = std::dynamic_pointer_cast<DataArray<T>>(ptr);
  T* p = data->getTuplePointer(pos);
  filestream.read(reinterpret_cast<char*>(p), bytes);
  // uint32_t b = 0;
}

inline void ReadStringArrayFromFile(std::ifstream& filestream, IDataArray::Pointer ptr, uint64_t pos, uint64_t bytes)
{
  StringDataArray::Pointer data = std::dynamic_pointer_cast<StringDataArray>(ptr);
  std::vector<uint32_t> offsets(bytes + 1, 0);
  filestream.read(reinterpret_cast<char*>(offsets.data() + 1), bytes * sizeof(uint32_t));
  for(uint64_t s = 0; s < bytes; s++)
  {
    std::vector<char> string(offsets[s + 1] - offsets[s]);
    filestream.read(string.data(), string.size());
    data->setValue(s, QString(string.data()));
  }
}

template <typename T> inline typename DataArray<T>::Pointer GenerateArray(uint64_t numTuples, std::string name)
{
  typename DataArray<T>::Pointer data = DataArray<T>::CreateArray(numTuples, QString::fromStdString(name), false);
  return data;
}

template <typename T> inline void AllocateArray(IDataArray::Pointer ptr)
{
  typename DataArray<T>::Pointer data = std::dynamic_pointer_cast<DataArray<T>>(ptr);
  if(!data)
  {
    return;
  }
  data->allocate();
  data->initializeWithZeros();
}

inline StringDataArray::Pointer GenerateStringArray(uint64_t numTuples, std::string name)
{
  StringDataArray::Pointer data = StringDataArray::CreateArray(numTuples, QString::fromStdString(name), false);
  return data;
}

inline void AllocateStringArray(IDataArray::Pointer ptr)
{
  StringDataArray::Pointer data = std::dynamic_pointer_cast<StringDataArray>(ptr);
  if(!data)
  {
    return;
  }
  data->resize(data->getNumberOfTuples());
  data->initializeWithValue(QString("??????????"));
}

template <typename T> inline std::function<IDataArray::Pointer(std::ifstream&, std::string)> ValueReaderFactory()
{
  return ReadValueFromFile<T>;
}

inline std::function<IDataArray::Pointer(std::ifstream&, std::string)> TimeStampReaderFactory()
{
  return ReadTimeStampFromFile;
}

inline std::function<IDataArray::Pointer(std::ifstream&, std::string)> StringReaderFactory()
{
  return ReadStringFromFile;
}

template <typename T> inline std::function<IDataArray::Pointer(uint64_t, std::string)> ArrayGeneratorFactory()
{
  return GenerateArray<T>;
}

template <typename T> inline std::function<void(IDataArray::Pointer)> ArrayAllocatorFactory()
{
  return AllocateArray<T>;
}

template <typename T> inline std::function<void(std::ifstream&, IDataArray::Pointer, uint64_t, uint64_t)> ArrayReaderFactory()
{
  return ReadArrayFromFile<T>;
}

inline std::function<IDataArray::Pointer(uint64_t, std::string)> StringArrayGeneratorFactory()
{
  return GenerateStringArray;
}

inline std::function<void(IDataArray::Pointer)> StringArrayAllocatorFactory()
{
  return AllocateStringArray;
}

inline std::function<void(std::ifstream&, IDataArray::Pointer, uint64_t, uint64_t)> StringArrayReaderFactory()
{
  return ReadStringArrayFromFile;
}
}

class TDMSDataType
{
public:
  virtual ~TDMSDataType()
  {
  }
  TDMSDataType(const TDMSDataType&) = delete;
  TDMSDataType& operator=(const TDMSDataType&) = delete;

  typedef std::shared_ptr<TDMSDataType> Pointer;

  static Pointer New(const std::string& name, size_t size, std::function<IDataArray::Pointer(std::ifstream&, std::string)> valueReader,
                     std::function<IDataArray::Pointer(uint64_t, std::string)> arrayGenerator, std::function<void(IDataArray::Pointer)> arrayAllocator,
                     std::function<void(std::ifstream&, IDataArray::Pointer, uint64_t, uint64_t)> arrayReader)
  {
    Pointer shared(new TDMSDataType(name, size, valueReader, arrayGenerator, arrayAllocator, arrayReader));
    return shared;
  }

  bool operator==(const TDMSDataType& rhs)
  {
    return (m_Name == rhs.name());
  }

  bool operator!=(const TDMSDataType& rhs)
  {
    return !(*this == rhs);
  }

  std::string name() const
  {
    return m_Name;
  }

  size_t size() const
  {
    return m_Size;
  }

  IDataArray::Pointer readSingleValueFromFile(std::ifstream& filestream, std::string name)
  {
    return m_ValueReader(filestream, name);
  }

  IDataArray::Pointer generateDataArray(uint64_t numTuples, std::string name)
  {
    IDataArray::Pointer data = m_ArrayGenerator(numTuples, name);
    return data;
  }

  void allocateDataArray(IDataArray::Pointer ptr)
  {
    m_ArrayAllocator(ptr);
  }

  void readArrayFromFile(std::ifstream& filestream, IDataArray::Pointer ptr, uint64_t pos, uint64_t numValues)
  {
    uint64_t bytes = numValues;
    if(m_Size > 0)
    {
      bytes *= m_Size;
    }
    m_ArrayReader(filestream, ptr, pos, bytes);
  }

private:
  TDMSDataType(const std::string& name, size_t size, std::function<IDataArray::Pointer(std::ifstream&, std::string)> valueReader,
               std::function<IDataArray::Pointer(uint64_t, std::string)> arrayGenerator, std::function<void(IDataArray::Pointer)> arrayAllocator,
               std::function<void(std::ifstream&, IDataArray::Pointer, uint64_t, uint64_t)> arrayReader)
  : m_Name(name)
  , m_Size(size)
  , m_ValueReader(valueReader)
  , m_ArrayGenerator(arrayGenerator)
  , m_ArrayAllocator(arrayAllocator)
  , m_ArrayReader(arrayReader)
  {
  }

  std::string m_Name;
  size_t m_Size;
  std::function<IDataArray::Pointer(std::ifstream&, std::string)> m_ValueReader;
  std::function<IDataArray::Pointer(uint64_t, std::string)> m_ArrayGenerator;
  std::function<void(IDataArray::Pointer)> m_ArrayAllocator;
  std::function<void(std::ifstream&, IDataArray::Pointer, uint64_t, uint64_t)> m_ArrayReader;
};

#endif
