#ifndef _tdmssegment_h
#define _tdmssegment_h

#include <memory>

#include "TDMSLeadIn.h"
#include "TDMSObject.h"

class TDMSFileProxy;

class TDMSSegment
{
public:
  virtual ~TDMSSegment();
  TDMSSegment(const TDMSSegment&) = delete;
  TDMSSegment& operator=(const TDMSSegment&) = delete;

  typedef std::shared_ptr<TDMSSegment> Pointer;

private:
  friend class TDMSFileProxy;

  TDMSSegment(std::ifstream& filestream, uint64_t currentSegment);
  static Pointer New(std::ifstream& filestream, uint64_t currentSegment);

  void initializeSegment();

  void readMetaData(std::unordered_map<std::string, TDMSObject::Pointer>& objects, std::vector<std::string>& order);

  void readRawData(std::unordered_map<std::string, TDMSObject::Pointer>& objects, std::vector<std::string>& order);

  void computeIncrementalChunks(std::unordered_map<std::string, TDMSObject::Pointer>& objects, std::vector<std::string>& order);

  std::ifstream& m_FileStream;
  uint64_t m_SegmentIndex;
  TDMSLeadIn::Pointer m_LeadIn;
  uint64_t m_RawDataPosition;
  uint64_t m_NextSegmentPosition;
  uint64_t m_NumberOfChunks;
  uint64_t m_TotalSegmentDataSize;
};

#endif
