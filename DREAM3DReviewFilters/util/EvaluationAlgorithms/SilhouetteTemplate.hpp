/* ============================================================================
* Copyright (c) 2009-2016 BlueQuartz Software, LLC
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
* Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
* contributors may be used to endorse or promote products derived from this software
* without specific prior written permission.
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
* The code contained herein was partially funded by the followig contracts:
*    United States Air Force Prime Contract FA8650-07-D-5800
*    United States Air Force Prime Contract FA8650-10-D-5210
*    United States Prime Contract Navy N00173-07-C-2068
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/DistanceTemplate.hpp"

template<typename T>
class SilhouetteTemplate
{
public:
  using Self = SilhouetteTemplate;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer()
  {
    return Pointer(static_cast<Self*>(nullptr));
  }

  /**
   * @brief Returns the name of the class for SilhouetteTemplate
   */
  /**
   * @brief Returns the name of the class for SilhouetteTemplate
   */
  QString getNameOfClass() const
  {
    return QString("SilhouetteTemplate");
  }

  /**
   * @brief Returns the name of the class for SilhouetteTemplate
   */
  QString ClassName()
  {
    return QString("SilhouetteTemplate");
  }

  SilhouetteTemplate(){}
  virtual ~SilhouetteTemplate(){}

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  bool operator()(IDataArray::Pointer p)
  {
    return (std::dynamic_pointer_cast<DataArray<T>>(p).get() != nullptr);
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void Execute(AbstractFilter* filter, IDataArray::Pointer inputIDataArray, DoubleArrayType::Pointer outputDataArray, BoolArrayType::Pointer maskDataArray, 
               size_t numClusters, Int32ArrayType::Pointer fIds, int distMetric)
  {
    typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inputIDataArray);
    T* inputData = inputDataPtr->getPointer(0);
    double* outputData = outputDataArray->getPointer(0);
    int32_t* fPtr = fIds->getPointer(0);
    bool* mask = maskDataArray->getPointer(0);

    size_t numTuples = inputDataPtr->getNumberOfTuples();
    size_t numCompDims = inputDataPtr->getNumberOfComponents();

    std::vector<double> inClusterDist(numTuples, 0.0);
    std::vector<double> outClusterMinDist(numTuples, 0.0);
    std::vector<double> numTuplesPerFeature(numClusters, 0.0);
    std::vector<std::vector<double>> clusterDist(numTuples, std::vector<double>(numClusters, 0.0));

    for(size_t i = 0; i < numTuples; i++)
    {
      if(mask[i]) { numTuplesPerFeature[fPtr[i]]++; }
    }

    int32_t cluster = 0;

    for(size_t i = 0; i < numTuples; i++)
    {
      if(mask[i])
      {
        for(size_t j = 0; j < numTuples; j++)
        {
          if(mask[j])
          {
            cluster = fPtr[j];
            clusterDist[i][cluster] += DistanceTemplate::GetDistance<T, T, double>(inputData + (numCompDims * i), 
                                                                                   inputData + (numCompDims * j), 
                                                                                   numCompDims,
                                                                                   distMetric);
          }
        }
      }
    }

    for(size_t i = 0; i < numTuples; i++)
    {
      if(mask[i])
      {
        for(size_t j = 0; j < numClusters; j++)
        {
          clusterDist[i][j] /= numTuplesPerFeature[j];
        }
      }
    }

    for(size_t i = 0; i < numTuples; i++)
    {
      if(mask[i])
      {
        cluster = fPtr[i];
        inClusterDist[i] = clusterDist[i][cluster];

        double dist = 0.0;
        double minDist = std::numeric_limits<double>::max();
        for(size_t j = 0; j < numClusters; j++)
        {
          if(cluster != j)
          {
            dist = clusterDist[i][j];
            if(dist < minDist)
            {
              minDist = dist;
              outClusterMinDist[i] = dist;
            }
          }
        }
      }
    }

    for(size_t i = 0; i < numTuples; i++)
    {
      if(mask[i])
      {
        outputData[i] = (outClusterMinDist[i] - inClusterDist[i]) / (std::max(outClusterMinDist[i], inClusterDist[i]));
      }
    }

  }

private:
  SilhouetteTemplate(const SilhouetteTemplate&); // Copy Constructor Not Implemented
  void operator=(const SilhouetteTemplate&);     // Move assignment Not Implemented
};

