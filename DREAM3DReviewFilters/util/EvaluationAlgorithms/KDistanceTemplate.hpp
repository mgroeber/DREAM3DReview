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
class KDistanceTemplate
{
public:
  using Self = KDistanceTemplate;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer()
  {
    return Pointer(static_cast<Self*>(nullptr));
  }

  /**
   * @brief Returns the name of the class for KDistanceTemplate
   */
  /**
   * @brief Returns the name of the class for KDistanceTemplate
   */
  QString getNameOfClass() const
  {
    return QString("KDistanceTemplate");
  }

  /**
   * @brief Returns the name of the class for KDistanceTemplate
   */
  QString ClassName()
  {
    return QString("KDistanceTemplate");
  }

  KDistanceTemplate() = default;
  virtual ~KDistanceTemplate() = default;

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
               int32_t minDist, int32_t distMetric)
  {
    typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inputIDataArray);
    T* inputData = inputDataPtr->getPointer(0);
    double* outputData = outputDataArray->getPointer(0);
    bool* mask = maskDataArray->getPointer(0);

    size_t numTuples = inputDataPtr->getNumberOfTuples();
    size_t cDims = inputDataPtr->getNumberOfComponents();

    double dist = 0.0;
    std::vector<double> neighbors(numTuples);

    int64_t progIncrement = static_cast<int64_t>(numTuples / 100);
    int64_t prog = 1;
    int64_t progressInt = 0;
    int64_t counter = 0;

    for(size_t i = 0; i < numTuples; i++)
    {
      if(filter->getCancel()) { return; }
      if(mask[i])
      {
        for (size_t j = 0; j < numTuples; j++)
        {
          if(mask[j])
          {
            dist = DistanceTemplate::GetDistance<T, T, double>(inputData + (cDims * j), inputData + (cDims * i), 
                                                               cDims, distMetric);
            neighbors[j] = dist;
          }
        }
        std::sort(neighbors.begin(), neighbors.end());
        if(minDist >= (numTuples - 1)) { outputData[i] = neighbors.back(); }
        else { outputData[i] = neighbors[minDist]; }
      }

      if(counter > prog)
      {
        progressInt = static_cast<int64_t>((static_cast<float>(counter) / numTuples) * 100.0f);
        QString ss = QObject::tr("Computing K Distances || Visited Point %1 of %2 || %3% Completed").arg(counter).arg(numTuples).arg(progressInt);
        filter->notifyStatusMessage(ss);
        prog = prog + progIncrement;
      }
      counter++;
    }
  }

private:
  KDistanceTemplate(const KDistanceTemplate&); // Copy Constructor Not Implemented
  void operator=(const KDistanceTemplate&);    // Move assignment Not Implemented
};

