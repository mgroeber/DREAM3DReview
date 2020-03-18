/*
 * Your License or Copyright Information can go here
 */

#pragma once

#include <QtCore/QFile>
#include <QtCore/QString>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The DistanceTemplate class contains a templated function getDistance to find the distance, via a variety of
 * metrics, between two vectors of arbirtrary dimensons.  The developer should ensure that the pointers passed to
 * getDistance do indeed contain vectors of the same component dimensions and start at the desired tuples.
 */
class DistanceTemplate
{
public:
  using Self = DistanceTemplate;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer()
  {
    return Pointer(static_cast<Self*>(nullptr));
  }

  /**
   * @brief Returns the name of the class for DistanceTemplate
   */
  /**
   * @brief Returns the name of the class for DistanceTemplate
   */
  QString getNameOfClass() const
  {
    return QString("DistanceTemplate");
  }

  /**
   * @brief Returns the name of the class for DistanceTemplate
   */
  QString ClassName()
  {
    return QString("DistanceTemplate");
  }

  DistanceTemplate(){}
  virtual ~DistanceTemplate(){}

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  static QVector<QString> GetDistanceMetricsOptions()
  {
    QVector<QString> distMetricOptions;

    QString distMetric1 = "Euclidean";
    distMetricOptions.push_back(distMetric1);
    QString distMetric2 = "Squared Euclidean";
    distMetricOptions.push_back(distMetric2);
    QString distMetric3 = "Manhattan";
    distMetricOptions.push_back(distMetric3);
    QString distMetric4 = "Cosine";
    distMetricOptions.push_back(distMetric4);
    QString distMetric5 = "Pearson";
    distMetricOptions.push_back(distMetric5);
    QString distMetric6 = "Squared Pearson";
    distMetricOptions.push_back(distMetric6);

    return distMetricOptions;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  template<typename leftDataType, typename rightDataType, typename outDataType>
  static outDataType GetDistance(leftDataType* leftVector, rightDataType* rightVector, size_t compDims, int distMetric)
  {
    double dist = 0.0;
    double lVal = 0.0;
    double rVal = 0.0;

    double epsilon = std::numeric_limits<double>::min();

    // Euclidean
    if (distMetric == 0)
    {
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        dist += (lVal - rVal) * (lVal - rVal);
      }

      dist = sqrt(dist);
    }
    // Squared Euclidean
    else if (distMetric == 1)
    {
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        dist += (lVal - rVal) * (lVal - rVal);
      }
    }
    // Manhattan
    else if (distMetric == 2)
    {
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        dist += fabs(lVal - rVal);
      }
    }
    // Cosine
    else if (distMetric == 3)
    {
      double r = 0;
      double x = 0;
      double y = 0;
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        r += lVal * rVal;
        x += lVal * lVal;
        y += rVal * rVal;
      }
      dist = 1 - (r / (sqrt(x * y) + epsilon));
    }
    // Pearson
    else if (distMetric == 4)
    {
      double r = 0;
      double x = 0;
      double y = 0;
      double xAvg = 0;
      double yAvg = 0;
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        xAvg += lVal;
        yAvg += rVal;
      }
      xAvg /= static_cast<double>(compDims);
      yAvg /= static_cast<double>(compDims);
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        r += (lVal - xAvg) * (rVal - yAvg);
        x += (lVal - xAvg) * (lVal - xAvg);
        y += (rVal - yAvg) * (rVal - yAvg);
      }
      dist = 1 - (r / (sqrt(x * y) + epsilon));
    }
    // Squared Pearson
    else if (distMetric == 5)
    {
      double r = 0;
      double x = 0;
      double y = 0;
      double xAvg = 0;
      double yAvg = 0;
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        xAvg += lVal;
        yAvg += rVal;
      }
      xAvg /= static_cast<double>(compDims);
      yAvg /= static_cast<double>(compDims);
      for (size_t i=0;i<compDims;i++)
      {
        lVal = static_cast<double>(leftVector[i]);
        rVal = static_cast<double>(rightVector[i]);
        r += (lVal - xAvg) * (rVal - yAvg);
        x += (lVal - xAvg) * (lVal - xAvg);
        y += (rVal - yAvg) * (rVal - yAvg);
      }
      dist = 1 - ((r * r) / ((x * y) + epsilon));
    }

    // Return the correct primitive type for distance
    return static_cast<outDataType>(dist);
  }

private:
  DistanceTemplate(const DistanceTemplate&); // Copy Constructor Not Implemented
  void operator=(const DistanceTemplate&);   // Move assignment Not Implemented
};

