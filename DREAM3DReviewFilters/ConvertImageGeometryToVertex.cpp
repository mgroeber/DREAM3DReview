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
#include "ConvertImageGeometryToVertex.h"

#include <QtCore/QTextStream>

#include "SIMPLib/SIMPLibVersion.h"
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

enum createdPathID : RenameDataPath::DataID_t
{
  DataContainerID = 1
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ConvertImageGeometryToVertex::ConvertImageGeometryToVertex() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ConvertImageGeometryToVertex::~ConvertImageGeometryToVertex() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ConvertImageGeometryToVertex::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  {
    AttributeMatrixSelectionFilterParameter::RequirementType req = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Cell Attribute Matrix", CellAttributeMatrixPath, FilterParameter::Category::RequiredArray, ConvertImageGeometryToVertex, req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ConvertImageGeometryToVertex::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ConvertImageGeometryToVertex::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

  DataContainer::Pointer dc = getDataContainerArray()->getPrereqDataContainer(this, getCellAttributeMatrixPath());
  AttributeMatrix::Pointer srcCellAttrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getCellAttributeMatrixPath(), -301);
  if(getErrorCode() < 0)
  {
    return;
  }

  QString newPath = getCellAttributeMatrixPath().getDataContainerName() + "_vertex";
  DataContainer::Pointer dc2 = getDataContainerArray()->createNonPrereqDataContainer(this, newPath, DataContainerID);

  size_t numElements = dc->getGeometry()->getNumberOfElements();
  DataArray<float>::Pointer verts = VertexGeom::CreateSharedVertexList(numElements, !getInPreflight());
  VertexGeom::Pointer vertex = VertexGeom::CreateGeometry(verts, SIMPL::Geometry::VertexGeometry);
  dc2->setGeometry(vertex);

  AttributeMatrix::Pointer destCellAttrMat = srcCellAttrMat->deepCopy(getInPreflight());
  destCellAttrMat->setType(AttributeMatrix::Type::Vertex);
  dc2->addOrReplaceAttributeMatrix(destCellAttrMat);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ConvertImageGeometryToVertex::execute()
{
  dataCheck();

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());
  ImageGeom::Pointer imageGeom = m->getGeometryAs<ImageGeom>();
  size_t numElements = imageGeom->getNumberOfElements();

  DataContainer::Pointer m2 = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName() + "_vertex");
  VertexGeom::Pointer vertexGeom = m2->getGeometryAs<VertexGeom>();
  DataArray<float>::Pointer vertsPtr = vertexGeom->getVertices();
  float* verts = vertsPtr->getPointer(0);

  double coords[3];
  for (size_t i = 0; i < numElements; i++)
  {
    imageGeom->getCoords(i, coords);
    verts[3 * i] = static_cast<float>(coords[0]);
    verts[3 * i + 1] = static_cast<float>(coords[1]);
    verts[3 * i + 2] = static_cast<float>(coords[2]);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ConvertImageGeometryToVertex::newFilterInstance(bool copyFilterParameters) const
{
  ConvertImageGeometryToVertex::Pointer filter = ConvertImageGeometryToVertex::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::getGroupName() const
{
  return SIMPL::FilterGroups::CoreFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ConvertImageGeometryToVertex::getUuid() const
{
  return QUuid("{9ac220b9-14f9-581a-9bac-5714466297cb}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::getHumanLabel() const
{
  return "Convert Image Geometry To Vertex Geometry";
}

// -----------------------------------------------------------------------------
ConvertImageGeometryToVertex::Pointer ConvertImageGeometryToVertex::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ConvertImageGeometryToVertex> ConvertImageGeometryToVertex::New()
{
  struct make_shared_enabler : public ConvertImageGeometryToVertex
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::getNameOfClass() const
{
  return QString("ConvertImageGeometryToVertex");
}

// -----------------------------------------------------------------------------
QString ConvertImageGeometryToVertex::ClassName()
{
  return QString("ConvertImageGeometryToVertex");
}

// -----------------------------------------------------------------------------
void ConvertImageGeometryToVertex::setCellAttributeMatrixPath(const DataArrayPath& value)
{
  m_CellAttributeMatrixPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ConvertImageGeometryToVertex::getCellAttributeMatrixPath() const
{
  return m_CellAttributeMatrixPath;
}
