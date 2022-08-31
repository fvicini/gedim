#include "MapHexahedron.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Gedim
{
  Eigen::MatrixXd MapHexahedron::ComputeFourVertices(const Eigen::MatrixXd& vertices,
                                                     const Eigen::MatrixXi& edges) const
  {
    Eigen::MatrixXd fourVertices(3, 4);

    fourVertices.col(0)<< vertices.col(0);
    unsigned int vertexIndex = 1;
    for(unsigned int e = 0; e < edges.cols(); e++)
    {
      const unsigned int& originId = edges(0, e);
      const unsigned int& endId = edges(1, e);

      bool first = (originId == 0);
      bool second = (endId == 0);
      if (first || second)
      {
        if (first)
          fourVertices.col(vertexIndex++)<< vertices.col(endId);
        else
          fourVertices.col(vertexIndex++)<< vertices.col(originId);
      }
    }

    if (vertexIndex != 4)
      throw runtime_error("Hexaedron connectivity is not correct");

    return fourVertices;
  }
  // ***************************************************************************
  bool MapHexahedron::TestMapConfiguration(const Eigen::MatrixXd& vertices,
                                           const Eigen::MatrixXd& referencePoints,
                                           const unsigned int& secondVertexIndex,
                                           const unsigned int& thirdVertexIndex,
                                           const unsigned int& fourthVertexIndex,
                                           MapHexahedron::MapHexahedronData& result) const
  {
    result.Q = Q(vertices.col(0),
                 vertices.col(secondVertexIndex),
                 vertices.col(thirdVertexIndex),
                 vertices.col(fourthVertexIndex));
    result.b = b(vertices.col(0));

    return (geometryUtility.IsValue1DZero((vertices -
                                           F(result,
                                             referencePoints)).norm()));
  }
  // ***************************************************************************
  MapHexahedron::MapHexahedronData MapHexahedron::Compute(const Eigen::MatrixXd& vertices,
                                                          const Eigen::MatrixXi& edges) const
  {
    MapHexahedronData result;

    MatrixXd referencePoints;
    referencePoints.resize(3, 8);
    referencePoints.col(0)<< 0.0, 0.0, 0.0;
    referencePoints.col(1)<< 1.0, 0.0, 0.0;
    referencePoints.col(2)<< 1.0, 1.0, 0.0;
    referencePoints.col(3)<< 0.0, 1.0, 0.0;
    referencePoints.col(4)<< 0.0, 0.0, 1.0;
    referencePoints.col(5)<< 1.0, 0.0, 1.0;
    referencePoints.col(6)<< 1.0, 1.0, 1.0;
    referencePoints.col(7)<< 0.0, 1.0, 1.0;

    Eigen::MatrixXd fourVertices = ComputeFourVertices(vertices,
                                                       edges);

    if (TestMapConfiguration(fourVertices,
                             referencePoints,
                             1, 2, 3,
                             result))
      return result;

    if (TestMapConfiguration(fourVertices,
                             referencePoints,
                             1, 3, 2,
                             result))
      return result;

    if (TestMapConfiguration(fourVertices,
                             referencePoints,
                             2, 1, 3,
                             result))
      return result;

    if (TestMapConfiguration(fourVertices,
                             referencePoints,
                             2, 3, 1,
                             result))
      return result;

    if (TestMapConfiguration(fourVertices,
                             referencePoints,
                             3, 1, 2,
                             result))
      return result;

    if (TestMapConfiguration(fourVertices,
                             referencePoints,
                             3, 2, 1,
                             result))
      return result;

    throw runtime_error("Hexahedron cannot be mapped");
  }
  // ***************************************************************************
  MatrixXd MapHexahedron::F(const MapHexahedronData& mapData,
                            const MatrixXd& x) const
  {
    return (mapData.Q * x).colwise() + mapData.b;
  }
  // ***************************************************************************
  MatrixXd MapHexahedron::J(const MapHexahedronData& mapData,
                            const MatrixXd& x) const
  {
    const unsigned int numPoints = x.cols();
    MatrixXd jacb(3, 3 * numPoints);
    Matrix3d q = mapData.Q;

    for (unsigned int p = 0; p < numPoints; p++)
      jacb.block(0, 3 * p, 3, 3) = q;

    return jacb;
  }
  // ***************************************************************************
  VectorXd MapHexahedron::DetJ(const MapHexahedronData& mapData,
                               const MatrixXd& x) const
  {
    MatrixXd jacb = J(mapData,
                      x);

    const unsigned int numPoints = x.cols();
    VectorXd detJacb(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
      detJacb[p] = jacb.block(0, 3 * p, 3, 3).determinant();

    return detJacb;
  }
  // ***************************************************************************
}
