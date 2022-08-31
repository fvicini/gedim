#include "MapHexahedron.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  MapHexahedron::MapHexahedronData MapHexahedron::Compute(const Eigen::MatrixXd& vertices,
                                                          const Eigen::MatrixXi& edges) const
  {
    MapHexahedron::MapHexahedronData result;

    result.Q = Q(vertices);
    result.b = b(vertices);

    return result;
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
