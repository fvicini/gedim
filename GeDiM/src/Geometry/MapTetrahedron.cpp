#include "MapTetrahedron.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  MapTetrahedron::MapTetrahedronData MapTetrahedron::Compute(const Eigen::MatrixXd& vertices) const
  {
    MapTetrahedronData result;

    result.Q = Q(vertices);
    result.b = b(vertices);

    return result;
  }
  // ***************************************************************************
  MatrixXd MapTetrahedron::J(const MapTetrahedronData& mapData,
                             const MatrixXd& x) const
  {
    const unsigned int numPoints = x.cols();
    MatrixXd jacb(3, 3 * numPoints);

    for (unsigned int p = 0; p < numPoints; p++)
      jacb.block(0, 3 * p, 3, 3) = mapData.Q;

    return jacb;
  }
  // ***************************************************************************
  VectorXd MapTetrahedron::DetJ(const MapTetrahedronData& mapData,
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
