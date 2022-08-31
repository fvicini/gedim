#include "MapTetrahedron.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  MatrixXd MapTetrahedron::F(const MatrixXd& vertices,
                             const MatrixXd& x) const
  {
    return (Q(vertices) * x).colwise() + b(vertices);
  }
  // ***************************************************************************
  MatrixXd MapTetrahedron::J(const MatrixXd& vertices,
                             const MatrixXd& x) const
  {
    const unsigned int numPoints = x.cols();
    MatrixXd jacb(3, 3 * numPoints);
    Matrix3d q = Q(vertices);

    for (unsigned int p = 0; p < numPoints; p++)
      jacb.block(0, 3 * p, 3, 3) = q;

    return jacb;
  }
  // ***************************************************************************
  VectorXd MapTetrahedron::DetJ(const MatrixXd& vertices,
                                const MatrixXd& x) const
  {
    MatrixXd jacb = J(vertices,
                      x);

    const unsigned int numPoints = x.cols();
    VectorXd detJacb(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
      detJacb[p] = jacb.block(0, 3 * p, 3, 3).determinant();

    return detJacb;
  }
  // ***************************************************************************
}
