#ifndef __Quadrature_Gauss3D_Hexahedron_H
#define __Quadrature_Gauss3D_Hexahedron_H

#include "Eigen/Eigen"

namespace Gedim
{
  /// Gauss quadrature rule for Hexahedrons
  class Quadrature_Gauss3D_Hexahedron
  {
    public:
      static void FillPointsAndWeights(const unsigned int& order,
                                       Eigen::MatrixXd& points,
                                       Eigen::VectorXd& weights);
  };
}

#endif // __Quadrature_Gauss3D_Hexahedron_H
