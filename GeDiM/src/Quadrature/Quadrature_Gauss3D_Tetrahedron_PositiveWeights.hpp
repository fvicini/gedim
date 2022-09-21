#ifndef __Quadrature_Gauss3D_Tetrahedron_PositiveWeights_H
#define __Quadrature_Gauss3D_Tetrahedron_PositiveWeights_H

#include "Eigen/Eigen"

namespace Gedim
{
  /// Gauss quadrature rule for Tetrahedrons
  class Quadrature_Gauss3D_Tetrahedron_PositiveWeights
  {
    public:
      static void FillPointsAndWeights(const unsigned int& order,
                                       Eigen::MatrixXd& points,
                                       Eigen::VectorXd& weights);
  };
}

#endif // __Quadrature_Gauss3D_Tetrahedron_PositiveWeights_H
