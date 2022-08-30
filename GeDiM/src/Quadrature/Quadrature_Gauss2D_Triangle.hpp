#ifndef __Quadrature_Gauss2D_Triangle_H
#define __Quadrature_Gauss2D_Triangle_H

#include "Eigen/Eigen"

namespace Gedim
{
  /// Gauss quadrature rule for triangles
  class Quadrature_Gauss2D_Triangle
  {
    public:
      static void FillPointsAndWeights(const unsigned int& order,
                                       Eigen::MatrixXd& points,
                                       Eigen::VectorXd& weights);
  };
}

#endif // __Quadrature_Gauss2D_Triangle_H
