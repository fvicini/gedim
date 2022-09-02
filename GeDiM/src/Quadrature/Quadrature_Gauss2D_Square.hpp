#ifndef __Quadrature_Gauss2D_Square_H
#define __Quadrature_Gauss2D_Square_H

#include "Eigen/Eigen"

namespace Gedim
{
  /// Gauss quadrature rule for Squares
  class Quadrature_Gauss2D_Square
  {
    public:
      static void FillPointsAndWeights(const unsigned int& order,
                                       Eigen::MatrixXd& points,
                                       Eigen::VectorXd& weights);
  };
}

#endif // __Quadrature_Gauss2D_Square_H
