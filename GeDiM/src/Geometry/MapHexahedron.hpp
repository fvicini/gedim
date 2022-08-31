#ifndef __MapHexahedron_H
#define __MapHexahedron_H

#include "Eigen/Eigen"

namespace Gedim
{
  class MapHexahedron
  {
    public:
      /// Matrix Q for linear map x = Q * x_r + b from reference Hexahedron [0,1]x[0,1]x[0,1] to Hexahedron with x points
      /// vertices the Hexahedron to map vertices, size 3 x 4
      /// return the resulting value, size 3 x 3
      inline Eigen::Matrix3d Q(const Eigen::MatrixXd& vertices) const
      {
        Eigen::Matrix3d Q;
        Q.row(0)<< vertices(0, 1) - vertices(0, 0), vertices(0, 2) - vertices(0, 0), vertices(0, 3) - vertices(0, 0);
        Q.row(1)<< vertices(1, 1) - vertices(1, 0), vertices(1, 2) - vertices(1, 0), vertices(1, 3) - vertices(1, 0);
        Q.row(2)<< vertices(2, 1) - vertices(2, 0), vertices(2, 2) - vertices(2, 0), vertices(2, 3) - vertices(2, 0);

        return Q;
      }

      /// Matrix Q for linear map x = Q * x_r + b from reference Hexahedron [0,1]x[0,1]x[0,1] to Hexahedron with x points
      /// vertices the Hexahedron to map vertices, size 3 x 4
      /// return the resulting value, size 3 x 3
      inline Eigen::Vector3d b(const Eigen::MatrixXd& vertices) const
      {
        return vertices.col(0);
      }

    public:
      MapHexahedron() {}
      ~MapHexahedron() {}

      /// Map from the Hexahedron reference element [0,1]x[0,1]x[0,1] to the polygon x = F(x_r) = Q * x_r + b
      /// vertices the Hexahedron to map vertices, size 3 x 4
      /// x points in reference Hexahedron, size 3 x numPoints
      /// \return the mapped points, size 3 x numPoints
      Eigen::MatrixXd F(const Eigen::MatrixXd& vertices,
                        const Eigen::MatrixXd& x) const;
      /// Compute the jacobian matrix of the transformation F
      /// x matrix of points between 0.0 and 1.0, size 2 x numPoints
      /// \return the Q matrix for each points, size 2 x (2 * numPoints)
      Eigen::MatrixXd J(const Eigen::MatrixXd& vertices,
                        const Eigen::MatrixXd& x) const;
      /// Compute the determinant of the jacobian matrix of the trasformation
      /// x matrix of points between 0.0 and 1.0, size 2 x numPoints
      /// \return the determinant of Jacobian matrix for each points, size 1 x numPoints
      Eigen::VectorXd DetJ(const Eigen::MatrixXd& vertices,
                           const Eigen::MatrixXd& x) const;
  };
}

#endif // __MapHexahedron_H
