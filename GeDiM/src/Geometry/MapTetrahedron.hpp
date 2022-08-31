#ifndef __MapTetrahedron_H
#define __MapTetrahedron_H

#include "Eigen/Eigen"

namespace Gedim
{
  class MapTetrahedron
  {
    public:
      struct MapTetrahedronData final
      {
          Eigen::Matrix3d Q;
          Eigen::Vector3d b;
      };

      /// Matrix Q for linear map x = Q * x_r + b from reference tetrahedron [0,1]x[0,1]x[0,1]/2 to tetrahedron with x points
      /// vertices the tetrahedron to map vertices, size 3 x 4
      /// return the resulting value, size 3 x 3
      inline Eigen::Matrix3d Q(const Eigen::MatrixXd& vertices) const
      {
        Eigen::Matrix3d Q;
        Q.row(0)<< vertices(0, 1) - vertices(0, 0), vertices(0, 2) - vertices(0, 0), vertices(0, 3) - vertices(0, 0);
        Q.row(1)<< vertices(1, 1) - vertices(1, 0), vertices(1, 2) - vertices(1, 0), vertices(1, 3) - vertices(1, 0);
        Q.row(2)<< vertices(2, 1) - vertices(2, 0), vertices(2, 2) - vertices(2, 0), vertices(2, 3) - vertices(2, 0);

        return Q;
      }

      /// Matrix Q for linear map x = Q * x_r + b from reference tetrahedron [0,1]x[0,1]x[0,1]/2 to tetrahedron with x points
      /// vertices the tetrahedron to map vertices, size 3 x 4
      /// return the resulting value, size 3 x 3
      inline Eigen::Vector3d b(const Eigen::MatrixXd& vertices) const
      {
        return vertices.col(0);
      }

    public:
      MapTetrahedron() {}
      ~MapTetrahedron() {}

      /// Map from the tetrahedron reference element [0,1]x[0,1]x[0,1]/2 to the polygon x = F(x_r) = Q * x_r + b
      /// \param vertices the tetrahedron to map vertices, size 3 x 4
      /// \return the map data
      MapTetrahedronData Compute(const Eigen::MatrixXd& vertices) const;

      /// Map from the tetrahedron reference element [0,1]x[0,1]x[0,1]/2 to the polygon x = F(x_r) = Q * x_r + b
      /// \param mapData the map data computed
      /// \param x points in reference tetrahedron, size 3 x numPoints
      /// \return the mapped points, size 3 x numPoints
      inline Eigen::MatrixXd F(const MapTetrahedronData& mapData,
                               const Eigen::MatrixXd& x) const
      { return (mapData.Q * x).colwise() + mapData.b; }
      /// Compute the jacobian matrix of the transformation F
      /// \param mapData the map data computed
      /// \param x points in reference tetrahedron, size 3 x numPoints
      /// \return the Q matrix for each points, size 2 x (2 * numPoints)
      Eigen::MatrixXd J(const MapTetrahedronData& mapData,
                        const Eigen::MatrixXd& x) const;
      /// Compute the determinant of the jacobian matrix of the trasformation
      /// \param mapData the map data computed
      /// \param x points in reference tetrahedron, size 3 x numPoints
      /// \return the determinant of Jacobian matrix for each points, size 1 x numPoints
      Eigen::VectorXd DetJ(const MapTetrahedronData& mapData,
                           const Eigen::MatrixXd& x) const;
  };
}

#endif // __MapTetrahedron_H
