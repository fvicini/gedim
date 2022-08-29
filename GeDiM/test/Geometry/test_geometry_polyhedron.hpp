#ifndef __TEST_GEOMETRY_POLYHEDRON_H
#define __TEST_GEOMETRY_POLYHEDRON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronBarycenter)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        ASSERT_TRUE(geometryUtility.PointsAreCoincident(geometryUtility.PolyhedronBarycenter(cube.Vertices),
                                                        Eigen::Vector3d(0.5, 0.5, 0.5)));
      }

      // check tetrahedron barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        ASSERT_TRUE(geometryUtility.PointsAreCoincident(geometryUtility.PolyhedronBarycenter(tetrahedron.Vertices),
                                                        Eigen::Vector3d(0.25, 0.25, 0.25)));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronEdgeTangents)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        Eigen::MatrixXd expectedEdgeTangents(3, 12);
        expectedEdgeTangents.col(0)<<   1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(1)<<   0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(2)<<  -1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(3)<<   0.0000000000000000e+00, -1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(4)<<   1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(5)<<   0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(6)<<  -1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(7)<<   0.0000000000000000e+00, -1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(8)<<   0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;
        expectedEdgeTangents.col(9)<<   0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;
        expectedEdgeTangents.col(10)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;
        expectedEdgeTangents.col(11)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;

        ASSERT_EQ(geometryUtility.PolyhedronEdgeTangents(cube.Vertices,
                                                         cube.Edges),
                  expectedEdgeTangents);
      }

      // check tetrahedron face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        Eigen::MatrixXd expectedEdgeTangents(3, 6);
        expectedEdgeTangents.col(0)<<  1.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00;
        expectedEdgeTangents.col(1)<<  0.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        expectedEdgeTangents.col(2)<< -1.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        expectedEdgeTangents.col(3)<<  0.0000000000000000e+00,  0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedEdgeTangents.col(4)<< -1.0000000000000000e+00,  0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedEdgeTangents.col(5)<<  0.0000000000000000e+00, -1.0000000000000000e+00, 1.0000000000000000e+00;

        ASSERT_EQ(geometryUtility.PolyhedronEdgeTangents(tetrahedron.Vertices,
                                                         tetrahedron.Edges),
                  expectedEdgeTangents);;

      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceVertices)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        vector<Eigen::MatrixXd> expectedFaceVertices(6, Eigen::MatrixXd(3, 4));
        expectedFaceVertices[0]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        expectedFaceVertices[1]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[2]<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[3]<< 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[4]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[5]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;

        ASSERT_EQ(geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                         cube.Faces),
                  expectedFaceVertices);
      }

      // check tetrahedron face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        vector<Eigen::MatrixXd> expectedFaceVertices(4, Eigen::MatrixXd(3, 3));

        expectedFaceVertices[0]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        expectedFaceVertices[1]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[2]<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[3]<< 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;

        ASSERT_EQ(geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                         tetrahedron.Faces),
                  expectedFaceVertices);

      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceNormals)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face normals
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::Vector3d barycenter(0.5, 0.5, 0.5);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);

        vector<Eigen::Vector3d> expectedFaceNormals(6);
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, -1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceNormals[2]<< -1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[4]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[5]<< +0.0000000000000000e+00, +1.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(geometryUtility.PolyhedronFaceNormals(faceVertices,
                                                        barycenter),
                  expectedFaceNormals);

      }

      // check tetrahedron face normals
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::Vector3d barycenter(0.25, 0.25, 0.25);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);

        vector<Eigen::Vector3d> expectedFaceNormals(4);
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, -1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[2]<< -1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +5.7735026918962584e-01, +5.7735026918962584e-01, +5.7735026918962584e-01;

        ASSERT_EQ(geometryUtility.PolyhedronFaceNormals(faceVertices,
                                                        barycenter),
                  expectedFaceNormals);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceEdgeDirections)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        vector<vector<bool>> expectedFaceEdgeDirections(6);
        expectedFaceEdgeDirections[0] = vector<bool>({ true, true, true, true });
        expectedFaceEdgeDirections[1] = vector<bool>({ true, true, true, true });
        expectedFaceEdgeDirections[2] = vector<bool>({ false, true, true, false });
        expectedFaceEdgeDirections[3] = vector<bool>({ true, true, false, false });
        expectedFaceEdgeDirections[4] = vector<bool>({ true, true, false, false });
        expectedFaceEdgeDirections[5] = vector<bool>({ false, true, true, false });

        ASSERT_EQ(expectedFaceEdgeDirections,
                  geometryUtility.PolyhedronFaceEdgeDirections(cube.Vertices,
                                                               cube.Edges,
                                                               cube.Faces));
      }

      // check tetrahedron face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        vector<vector<bool>> expectedFaceEdgeDirections(4);
        expectedFaceEdgeDirections[0] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[1] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[2] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[3] = vector<bool>({ true, true, false });

        ASSERT_EQ(expectedFaceEdgeDirections,
                  geometryUtility.PolyhedronFaceEdgeDirections(tetrahedron.Vertices,
                                                               tetrahedron.Edges,
                                                               tetrahedron.Faces));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_POLYHEDRON_H
