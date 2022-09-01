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
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceNormals[2]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[4]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[5]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(geometryUtility.PolyhedronFaceNormals(faceVertices),
                  expectedFaceNormals);
        ASSERT_EQ(geometryUtility.PolyhedronFaceNormalDirections(faceVertices,
                                                                 barycenter,
                                                                 expectedFaceNormals),
                  vector<bool>({ false, true, false, true, true, false }));

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
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[2]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +5.7735026918962584e-01, +5.7735026918962584e-01, +5.7735026918962584e-01;

        ASSERT_EQ(geometryUtility.PolyhedronFaceNormals(faceVertices),
                  expectedFaceNormals);
        ASSERT_EQ(geometryUtility.PolyhedronFaceNormalDirections(faceVertices,
                                                                 barycenter,
                                                                 expectedFaceNormals),
                  vector<bool>({ false, true, false, true }));
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

      // check cube face edge directions
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

      // check tetrahedron face edge directions
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

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceEdgeTangents)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::MatrixXd edgeTangents = geometryUtility.PolyhedronEdgeTangents(cube.Vertices,
                                                                                    cube.Edges);
        const vector<vector<bool>> faceEdgeDirections = geometryUtility.PolyhedronFaceEdgeDirections(cube.Vertices,
                                                                                                     cube.Edges,
                                                                                                     cube.Faces);

        vector<Eigen::MatrixXd> expectedFaceEdgeTangents(6);
        expectedFaceEdgeTangents[0] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1,  0,
                                       0,  1,  0, -1,
                                       0,  0,  0,  0).finished();
        expectedFaceEdgeTangents[1] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1,  0,
                                       0,  1,  0, -1,
                                       0,  0,  0,  0).finished();
        expectedFaceEdgeTangents[2] = (Eigen::MatrixXd(3, 4)<<
                                       -0,  0,  0, -0,
                                       1,  0, -1, -0,
                                       -0,  1,  0, -1).finished();
        expectedFaceEdgeTangents[3] = (Eigen::MatrixXd(3, 4)<<
                                       0,  0, -0, -0,
                                       1,  0, -1, -0,
                                       0,  1, -0, -1).finished();
        expectedFaceEdgeTangents[4] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1, -0,
                                       0,  0, -0, -0,
                                       0,  1, -0, -1).finished();
        expectedFaceEdgeTangents[5] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1, -0,
                                       -0,  0,  0, -0,
                                       -0,  1,  0, -1).finished();

        ASSERT_EQ(geometryUtility.PolyhedronFaceEdgeTangents(cube.Vertices,
                                                             cube.Edges,
                                                             cube.Faces,
                                                             faceEdgeDirections,
                                                             edgeTangents),
                  expectedFaceEdgeTangents);

      }

      // check tetrahedron face edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::MatrixXd edgeTangents = geometryUtility.PolyhedronEdgeTangents(tetrahedron.Vertices,
                                                                                    tetrahedron.Edges);
        const vector<vector<bool>> faceEdgeDirections = geometryUtility.PolyhedronFaceEdgeDirections(tetrahedron.Vertices,
                                                                                                     tetrahedron.Edges,
                                                                                                     tetrahedron.Faces);


        vector<Eigen::MatrixXd> expectedFaceEdgeTangents(4);
        expectedFaceEdgeTangents[0] = (Eigen::MatrixXd(3, 3)<<
                                       1, -1, -0,
                                       0,  1, -1,
                                       0,  0, -0).finished();
        expectedFaceEdgeTangents[1] = (Eigen::MatrixXd(3, 3)<<
                                       1, -1, -0,
                                       0,  0, -0,
                                       0,  1, -1).finished();
        expectedFaceEdgeTangents[2] = (Eigen::MatrixXd(3, 3)<<
                                       0,  0, -0,
                                       1, -1, -0,
                                       0,  1, -1).finished();
        expectedFaceEdgeTangents[3] = (Eigen::MatrixXd(3, 3)<<
                                       -1,  0,  1,
                                       1, -1, -0,
                                       0,  1, -1).finished();

        ASSERT_EQ(geometryUtility.PolyhedronFaceEdgeTangents(tetrahedron.Vertices,
                                                             tetrahedron.Edges,
                                                             tetrahedron.Faces,
                                                             faceEdgeDirections,
                                                             edgeTangents),
                  expectedFaceEdgeTangents);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceTranslations)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face translations
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);

        vector<Eigen::Vector3d> expectedFaceTranslations(6);
        expectedFaceTranslations[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceTranslations[2]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[4]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[5]<< +0.0000000000000000e+00, +1.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(expectedFaceTranslations,
                  geometryUtility.PolyhedronFaceTranslations(faceVertices));

      }

      // check tetrahedron face translations
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);

        vector<Eigen::Vector3d> expectedFaceTranslations(4);
        expectedFaceTranslations[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[2]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(expectedFaceTranslations,
                  geometryUtility.PolyhedronFaceTranslations(faceVertices));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceRotationMatrices)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face rotation matrices
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);
        const Eigen::Vector3d barycenter = geometryUtility.PolyhedronBarycenter(cube.Vertices);
        const vector<Eigen::Vector3d> faceNormals = geometryUtility.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtility.PolyhedronFaceTranslations(faceVertices);

        vector<Eigen::Matrix3d> expectedFaceRotationMatrices(6);
        expectedFaceRotationMatrices[0] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 9.9999999999999978e-01, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[1] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 9.9999999999999978e-01, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[2] = (Eigen::Matrix3d()<<
                                           0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00,
                                           9.9999999999999989e-01, -1.1102230246251565e-16,  0.0000000000000000e+00,
                                           1.1102230246251565e-16,  9.9999999999999989e-01,  0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[3] = (Eigen::Matrix3d()<<
                                           0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00,
                                           9.9999999999999989e-01, -1.1102230246251565e-16,  0.0000000000000000e+00,
                                           1.1102230246251565e-16,  9.9999999999999989e-01,  0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[4] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01,  0.0000000000000000e+00,  4.9650683064945600e-17,
                                           4.9650683064945576e-17,  2.4825341532472850e-17, -1.0000000000000000e+00,
                                           0.0000000000000000e+00,  9.9999999999999978e-01,  2.4825341532472825e-17).finished();
        expectedFaceRotationMatrices[5] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01,  0.0000000000000000e+00,  4.9650683064945600e-17,
                                           4.9650683064945576e-17,  2.4825341532472850e-17, -1.0000000000000000e+00,
                                           0.0000000000000000e+00,  9.9999999999999978e-01,  2.4825341532472825e-17).finished();

        ASSERT_EQ(expectedFaceRotationMatrices,
                  geometryUtility.PolyhedronFaceRotationMatrices(faceVertices,
                                                                 faceNormals,
                                                                 faceTranslations));

      }

      // check tetrahedron face rotation matrices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);

        const Eigen::Vector3d barycenter = geometryUtility.PolyhedronBarycenter(tetrahedron.Vertices);
        const vector<Eigen::Vector3d> faceNormals = geometryUtility.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtility.PolyhedronFaceTranslations(faceVertices);

        vector<Eigen::Matrix3d> expectedFaceRotationMatrices(4);
        expectedFaceRotationMatrices[0] = (Eigen::Matrix3d()<<
                                           1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[1] = (Eigen::Matrix3d()<<
                                           1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
                                           0.0000000000000000e+00,  0.0000000000000000e+00, -1.0000000000000000e+00,
                                           0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[2] = (Eigen::Matrix3d()<<
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                                           1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[3] = (Eigen::Matrix3d()<<
                                           -7.0710678118654724e-01, -4.0824829046386296e-01,  5.7735026918962595e-01,
                                           7.0710678118654724e-01, -4.0824829046386313e-01,  5.7735026918962573e-01,
                                           0.0000000000000000e+00,  8.1649658092772581e-01,  5.7735026918962651e-01).finished();

        ASSERT_EQ(expectedFaceRotationMatrices,
                  geometryUtility.PolyhedronFaceRotationMatrices(faceVertices,
                                                                 faceNormals,
                                                                 faceTranslations));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronCoordinateSystem)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube coordinate system
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        ASSERT_EQ(geometryUtility.PolyhedronCoordinateSystem(cube.Vertices,
                                                             cube.Edges),
                  vector<unsigned int>({ 0, 1, 3, 4 }));
      }

      // check tetrahedron face normals
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        ASSERT_EQ(geometryUtility.PolyhedronCoordinateSystem(tetrahedron.Vertices,
                                                             tetrahedron.Edges),
                  vector<unsigned int>({ 0, 1, 2, 3 }));
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
