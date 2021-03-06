#ifndef __TEST_GEOMETRY_H
#define __TEST_GEOMETRY_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestEquispaceCoordinates)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      ASSERT_EQ(geometryUtility.EquispaceCoordinates(1.0, true), vector<double>({ 0.0, 1.0 }));
      ASSERT_EQ(geometryUtility.EquispaceCoordinates(1.0, false), vector<double>({ }));
      ASSERT_EQ(geometryUtility.EquispaceCoordinates(0.5, true), vector<double>({ 0.0, 0.5, 1.0 }));
      ASSERT_EQ(geometryUtility.EquispaceCoordinates(0.5, false), vector<double>({ 0.5 }));
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestCompareValues)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check Compare1DValues
      {
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, -7.6547685574189633e-17), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare1DValues(-7.6547685574189633e-17, 0.0), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare1DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check Compare2DValues
      {
        ASSERT_EQ(geometryUtility.Compare2DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare2DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare2DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check Compare3DValues
      {
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare3DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check IsLenghtPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(0.0));
				ASSERT_FALSE(geometryUtility.IsValue1DPositive(geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue1DPositive(2 * geometryUtilityConfig.Tolerance));
        ASSERT_TRUE(geometryUtility.IsValue1DPositive(10.0));
      }

      // check IsAreaPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(0.0));
				ASSERT_TRUE(geometryUtility.IsValue2DPositive(geometryUtilityConfig.Tolerance));
				ASSERT_FALSE(geometryUtility.IsValue2DPositive(geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue2DPositive(10.0));
      }

      // check IsLenghtPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(0.0));
				ASSERT_TRUE(geometryUtility.IsValue3DPositive(geometryUtilityConfig.Tolerance));
				ASSERT_TRUE(geometryUtility.IsValue3DPositive(geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance));
				ASSERT_FALSE(geometryUtility.IsValue3DPositive(geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue3DPositive(10.0));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPlaneRotationMatrix)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check rotation matrix of plane 2D
      {
        const Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
        const Eigen::Vector3d planeOrigin(1.0, 2.0, 3.0);
        const Eigen::Matrix3d planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtility.PlaneTranslation(planeNormal,
                                                                                  planeOrigin);
        const Eigen::Vector3d planeOrigin2D = geometryUtility.RotatePointsFrom3DTo2D(planeOrigin,
                                                                                     planeRotationMatrix.transpose(),
                                                                                     planeTranslation);

        const Eigen::Vector3d point2D = geometryUtility.RotatePointsFrom3DTo2D(Eigen::Vector3d(0.0, 0.0, 3.0),
                                                                               planeRotationMatrix.transpose(),
                                                                               planeTranslation);


        ASSERT_DOUBLE_EQ(planeRotationMatrix(0, 0), 1.0);
        ASSERT_DOUBLE_EQ(planeRotationMatrix(1, 1), 1.0);
        ASSERT_DOUBLE_EQ(planeRotationMatrix(2, 2), 1.0);
        ASSERT_DOUBLE_EQ(planeTranslation[0], 1.0);
        ASSERT_DOUBLE_EQ(planeTranslation[1], 2.0);
        ASSERT_DOUBLE_EQ(planeTranslation[2], 3.0);
        ASSERT_DOUBLE_EQ(planeOrigin2D[0], 0.0);
        ASSERT_DOUBLE_EQ(planeOrigin2D[1], 0.0);
        ASSERT_DOUBLE_EQ(planeOrigin2D[2], 0.0);
        ASSERT_DOUBLE_EQ(point2D[0], -1.0);
        ASSERT_DOUBLE_EQ(point2D[1], -2.0);
        ASSERT_DOUBLE_EQ(point2D[2], 0.0);
      }

      // check rotation matrix of polygon 3D
      {
        Eigen::Vector3d planeNormal(1.0, 1.0, 1.0);
        planeNormal.normalize();
        const Eigen::Vector3d planeOrigin(1.0, 2.0, 3.0);
        const Eigen::Matrix3d planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtility.PlaneTranslation(planeNormal,
                                                                                  planeOrigin);
        const Eigen::Vector3d planeOrigin2D = geometryUtility.RotatePointsFrom3DTo2D(planeOrigin,
                                                                                     planeRotationMatrix.transpose(),
                                                                                     planeTranslation);
        const Eigen::Vector3d point2D = geometryUtility.RotatePointsFrom3DTo2D(Eigen::Vector3d(0.0, 0.0, 6.0),
                                                                               planeRotationMatrix.transpose(),
                                                                               planeTranslation);

        ASSERT_DOUBLE_EQ(planeRotationMatrix(0, 0), -7.0710678118654757e-01);
        ASSERT_DOUBLE_EQ(planeRotationMatrix(1, 1), -4.0824829046386307e-01);
        ASSERT_DOUBLE_EQ(planeRotationMatrix(2, 2), 5.7735026918962584e-01);
        ASSERT_DOUBLE_EQ(planeTranslation[0], 1.0);
        ASSERT_DOUBLE_EQ(planeTranslation[1], 2.0);
        ASSERT_DOUBLE_EQ(planeTranslation[2], 3.0);
        ASSERT_DOUBLE_EQ(planeOrigin2D[0], 0.0);
        ASSERT_DOUBLE_EQ(planeOrigin2D[1], 0.0);
        ASSERT_DOUBLE_EQ(planeOrigin2D[2], 0.0);
        ASSERT_DOUBLE_EQ(point2D[0], -7.071067811865476e-01);
        ASSERT_DOUBLE_EQ(point2D[1], 3.674234614174768e+00);
        ASSERT_DOUBLE_EQ(point2D[2], 0.0);
      }

      // check rotation matrix of other polygon 3D
      {
        Eigen::Vector3d normal(-1.0, -1.0, -1.0);
        normal.normalize();
        Eigen::Matrix3d rotationMatrix = geometryUtility.PlaneRotationMatrix(normal);

        ASSERT_DOUBLE_EQ(rotationMatrix(0, 0), 7.0710678118654757e-01);
        ASSERT_DOUBLE_EQ(rotationMatrix(1, 1), -4.0824829046386307e-01);
        ASSERT_DOUBLE_EQ(rotationMatrix(2, 2), -5.7735026918962584e-01);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestConvexHull)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check simple convex hull
      {
        Eigen::MatrixXd points(3, 4);
        points.col(0)<< 1.0, 1.0, 0.0;
        points.col(1)<< 0.0, 0.0, 0.0;
        points.col(2)<< 0.0, 1.0, 0.0;
        points.col(3)<< 1.0, 0.0, 0.0;

        vector<unsigned int> convexHull = geometryUtility.ConvexHull(points);
        ASSERT_EQ(convexHull, vector<unsigned int>({ 1, 3, 0, 2 }));

        Eigen::MatrixXd result(3, 4);
        result.col(0)<< 0.0, 0.0, 0.0;
        result.col(1)<< 1.0, 0.0, 0.0;
        result.col(2)<< 1.0, 1.0, 0.0;
        result.col(3)<< 0.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtility.ExtractPoints(points,
                                                convexHull), result);
      }

      // check complex convex hull
      {
        Eigen::MatrixXd points(3, 8);
        points.col(0)<< 20.0, 0.0, 0.0;
        points.col(1)<< 30.0, 60.0, 0.0;
        points.col(2)<< 50.0, 40.0, 0.0;
        points.col(3)<< 70.0, 30.0, 0.0;
        points.col(4)<< 55.0, 20.0, 0.0;
        points.col(5)<< 50.0, 10.0, 0.0;
        points.col(6)<< 0.0, 30.0, 0.0;
        points.col(7)<< 15.0, 25.0, 0.0;

        vector<unsigned int> convexHull = geometryUtility.ConvexHull(points);
        ASSERT_EQ(convexHull, vector<unsigned int>({ 6, 0, 5, 3, 1 }));

        Eigen::MatrixXd result(3, 5);
        result.col(0)<< 0.0, 30.0, 0.0;
        result.col(1)<< 20.0, 0.0, 0.0;
        result.col(2)<< 50.0, 10.0, 0.0;
        result.col(3)<< 70.0, 30.0, 0.0;
        result.col(4)<< 30.0, 60.0, 0.0;

        ASSERT_EQ(geometryUtility.ExtractPoints(points,
                                                convexHull), result);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestUnalignedPoints)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check two aligned points
      {
        Eigen::MatrixXd points(3, 2);
        points.col(0)<< 0.0, 0.0, 0.0;
        points.col(1)<< 1.0, 1.0, 0.0;

        vector<unsigned int> unalignedPoints = geometryUtility.UnalignedPoints(points);
        ASSERT_EQ(unalignedPoints, vector<unsigned int>({ 0, 1 }));

        Eigen::MatrixXd result(3, 2);
        result.col(0)<< 0.0, 0.0, 0.0;
        result.col(1)<< 1.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtility.ExtractPoints(points,
                                                unalignedPoints), result);
      }

      // check triangle points
      {
        Eigen::MatrixXd points(3, 3);
        points.col(0)<< 0.0, 0.0, 0.0;
        points.col(1)<< 1.0, 0.0, 0.0;
        points.col(2)<< 0.0, 1.0, 0.0;

        vector<unsigned int> unalignedPoints = geometryUtility.UnalignedPoints(points);
        ASSERT_EQ(unalignedPoints, vector<unsigned int>({ 0, 1, 2 }));

        Eigen::MatrixXd result(3, 3);
        result.col(0)<< 0.0, 0.0, 0.0;
        result.col(1)<< 1.0, 0.0, 0.0;
        result.col(2)<< 0.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtility.ExtractPoints(points,
                                                unalignedPoints), result);
      }

      // check triangle with aligned points
      {
        Eigen::MatrixXd points(3, 4);
        points.col(0)<< 0.0, 0.0, 0.0;
        points.col(1)<< 1.0, 0.0, 0.0;
        points.col(2)<< 0.5, 0.5, 0.0;
        points.col(3)<< 0.0, 1.0, 0.0;

        vector<unsigned int> unalignedPoints = geometryUtility.UnalignedPoints(points);
        ASSERT_EQ(unalignedPoints, vector<unsigned int>({ 0, 1, 3 }));

        Eigen::MatrixXd result(3, 3);
        result.col(0)<< 0.0, 0.0, 0.0;
        result.col(1)<< 1.0, 0.0, 0.0;
        result.col(2)<< 0.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtility.ExtractPoints(points,
                                                unalignedPoints), result);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_H
