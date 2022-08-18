#ifndef __TEST_VTPUtilities_H
#define __TEST_VTPUtilities_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "IOUtilities.hpp"
#include "VTPUtilities.hpp"

namespace GedimUnitTesting
{
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test0D)
  {
    unsigned int numPoints = 2;

    Gedim::VTPUtilities vtpUtilities;
    vtpUtilities.SetExportFormat(Gedim::VTPUtilities::Ascii);

    Eigen::MatrixXd points;
    vector<vector<double>> properties;

    points.setZero(3, 3 * numPoints);
    properties.resize(numPoints);

    vtpUtilities.AddProperty("TestSolution",
                             Gedim::VTPProperty::Formats::Points);

    for (unsigned p = 0; p < numPoints; p++)
    {
      points.col(p)<< 1.0 + p, 0.0 + p, 0.0 + p, 0;
      vtpUtilities.AddPoint(points.col(p));

      /// <li> Create Solution
      vector<double>& property = properties[p];
      property.resize(1, 10.8);

      vtpUtilities.AddGeometrySolution(property.size(),
                                       property.data());
    }

    std::string exportFolder = "./Export/Paraview";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry1D.vtu");
  }
  // ***************************************************************************
  //  TEST(TestVTPUtilities, VTPUtilities_Test1D)
  //  {
  //    unsigned int numberDomains = 2;

  //    VTPUtilities VTPUtilities;
  //    Output::AssertTest(VTPUtilities.Initialize(numberDomains), "%s: Initialize", __func__);
  //    Output::AssertTest(VTPUtilities.InitializeSolutions(1), "%s: InitializeSolutions", __func__);
  //    VTPUtilities.SetExportFormat(VTPUtilities::Ascii);

  //    vector<Point*> points;
  //    vector<Segment*> segments;
  //    vector<vector<double>*> solutions;

  //    segments.reserve(numberDomains);
  //    points.reserve(2 * numberDomains);
  //    solutions.reserve(numberDomains);

  //    for (unsigned d = 0; d < numberDomains; d++)
  //    {
  //      unsigned int id = d + 1;

  //      points.push_back(new Point(1.0 + d, 0.0 + d, 0.0 + d, 0));
  //      points.push_back(new Point(0.0 + d, 1.0 + d, 1.0 + d, 1));
  //      segments.push_back(new Segment(id));
  //      solutions.push_back(new vector<double>());

  //      Segment& segment = *segments[d];
  //      segment.AddVertex(*points[2*d]);
  //      segment.AddVertex(*points[2*d + 1]);

  //      Output::AssertTest(VTPUtilities.AddGeometry(id, segment), "%s: AddGeometry", __func__);

  //      /// <li> Create Solution
  //      vector<double>& solution = *solutions[d];
  //      solution.resize(segment.NumberOfVertices(), 10.8);

  //      Output::AssertTest(VTPUtilities.AddSolution(id, 0, solution.size(), solution.data()), "%s: AddSolution", __func__);
  //      Output::AssertTest(VTPUtilities.SetSolutionOptions(0, "TestSolution", VTPUtilities::Points), "%s: SetSolutionOptions", __func__);
  //    }

  //    string exportFolder = "./Export/Paraview";
  //    Output::CreateFolder(exportFolder);

  //    Output::AssertTest(VTPUtilities.Export(exportFolder + "/Geometry1D.vtu"), "%s: Export", __func__);
  //    for (unsigned d = 0; d < numberDomains; d++)
  //    {
  //      delete points[2*d];
  //      delete points[2*d + 1];
  //      delete segments[d];
  //      delete solutions[d];
  //    }
  //  }
  //  // ***************************************************************************
  //  TEST(TestVTPUtilities, VTPUtilities_Test2D)
  //  {
  //    double angle = 0.785398; // 45Â°
  //    Matrix3d rotationMatrix;
  //    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
  //    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
  //    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

  //    Vector3d translation(0.0, 2.0, 0.0);

  //    GeometryFactory geometryFactory;
  //    ShapeCreator shapeCreator(geometryFactory);

  //    vector<Vector3d> points {
  //      Vector3d {0.0, 0.0, 0.0},
  //      Vector3d {1.0, 0.0, 0.0},
  //      Vector3d {1.0, 1.0, 0.0},
  //      Vector3d {0.0, 1.0, 0.0},
  //    };

  //    Polygon& geometry = shapeCreator.CreatePolygon(points);
  //    geometry.Set2DPolygon();

  //    unsigned int id = geometry.Id();

  //    VTPUtilities VTPUtilities;
  //    Output::AssertTest(VTPUtilities.Initialize(1), "%s: Initialize", __func__);
  //    Output::AssertTest(VTPUtilities.InitializeSolutions(1), "%s: InitializeSolutions", __func__);
  //    VTPUtilities.SetExportFormat(VTPUtilities::Ascii);

  //    string exportFolder = "./Export/Paraview";
  //    Output::CreateFolder(exportFolder);

  //    Output::AssertTest(VTPUtilities.AddGeometry(id, geometry), "%s: AddGeometry", __func__);
  //    Output::AssertTest(VTPUtilities.AddRotationAndTranslation(id,
  //                                                              rotationMatrix,
  //                                                              translation), "%s: AddRotationAndTranslation", __func__);

  //    vector<double> solution(geometry.NumberOfVertices(), 10.8);
  //    Output::AssertTest(VTPUtilities.AddSolution(id, 0, solution.size(), solution.data()), "%s: AddSolution", __func__);
  //    Output::AssertTest(VTPUtilities.SetSolutionOptions(0, "TestSolution", VTPUtilities::Points), "%s: SetSolutionOptions", __func__);

  //    Output::AssertTest(VTPUtilities.Export(exportFolder + "/Geometry2D.vtu"), "%s: Export", __func__);
  //  }

  //  // ***************************************************************************
  //  TEST(TestVTPUtilities, VTPUtilities_Test3D)
  //  {
  //    double angle = 0.785398; // 45Â°
  //    Matrix3d rotationMatrix;
  //    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
  //    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
  //    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

  //    Vector3d translation(0.0, 2.0, 0.0);

  //    /// <li> Create Domain
  //    GeometryFactory geometryFactory;
  //    ShapeCreator shapeCreator(geometryFactory);

  //    Vector3d origin(0.0, 0.0, 0.0);
  //    Polyhedron& geometry = shapeCreator.CreateCube(origin,
  //                                                   1.0);

  //    unsigned int id = geometry.Id();

  //    VTPUtilities VTPUtilities;
  //    Output::AssertTest(VTPUtilities.Initialize(1), "%s: Initialize", __func__);
  //    Output::AssertTest(VTPUtilities.InitializeSolutions(1), "%s: InitializeSolutions", __func__);
  //    VTPUtilities.SetExportFormat(VTPUtilities::Ascii);

  //    string exportFolder = "./Export/Paraview";
  //    Output::CreateFolder(exportFolder);

  //    Output::AssertTest(VTPUtilities.AddGeometry(id, geometry), "%s: AddGeometry", __func__);
  //    Output::AssertTest(VTPUtilities.AddRotationAndTranslation(id,
  //                                                              rotationMatrix,
  //                                                              translation), "%s: AddRotationAndTranslation", __func__);

  //    vector<double> solution(geometry.NumberOfVertices(), 10.8);
  //    Output::AssertTest(VTPUtilities.AddSolution(id, 0, solution.size(), solution.data()), "%s: AddSolution", __func__);
  //    Output::AssertTest(VTPUtilities.SetSolutionOptions(0, "TestSolution", VTPUtilities::Points), "%s: SetSolutionOptions", __func__);

  //    Output::AssertTest(VTPUtilities.Export(exportFolder + "/Geometry3D.vtu"), "%s: Export", __func__);
  //  }
  // ***************************************************************************
}

#endif // __TEST_VTPUtilities_H
