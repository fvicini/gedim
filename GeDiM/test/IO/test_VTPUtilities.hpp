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
    const unsigned int numGeometries = 4;
    const unsigned int numProperties = 2;

    Eigen::MatrixXd points(3, numGeometries);
    vector<vector<vector<double>>> properties(numProperties, vector<vector<double>>(numGeometries));
    vector<string> propertyLabels = { "Id", "Data" };

    vector<vector<double>>& id = properties[0];
    vector<vector<double>>& test = properties[1];

    // Initialize data
    for (unsigned g = 0; g < numGeometries; g++)
    {
      points.col(g)<< 1.0 + g, 0.0 + g, 0.0 + g;
      id[g].resize(1, g + 1);
      test[g].resize(1, 10.8 + g);
    }

    Gedim::VTKUtilities vtpUtilities;
    vtpUtilities.SetExportFormat(Gedim::VTKUtilities::Ascii);

    // Export to VTK
    //    for (unsigned int p = 0; p < numProperties; p++)
    //      vtpUtilities.AddProperty(propertyLabels[p],
    //                               Gedim::VTPProperty::Formats::Points);

    for (unsigned g = 0; g < numGeometries; g++)
    {
      vtpUtilities.AddPoint(points.col(g));

      //      for (unsigned int p = 0; p < numProperties; p++)
      //      {
      //        vtpUtilities.AddGeometryProperty(propertyLabels[p],
      //                                         properties[p][g].size(),
      //                                         properties[p][g].data());
      //      }
    }

    std::string exportFolder = "./Export/TestVTPUtilities/Test0D";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry0D.vtu");
  }
  // ***************************************************************************
  //  TEST(TestVTPUtilities, VTPUtilities_Test0D)
  //  {
  //    const unsigned int numGeometries = 4;
  //    const unsigned int numProperties = 2;

  //    vector<Eigen::Vector3d> points(numGeometries);
  //    vector<vector<vector<double>>> properties(numProperties, vector<vector<double>>(numGeometries));
  //    vector<string> propertyLabels = { "Id", "Data" };

  //    vector<vector<double>>& id = properties[0];
  //    vector<vector<double>>& test = properties[1];

  //    // Initialize data
  //    for (unsigned g = 0; g < numGeometries; g++)
  //    {
  //      points[g]<< 1.0 + g, 0.0 + g, 0.0 + g;
  //      id[g].resize(1, g + 1);
  //      test[g].resize(1, 10.8 + g);
  //    }

  //    Gedim::VTPUtilities vtpUtilities;
  //    vtpUtilities.SetExportFormat(Gedim::VTPUtilities::Ascii);

  //    // Export to VTK
  //    for (unsigned int p = 0; p < numProperties; p++)
  //      vtpUtilities.AddProperty(propertyLabels[p],
  //                               Gedim::VTPProperty::Formats::Points);

  //    for (unsigned g = 0; g < numGeometries; g++)
  //    {
  //      vtpUtilities.AddPoint(points[g]);

  //      for (unsigned int p = 0; p < numProperties; p++)
  //      {
  //        vtpUtilities.AddGeometryProperty(propertyLabels[p],
  //                                         properties[p][g].size(),
  //                                         properties[p][g].data());
  //      }
  //    }

  //    std::string exportFolder = "./Export/TestVTPUtilities/Test0D";
  //    Gedim::Output::CreateFolder(exportFolder);

  //    vtpUtilities.Export(exportFolder + "/Geometry0D.vtu");
  //  }
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
