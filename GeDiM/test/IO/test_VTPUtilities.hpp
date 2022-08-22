#ifndef __TEST_VTPUtilities_H
#define __TEST_VTPUtilities_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "IOUtilities.hpp"
#include "VTPUtilities.hpp"

namespace GedimUnitTesting
{
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test0D)
  {
    const unsigned int numGeometries = 4;

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::Vector3d geometry(1.0 + g,
                               0.0 + g,
                               0.0 + g);

      vector<double> id(1, g + 1);
      vector<double> data(1, 10.8 + g);

      vtpUtilities.AddPoint(geometry,
                            {
                              {
                                "Id",
                                Gedim::VTPProperty::Formats::Points,
                                static_cast<unsigned int>(id.size()),
                                id.data()
                              },
                              {
                                "Data",
                                Gedim::VTPProperty::Formats::Points,
                                static_cast<unsigned int>(data.size()),
                                data.data()
                              }
                            });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry0D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test1D)
  {
    const unsigned int numGeometries = 4;

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::MatrixXd geometry(3, 2);

      geometry.col(0)<< 1.0 + g, 0.0 + g, 0.0 + g;
      geometry.col(1)<< 0.0 + g, 1.0 + g, 1.0 + g;

      vector<double> id(geometry.cols());
      vector<double> data(geometry.cols());

      for (unsigned int v = 0; v < geometry.cols(); v++)
      {
        id[v] = g + v;
        data[v] = 10.8 + g + v;
      }

      vtpUtilities.AddSegment(geometry,
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Points,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Data",
                                  Gedim::VTPProperty::Formats::Points,
                                  static_cast<unsigned int>(data.size()),
                                  data.data()
                                }
                              });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry1D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test2D)
  {
    const unsigned int numGeometries = 4;


    double angle = M_PI / 4.0; // 45 degrees
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

    Eigen::Vector3d translation(0.0, 2.0, 0.0);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::MatrixXd geometry(3, 4);
      geometry.col(0)<< 0.0 + g, 0.0 + g, 0.0;
      geometry.col(1)<< 1.0 + g, 0.0 + g, 0.0;
      geometry.col(2)<< 1.0 + g, 1.0 + g, 0.0;
      geometry.col(3)<< 0.0 + g, 1.0 + g, 0.0;

      geometry = geometryUtilities.RotatePointsFrom2DTo3D(geometry,
                                                          rotationMatrix,
                                                          translation);

      vector<double> id(geometry.cols());
      vector<double> data(geometry.cols());

      for (unsigned int v = 0; v < geometry.cols(); v++)
      {
        id[v] = g + v;
        data[v] = 10.8 + g + v;
      }

      vtpUtilities.AddPolygon(geometry,
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Points,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Data",
                                  Gedim::VTPProperty::Formats::Points,
                                  static_cast<unsigned int>(data.size()),
                                  data.data()
                                }
                              });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry2D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test3D)
  {
    const unsigned int numGeometries = 4;


    double angle = M_PI / 4.0; // 45 degrees
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

    Eigen::Vector3d translation(0.0, 2.0, 0.0);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Gedim::GeometryUtilities::Polyhedron geometry = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0 + g,
                                                                                                             0.0 + g,
                                                                                                             0.0 + g),
                                                                                             1.0);

      geometry.Vertices = geometryUtilities.RotatePoints(geometry.Vertices,
                                                         rotationMatrix,
                                                         translation);

      vector<double> id(geometry.Vertices.cols());
      vector<double> data(geometry.Vertices.cols());

      for (unsigned int v = 0; v < geometry.Vertices.cols(); v++)
      {
        id[v] = g + v;
        data[v] = 10.8 + g + v;
      }

      vtpUtilities.AddPolyhedron(geometry.Vertices,
                                 geometry.Edges,
                                 geometry.Faces,
                                 {
                                   {
                                     "Id",
                                     Gedim::VTPProperty::Formats::Points,
                                     static_cast<unsigned int>(id.size()),
                                     id.data()
                                   },
                                   {
                                     "Data",
                                     Gedim::VTPProperty::Formats::Points,
                                     static_cast<unsigned int>(data.size()),
                                     data.data()
                                   }
                                 });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry3D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
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
