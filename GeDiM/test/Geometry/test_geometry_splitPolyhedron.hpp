#ifndef __TEST_GEOMETRY_SPLITPOLYHEDRON_H
#define __TEST_GEOMETRY_SPLITPOLYHEDRON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "VTPUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      geometryUtilityConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // split tetrahedron with plane
      {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                             Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                             Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                                                             Eigen::Vector3d(0.0, 1.0, 0.0));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtility.PolyhedronBarycenter(polyhedron.Vertices);
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtility.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                              polyhedron.Edges);
        const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtility.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                      polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtility.PolyhedronFaceNormals(polyhedronFaceVertices,
                                                                                                    polyhedronBarycenter);
        const vector<vector<bool>> polyhedronFaceEdgeDirections = geometryUtility.PolyhedronFaceEdgeDirections(polyhedron.Vertices,
                                                                                                               polyhedron.Edges,
                                                                                                               polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents = geometryUtility.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                                                                               polyhedron.Edges,
                                                                                                               polyhedron.Faces,
                                                                                                               polyhedronFaceEdgeDirections,
                                                                                                               polyhedronEdgeTangents);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtility.PolyhedronFaceTranslations(polyhedronFaceVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtility.PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                      polyhedronFaceNormals,
                                                                                                                      polyhedronFaceTranslations);

        const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.5);
        const Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
        const Eigen::MatrixXd planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtility.PlaneTranslation(planeNormal,
                                                                                  planeOrigin);

        const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result = geometryUtility.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                                                                                         polyhedron.Edges,
                                                                                                                         polyhedron.Faces,
                                                                                                                         polyhedronFaceVertices,
                                                                                                                         polyhedronFaceEdgesTangents,
                                                                                                                         polyhedronFaceTranslations,
                                                                                                                         polyhedronFaceRotationMatrices,
                                                                                                                         planeNormal,
                                                                                                                         planeOrigin,
                                                                                                                         planeRotationMatrix,
                                                                                                                         planeTranslation);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
        ASSERT_EQ(result.Vertices.Vertices, (Eigen::MatrixXd(3, 7)<<
                                             0,   1,   0,   0, 0.5,   0,   0,
                                             0,   0,   1,   0,   0,   0, 0.5,
                                             0,   0,   0,   1, 0.5, 0.5, 0.5).finished());
        ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({ 4, 3, 5 }));
        ASSERT_EQ(result.Edges.Edges, (Eigen::MatrixXi(2, 12)<<
                                       5, 4, 6, 3, 4, 6, 0, 1, 0, 1, 5, 2,
                                       4, 6, 5, 5, 3, 3, 1, 2, 2, 4, 0, 6).finished());
        ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({ -1,-1,-1,3,4,5,0,2,1,4,3,5 }));
        ASSERT_EQ(result.Faces.Faces.size(), 8);
        ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3)<<
                                          3, 5, 4,
                                          3, 0, 4).finished());
        ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3)<<
                                          3, 5, 6,
                                          3, 2, 5).finished());
        ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3)<<
                                          4, 6, 3,
                                          1, 5, 4).finished());
        ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3)<<
                                          0, 1, 2,
                                          6, 7, 8).finished());
        ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 4)<<
                                          0, 1, 4, 5,
                                          6, 9, 0, 10).finished());
        ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 4)<<
                                          0, 2, 6 , 5,
                                          8, 11, 2, 10).finished());
        ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4)<<
                                          1, 2, 6, 4,
                                          7, 11, 1, 9).finished());
        ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3)<<
                                          5, 4, 6,
                                          0, 1, 2).finished());
        ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({ 1,2,3,0,1,2,3,-1 }));

        ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({ 3,4,5,6 }));
        ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,3,4,5 }));
        ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({ 0,1,2,7 }));

        ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({ 0,1,2,4,5,6 }));
        ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,6,7,8,9,10,11 }));
        ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({ 3,4,5,6,7 }));

        const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra = geometryUtility.SplitPolyhedronWithPlaneResultToPolyhedra(result);

        ASSERT_EQ(splitPolyhedra.size(), 2);

        // Export to VTK
        std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/Tetra";
        Gedim::Output::CreateFolder(exportFolder);

        {
          Gedim::VTKUtilities vtpUtilities;

          //  original polyhedron
          vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                     polyhedron.Edges,
                                     polyhedron.Faces);


          vtpUtilities.Export(exportFolder + "/Original.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices,
              splitPolyhedra[0].Edges,
              splitPolyhedra[0].Faces);


          vtpUtilities.Export(exportFolder + "/Positive.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices,
              splitPolyhedra[1].Edges,
              splitPolyhedra[1].Faces);


          vtpUtilities.Export(exportFolder + "/Negative.vtu",
                              Gedim::VTKUtilities::Ascii);
        }
      }

      // split tetrahedron with plane
      {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                             Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                             Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                                                             Eigen::Vector3d(0.0, 1.0, 0.0));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtility.PolyhedronBarycenter(polyhedron.Vertices);
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtility.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                              polyhedron.Edges);
        const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtility.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                      polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtility.PolyhedronFaceNormals(polyhedronFaceVertices,
                                                                                                    polyhedronBarycenter);
        const vector<vector<bool>> polyhedronFaceEdgeDirections = geometryUtility.PolyhedronFaceEdgeDirections(polyhedron.Vertices,
                                                                                                               polyhedron.Edges,
                                                                                                               polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents = geometryUtility.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                                                                               polyhedron.Edges,
                                                                                                               polyhedron.Faces,
                                                                                                               polyhedronFaceEdgeDirections,
                                                                                                               polyhedronEdgeTangents);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtility.PolyhedronFaceTranslations(polyhedronFaceVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtility.PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                      polyhedronFaceNormals,
                                                                                                                      polyhedronFaceTranslations);

        const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d planeNormal(-1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0));
        const Eigen::MatrixXd planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtility.PlaneTranslation(planeNormal,
                                                                                  planeOrigin);

        const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result = geometryUtility.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                                                                                         polyhedron.Edges,
                                                                                                                         polyhedron.Faces,
                                                                                                                         polyhedronFaceVertices,
                                                                                                                         polyhedronFaceEdgesTangents,
                                                                                                                         polyhedronFaceTranslations,
                                                                                                                         polyhedronFaceRotationMatrices,
                                                                                                                         planeNormal,
                                                                                                                         planeOrigin,
                                                                                                                         planeRotationMatrix,
                                                                                                                         planeTranslation);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
        ASSERT_EQ(result.Vertices.Vertices, (Eigen::MatrixXd(3, 6)<<
                                             0,   1,   0,   0, 0.5,   0,
                                             0,   0,   1,   0,   0, 0.5,
                                             0,   0,   0,   1, 0.5, 0.5).finished());
        ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({ 4, 5 }));
        ASSERT_EQ(result.Edges.Edges, (Eigen::MatrixXi(2, 11)<<
                                       0, 4, 0, 0, 5, 4, 0, 1, 0, 1, 2,
                                       4, 3, 3, 5, 3, 5, 1, 2, 2, 4, 5).finished());
        ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({ -1, 4, 3, -1, 5, -1, 0, 2, 1, 4, 5 }));
        ASSERT_EQ(result.Faces.Faces.size(), 8);
        ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3)<<
                                          0, 4, 3,
                                          0, 1, 2).finished());
        ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3)<<
                                          0, 5, 3,
                                          3, 4, 2).finished());
        ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3)<<
                                          4, 5, 3,
                                          5, 4, 1).finished());
        ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3)<<
                                          0, 1, 2,
                                          6, 7, 8).finished());
        ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 3)<<
                                          0, 1, 4,
                                          6, 9, 0).finished());
        ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 3)<<
                                          0, 2, 5,
                                          8, 10, 3).finished());
        ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4)<<
                                          1, 2, 5, 4,
                                          7, 10, 5, 9).finished());
        ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3)<<
                                          5, 4, 0,
                                          5, 0, 3).finished());
        ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({ 1,2,3,0,1,2,3,-1 }));

        ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({ 0,3,4,5 }));
        ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,3,4,5 }));
        ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({ 0,1,2,7 }));

        ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({ 0,1,2,4,5 }));
        ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({ 0,3,5,6,7,8,9,10 }));
        ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({ 3,4,5,6,7 }));

        const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra = geometryUtility.SplitPolyhedronWithPlaneResultToPolyhedra(result);

        // Export to VTK
        std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/TetraTwo";
        Gedim::Output::CreateFolder(exportFolder);

        {
          Gedim::VTKUtilities vtpUtilities;

          //  original polyhedron
          vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                     polyhedron.Edges,
                                     polyhedron.Faces);


          vtpUtilities.Export(exportFolder + "/Original.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices,
              splitPolyhedra[0].Edges,
              splitPolyhedra[0].Faces);


          vtpUtilities.Export(exportFolder + "/Positive.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices,
              splitPolyhedra[1].Edges,
              splitPolyhedra[1].Faces);


          vtpUtilities.Export(exportFolder + "/Negative.vtu",
                              Gedim::VTKUtilities::Ascii);
        }
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_SPLIT_H
