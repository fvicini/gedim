#ifndef __TEST_MESH_UTILITIES3D_H
#define __TEST_MESH_UTILITIES3D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_3D_22Cells_Mock.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatrices_3D_1Cells_Mock.hpp"
#include "MeshMatrices_3D_68Cells_Mock.hpp"

#include "MeshFromCsvUtilities.hpp"
#include "MeshDAOImporterFromCsv.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestMeshUtilities, TestMesh3DFromPolyhedron)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);

    const vector<unsigned int> vertexMarkers = { 1, 2, 3, 4, 5, 6, 7, 8 };
    const vector<unsigned int> edgeMarkers = { 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
    const vector<unsigned int> faceMarkers = { 21, 22, 23, 24, 25, 26 };

    meshUtilities.Mesh3DFromPolyhedron(cube.Vertices,
                                       cube.Edges,
                                       cube.Faces,
                                       vertexMarkers,
                                       edgeMarkers,
                                       faceMarkers,
                                       meshDao);

    std::string exportFolder = "./Export/TestMesh3DFromPolyhedron/";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    EXPECT_EQ(meshDao.Dimension(), 3);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 8);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 12);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 6);
  }

  TEST(TestMeshUtilities, TestCreateTetrahedralMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                   1.0);

    meshUtilities.CreateTetrahedralMesh(polyhedron.Vertices,
                                        polyhedron.Edges,
                                        polyhedron.Faces,
                                        0.03,
                                        meshDao,
                                        "Qpqfezna");

    std::string exportFolder = "./Export/TestMeshUtilities/TestCreateTetrahedralMesh";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "CreatedTetrahedralMesh");

    EXPECT_EQ(3,
              meshDao.Dimension());
    EXPECT_EQ(28,
              meshDao.Cell0DTotalNumber());
    EXPECT_EQ(103,
              meshDao.Cell1DTotalNumber());
    EXPECT_EQ(127,
              meshDao.Cell2DTotalNumber());
    EXPECT_EQ(51,
              meshDao.Cell3DTotalNumber());
  }

  TEST(TestMeshUtilities, TestCheckMesh3D)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_68Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshUtilities::CheckMesh3DConfiguration config;
    ASSERT_NO_THROW(meshUtilities.CheckMesh3D(config,
                                              geometryUtilities,
                                              meshDao));
  }

  TEST(TestMeshUtilities, TestComputeCell2DCell3DNeighbours)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mesh;
    mesh.Mesh.NumberCell2DNeighbourCell3D.clear();
    mesh.Mesh.Cell2DNeighbourCell3Ds.clear();
    mesh.Mesh.NumberCell2DNeighbourCell3D.resize(mesh.Mesh.NumberCell2D + 1, 0);
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDao);

    EXPECT_EQ(mesh.Mesh.NumberCell2DNeighbourCell3D,
              vector<unsigned int>({ 0,1,2,3,4,5, 6 }));
    EXPECT_EQ(mesh.Mesh.Cell2DNeighbourCell3Ds,
              vector<unsigned int>({ 0,0,0,0,0,0 }));
  }

  TEST(TestMeshUtilities, TestFillMesh3DGeometricData_Convex)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-14;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);

    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::MeshUtilities::MeshGeometricData3D result = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                   meshDao);

    Gedim::MeshUtilities::MeshGeometricData3D expectedResult;
    expectedResult.Cell3DsVolumes = { 1.0 };
    expectedResult.Cell3DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.5) };
    expectedResult.Cell3DsDiameters = { sqrt(3.0) };
    expectedResult.Cell3DsVertices = { cube.Vertices };
    expectedResult.Cell3DsEdges = { cube.Edges };
    expectedResult.Cell3DsFaces = { cube.Faces };

    EXPECT_EQ(result.Cell3DsVertices, expectedResult.Cell3DsVertices);
    EXPECT_EQ(result.Cell3DsEdges, expectedResult.Cell3DsEdges);
    EXPECT_EQ(result.Cell3DsFaces, expectedResult.Cell3DsFaces);
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsVolumes[0], expectedResult.Cell3DsVolumes[0]));
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsCentroids[0].x(), expectedResult.Cell3DsCentroids[0].z()));
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsCentroids[0].y(), expectedResult.Cell3DsCentroids[0].y()));
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsCentroids[0].z(), expectedResult.Cell3DsCentroids[0].x()));
    EXPECT_EQ(result.Cell3DsDiameters, expectedResult.Cell3DsDiameters);
  }

  TEST(TestMeshUtilities, TestSetMeshMarkersOnPlane)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestSetMeshMarkersOnPlane";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh_Original");

    meshUtilities.SetMeshMarkersOnPlane(geometryUtilities,
                                        Eigen::Vector3d { -1.0, 0.0, 0.0 },
                                        Eigen::Vector3d { 0.0, 0.0, 0.0 },
                                        100,
                                        meshDao);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh_Modified");

    ASSERT_EQ(100, meshDao.Cell0DMarker(0));
    ASSERT_EQ(100, meshDao.Cell0DMarker(3));
    ASSERT_EQ(100, meshDao.Cell0DMarker(4));
    ASSERT_EQ(100, meshDao.Cell0DMarker(7));
    ASSERT_EQ(100, meshDao.Cell0DMarker(11));
    ASSERT_EQ(100, meshDao.Cell0DMarker(15));
    ASSERT_EQ(100, meshDao.Cell0DMarker(16));
    ASSERT_EQ(100, meshDao.Cell0DMarker(17));

    ASSERT_EQ(100, meshDao.Cell1DMarker(27));
    ASSERT_EQ(100, meshDao.Cell1DMarker(28));
    ASSERT_EQ(100, meshDao.Cell1DMarker(32));
    ASSERT_EQ(100, meshDao.Cell1DMarker(33));
    ASSERT_EQ(100, meshDao.Cell1DMarker(34));
    ASSERT_EQ(100, meshDao.Cell1DMarker(39));
    ASSERT_EQ(100, meshDao.Cell1DMarker(45));
    ASSERT_EQ(100, meshDao.Cell1DMarker(46));
    ASSERT_EQ(100, meshDao.Cell1DMarker(47));
    ASSERT_EQ(100, meshDao.Cell1DMarker(48));
    ASSERT_EQ(100, meshDao.Cell1DMarker(50));
    ASSERT_EQ(100, meshDao.Cell1DMarker(52));
    ASSERT_EQ(100, meshDao.Cell1DMarker(53));

    ASSERT_EQ(100, meshDao.Cell2DMarker(25));
    ASSERT_EQ(100, meshDao.Cell2DMarker(28));
    ASSERT_EQ(100, meshDao.Cell2DMarker(42));
    ASSERT_EQ(100, meshDao.Cell2DMarker(47));
    ASSERT_EQ(100, meshDao.Cell2DMarker(49));
    ASSERT_EQ(100, meshDao.Cell2DMarker(52));
  }
}

#endif // __TEST_MESH_UTILITIES3D_H