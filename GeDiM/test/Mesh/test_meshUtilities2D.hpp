#ifndef __TEST_MESH_UTILITIES2D_H
#define __TEST_MESH_UTILITIES2D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "CommonUtilities.hpp"

#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatrices_2D_1Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "OVM_Mesh_Mock.hpp"
#include "OpenVolumeMeshInterface.hpp"

#include "MeshFromCsvUtilities.hpp"
#include "MeshDAOImporterFromCsv.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestMeshUtilities, TestMesh2DFromPolygon)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    Eigen::MatrixXd polygon(3, 4);
    polygon.col(0)<< 0.0, 0.0, 0.0;
    polygon.col(1)<< 1.0, 0.0, 0.0;
    polygon.col(2)<< 1.0, 1.0, 0.0;
    polygon.col(3)<< 0.0, 1.0, 0.0;
    vector<unsigned int> vertexMarkers = { 1, 2, 3, 4 };
    vector<unsigned int> edgeMarkers = { 5, 6, 7, 8 };

    meshUtilities.Mesh2DFromPolygon(polygon,
                                    vertexMarkers,
                                    edgeMarkers,
                                    meshDao);

    std::string exportFolder = "./Export/TestMesh2DFromPolygon/";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    EXPECT_EQ(meshDao.Dimension(), 2);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 1);
  }

  TEST(TestMeshUtilities, TestCreateRectangleMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Eigen::Vector3d rectangleOrigin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleBaseTangent = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleHeightTangent = Eigen::Vector3d(0.0, 1.0, 0.0);
    const vector<double> baseCoordinates = geometryUtilities.EquispaceCoordinates(0.25, true);
    const vector<double> heightCoordinates = geometryUtilities.EquispaceCoordinates(0.5, true);

    meshUtilities.CreateRectangleMesh(rectangleOrigin,
                                      rectangleBaseTangent,
                                      rectangleHeightTangent,
                                      baseCoordinates,
                                      heightCoordinates,
                                      meshDao);

    Eigen::MatrixXd cell0DCoordinates(3, 15);
    cell0DCoordinates.col(0)<< 0.0, 0.0, 0.0;
    cell0DCoordinates.col(1)<< 0.25, 0.0, 0.0;
    cell0DCoordinates.col(2)<< 0.5, 0.0, 0.0;
    cell0DCoordinates.col(3)<< 0.75, 0.0, 0.0;
    cell0DCoordinates.col(4)<< 1.0, 0.0, 0.0;
    cell0DCoordinates.col(5)<< 0.0, 0.5, 0.0;
    cell0DCoordinates.col(6)<< 0.25, 0.5, 0.0;
    cell0DCoordinates.col(7)<< 0.5, 0.5, 0.0;
    cell0DCoordinates.col(8)<< 0.75, 0.5, 0.0;
    cell0DCoordinates.col(9)<< 1.0, 0.5, 0.0;
    cell0DCoordinates.col(10)<< 0.0, 1.0, 0.0;
    cell0DCoordinates.col(11)<< 0.25, 1.0, 0.0;
    cell0DCoordinates.col(12)<< 0.5, 1.0, 0.0;
    cell0DCoordinates.col(13)<< 0.75, 1.0, 0.0;
    cell0DCoordinates.col(14)<< 1.0, 1.0, 0.0;

    Eigen::MatrixXi cell1DExtremes(2, 22);
    cell1DExtremes.col(0)<< 0, 1;
    cell1DExtremes.col(1)<< 1, 2;
    cell1DExtremes.col(2)<< 2, 3;
    cell1DExtremes.col(3)<< 3, 4;
    cell1DExtremes.col(4)<< 5, 6;
    cell1DExtremes.col(5)<< 6, 7;
    cell1DExtremes.col(6)<< 7, 8;
    cell1DExtremes.col(7)<< 8, 9;
    cell1DExtremes.col(8)<< 10, 11;
    cell1DExtremes.col(9)<< 11, 12;
    cell1DExtremes.col(10)<< 12, 13;
    cell1DExtremes.col(11)<< 13, 14;
    cell1DExtremes.col(12)<< 0, 5;
    cell1DExtremes.col(13)<< 1, 6;
    cell1DExtremes.col(14)<< 2, 7;
    cell1DExtremes.col(15)<< 3, 8;
    cell1DExtremes.col(16)<< 4, 9;
    cell1DExtremes.col(17)<< 5, 10;
    cell1DExtremes.col(18)<< 6, 11;
    cell1DExtremes.col(19)<< 7, 12;
    cell1DExtremes.col(20)<< 8, 13;
    cell1DExtremes.col(21)<< 9, 14;

    EXPECT_EQ(meshDao.Dimension(), 2);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 15);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 22);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 8);
    EXPECT_EQ(meshDao.Cell0DsCoordinates(),
              cell0DCoordinates);
    EXPECT_EQ(meshDao.Cell1DsExtremes(),
              cell1DExtremes);
  }

  TEST(TestMeshUtilities, TestFillMesh2D_TriangularMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO expectedMesh(mockMesh.Mesh);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Eigen::MatrixXd vertices = (Eigen::MatrixXd(3, 5)<<
                                      0.0, 1.0, 1.0, 0.0, 0.5,
                                      0.0, 0.0, 1.0, 1.0, 0.5,
                                      0.0, 0.0, 0.0, 0.0, 0.0).finished();
    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 8)<<
                                   1, 2, 4, 3, 0, 4, 2, 0,
                                   2, 4, 1, 0, 4, 3, 3, 1).finished();

    vector<Eigen::MatrixXi> triangles(4, Eigen::MatrixXi(2, 3));
    triangles[0] = (Eigen::MatrixXi(2, 3)<<
                    1,2,4,
                    0,1,2).finished();
    triangles[1] = (Eigen::MatrixXi(2, 3)<<
                    3,0,4,
                    3,4,5).finished();
    triangles[2] = (Eigen::MatrixXi(2, 3)<<
                    4,2,3,
                    1,6,5).finished();
    triangles[3] = (Eigen::MatrixXi(2, 3)<<
                    0,1,4,
                    7,2,4).finished();

    meshUtilities.FillMesh2D(vertices,
                             edges,
                             triangles,
                             meshDao);

    ASSERT_EQ(expectedMesh.Cell0DsCoordinates(),
              meshDao.Cell0DsCoordinates());
    ASSERT_EQ(expectedMesh.Cell1DsExtremes(),
              meshDao.Cell1DsExtremes());
    for (unsigned int c = 0; c < 4; c++)
    {
      ASSERT_EQ(expectedMesh.Cell2DVertices(c),
                meshDao.Cell2DVertices(c));
      ASSERT_EQ(expectedMesh.Cell2DEdges(c),
                meshDao.Cell2DEdges(c));
    }

  }

  TEST(TestMeshUtilities, TestCreateTriangleMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO expectedMesh(mockMesh.Mesh);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Eigen::MatrixXd polygon = geometryUtilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                   1.0);

    meshUtilities.CreateTriangularMesh(polygon,
                                       0.033,
                                       meshDao,
                                       "-QDzpqnea");

    std::string exportFolder = "./Export/TestMeshUtilities/TestCreateTriangleMesh";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "CreatedTriangleMesh");
    meshUtilities.ExportMeshToVTU(expectedMesh,
                                  exportFolder,
                                  "ExpectedTriangleMesh");

    Gedim::MeshUtilities::MeshGeometricData2D cell2DsGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                           meshDao);
    const unsigned int cell2DToExportIndex = 0;
    meshUtilities.ExportCell2DToVTU(meshDao,
                                    cell2DToExportIndex,
                                    cell2DsGeometricData.Cell2DsVertices[cell2DToExportIndex],
                                    cell2DsGeometricData.Cell2DsTriangulations[cell2DToExportIndex],
                                    cell2DsGeometricData.Cell2DsAreas[cell2DToExportIndex],
                                    cell2DsGeometricData.Cell2DsCentroids[cell2DToExportIndex],
                                    exportFolder);


    EXPECT_EQ(expectedMesh.Dimension(),
              meshDao.Dimension());
    EXPECT_EQ(expectedMesh.Cell0DTotalNumber(),
              meshDao.Cell0DTotalNumber());
    EXPECT_EQ(expectedMesh.Cell1DTotalNumber(),
              meshDao.Cell1DTotalNumber());
    EXPECT_EQ(expectedMesh.Cell2DTotalNumber(),
              meshDao.Cell2DTotalNumber());
    EXPECT_EQ(expectedMesh.Cell0DsCoordinates(),
              meshDao.Cell0DsCoordinates());
    EXPECT_EQ(expectedMesh.Cell1DsExtremes(),
              meshDao.Cell1DsExtremes());

    const vector<unsigned int> cell0DMarkers = { 11, 12, 13, 14 };
    const vector<unsigned int> cell1DMarkers = { 15, 16, 17, 18 };
    meshUtilities.ChangePolygonMeshMarkers(polygon,
                                           cell0DMarkers,
                                           cell1DMarkers,
                                           meshDao);

    EXPECT_EQ(cell0DMarkers[0],
        meshDao.Cell0DMarker(0));
    EXPECT_EQ(cell0DMarkers[1],
        meshDao.Cell0DMarker(1));
    EXPECT_EQ(cell0DMarkers[2],
        meshDao.Cell0DMarker(2));
    EXPECT_EQ(cell0DMarkers[3],
        meshDao.Cell0DMarker(3));
    EXPECT_EQ(cell1DMarkers[1],
        meshDao.Cell1DMarker(20));
    EXPECT_EQ(cell1DMarkers[0],
        meshDao.Cell1DMarker(23));
  }

  TEST(TestMeshUtilities, TestCheckMesh2D)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    ASSERT_NO_THROW(meshUtilities.CheckMesh2D(config,
                                              geometryUtilities,
                                              meshDao));
  }

  TEST(TestMeshUtilities, TestComputeCell1DCell2DNeighbours)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    mesh.Mesh.NumberCell1DNeighbourCell2D.clear();
    mesh.Mesh.Cell1DNeighbourCell2Ds.clear();
    mesh.Mesh.NumberCell1DNeighbourCell2D.resize(mesh.Mesh.NumberCell1D + 1, 0);
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    meshUtilities.ComputeCell1DCell2DNeighbours(meshDao);

    EXPECT_EQ(mesh.Mesh.NumberCell1DNeighbourCell2D,
              vector<unsigned int>({ 0,2,4,6,8 }));
    EXPECT_EQ(mesh.Mesh.Cell1DNeighbourCell2Ds,
              vector<unsigned int>({ 1,0,1,0,1,0,1,0 }));
  }

  TEST(TestMeshUtilities, TestFillMesh2DGeometricData_NonConvex)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    GedimUnitTesting::MeshMatrices_2D_2Cells_Mock convexMesh;
    Gedim::MeshMatricesDAO convexMeshDao(convexMesh.Mesh);

    vector<vector<unsigned int>> meshCell2DToConvexCell2DIndices = { { 0, 1 } };

    Gedim::MeshUtilities meshUtilities;

    const Gedim::MeshUtilities::MeshGeometricData2D result = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                   meshDao,
                                                                                                   convexMeshDao,
                                                                                                   meshCell2DToConvexCell2DIndices);

    Gedim::MeshUtilities::MeshGeometricData2D expectedResult;
    expectedResult.Cell2DsAreas = { 1.0 };
    expectedResult.Cell2DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.0) };
    expectedResult.Cell2DsDiameters = { sqrt(2.0) };
    expectedResult.Cell2DsEdgeDirections = { { true, true, true, true } };
    Eigen::VectorXd edgeLengths(4);
    edgeLengths<< 1.0, 1.0, 1.0, 1.0;
    expectedResult.Cell2DsEdgeLengths = { edgeLengths };
    Eigen::MatrixXd edgeNormals(3, 4);
    edgeNormals.col(0)<< -1.0, 0.0, 0.0;
    edgeNormals.col(1)<< 0.0, -1.0, 0.0;
    edgeNormals.col(2)<< 1.0, 0.0, 0.0;
    edgeNormals.col(3)<< 0.0, 1.0, 0.0;
    expectedResult.Cell2DsEdgeNormals = { edgeNormals };
    Eigen::MatrixXd edgeTangents(3, 4);
    edgeTangents.col(0)<< 0.0, -1.0, 0.0;
    edgeTangents.col(1)<< 1.0, 0.0, 0.0;
    edgeTangents.col(2)<< 0.0, 1.0, 0.0;
    edgeTangents.col(3)<< -1.0, 0.0, 0.0;
    expectedResult.Cell2DsEdgeTangents = { edgeTangents };
    Eigen::Matrix3d triangleOne;
    triangleOne.col(0)<< 0.0, 1.0, 0.0;
    triangleOne.col(1)<< 0.0, 0.0, 0.0;
    triangleOne.col(2)<< 1.0, 0.0, 0.0;
    Eigen::Matrix3d triangleTwo;
    triangleTwo.col(0)<< 1.0, 0.0, 0.0;
    triangleTwo.col(1)<< 1.0, 1.0, 0.0;
    triangleTwo.col(2)<< 0.0, 1.0, 0.0;
    expectedResult.Cell2DsTriangulations = { { triangleOne, triangleTwo } };
    Eigen::MatrixXd vertices(3, 4);
    vertices.col(0)<< 0.0, 1.0, 0.0;
    vertices.col(1)<< 0.0, 0.0, 0.0;
    vertices.col(2)<< 1.0, 0.0, 0.0;
    vertices.col(3)<< 1.0, 1.0, 0.0;
    expectedResult.Cell2DsVertices = { vertices };

    EXPECT_EQ(result.Cell2DsAreas, expectedResult.Cell2DsAreas);
    EXPECT_EQ(result.Cell2DsCentroids, expectedResult.Cell2DsCentroids);
    EXPECT_EQ(result.Cell2DsDiameters, expectedResult.Cell2DsDiameters);
    EXPECT_EQ(result.Cell2DsEdgeDirections, expectedResult.Cell2DsEdgeDirections);
    EXPECT_EQ(result.Cell2DsEdgeLengths, expectedResult.Cell2DsEdgeLengths);
    EXPECT_EQ(result.Cell2DsEdgeNormals, expectedResult.Cell2DsEdgeNormals);
    EXPECT_EQ(result.Cell2DsEdgeTangents, expectedResult.Cell2DsEdgeTangents);
    EXPECT_EQ(result.Cell2DsTriangulations, expectedResult.Cell2DsTriangulations);
    EXPECT_EQ(result.Cell2DsVertices, expectedResult.Cell2DsVertices);
  }

  TEST(TestMeshUtilities, TestRefineMesh2D)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestRefineMesh2D/";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-14;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ComputeCell1DCell2DNeighbours(meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_R0");

    // first refine
    {
      const unsigned int newCell0DIndex = meshDAO.Cell0DAppend(1);
      meshDAO.Cell0DInsertCoordinates(newCell0DIndex, Eigen::Vector3d(0.5, 1.0, 0.0));
      meshDAO.Cell0DSetMarker(newCell0DIndex, meshDAO.Cell1DMarker(3));
      meshDAO.Cell0DSetState(newCell0DIndex, true);

      EXPECT_EQ(4, newCell0DIndex);

      Eigen::MatrixXi newCell1DsExtreme(2, 2);
      newCell1DsExtreme.col(0)<< 2, 4;
      newCell1DsExtreme.col(1)<< 4, 3;
      const std::vector<unsigned int> newCell1DsIndex = meshUtilities.SplitCell1D(3,
                                                                                  newCell1DsExtreme,
                                                                                  meshDAO);

      EXPECT_EQ(std::vector<unsigned int>({ 4, 5 }), newCell1DsIndex);
      EXPECT_FALSE(meshDAO.Cell1DIsActive(3));
      EXPECT_TRUE(meshDAO.Cell1DIsActive(4));
      EXPECT_TRUE(meshDAO.Cell1DIsActive(5));
      EXPECT_EQ(meshDAO.Cell1DMarker(3), meshDAO.Cell1DMarker(4));
      EXPECT_EQ(meshDAO.Cell1DMarker(3), meshDAO.Cell1DMarker(5));
      EXPECT_EQ(2, meshDAO.Cell1DNumberNeighbourCell2D(4));
      EXPECT_FALSE(meshDAO.Cell1DHasNeighbourCell2D(4, 0));
      EXPECT_TRUE(meshDAO.Cell1DHasNeighbourCell2D(4, 1));
      EXPECT_EQ(0, meshDAO.Cell1DNeighbourCell2D(4, 1));
      EXPECT_EQ(2, meshDAO.Cell1DNumberNeighbourCell2D(5));
      EXPECT_FALSE(meshDAO.Cell1DHasNeighbourCell2D(5, 0));
      EXPECT_TRUE(meshDAO.Cell1DHasNeighbourCell2D(5, 1));
      EXPECT_EQ(0, meshDAO.Cell1DNeighbourCell2D(5, 1));
      EXPECT_EQ(3, meshDAO.Cell1DOriginalCell1D(4));
      EXPECT_EQ(3, meshDAO.Cell1DOriginalCell1D(5));
      std::list<unsigned int> updatedCell1Ds;
      EXPECT_FALSE(meshDAO.Cell1DUpdatedCell1Ds(3,
                                                updatedCell1Ds));
      updatedCell1Ds.sort();
      EXPECT_EQ(std::list<unsigned int>({4, 5}), updatedCell1Ds);

      Eigen::MatrixXi subCell(2, 5);
      subCell.row(0)<< 0, 1, 2, 4, 3;
      subCell.row(1)<< 1, 2, 4, 5, 0;

      const std::vector<unsigned int> newCell2DsIndex = meshUtilities.SplitCell2D(0,
                                                                                  { subCell },
                                                                                  meshDAO);
      EXPECT_EQ(std::vector<unsigned int>({ 1 }), newCell2DsIndex);
      EXPECT_FALSE(meshDAO.Cell2DIsActive(0));
      EXPECT_TRUE(meshDAO.Cell2DIsActive(1));
      EXPECT_EQ(meshDAO.Cell2DMarker(0), meshDAO.Cell2DMarker(1));
      std::list<unsigned int> updatedCell2Ds;
      EXPECT_FALSE(meshDAO.Cell2DUpdatedCell2Ds(0,
                                                updatedCell2Ds));
      EXPECT_EQ(0, meshDAO.Cell2DOriginalCell2D(1));
      EXPECT_EQ(std::list<unsigned int>({1}), updatedCell2Ds);
      EXPECT_EQ(1, meshDAO.Cell1DNeighbourCell2D(0, 1));
      EXPECT_EQ(1, meshDAO.Cell1DNeighbourCell2D(1, 1));
      EXPECT_EQ(1, meshDAO.Cell1DNeighbourCell2D(2, 1));
      EXPECT_EQ(1, meshDAO.Cell1DNeighbourCell2D(4, 1));
      EXPECT_EQ(1, meshDAO.Cell1DNeighbourCell2D(5, 1));

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R1");
    }

    // second refine
    {
      const unsigned int newCell1DIndex = meshDAO.Cell1DAppend(1);
      meshDAO.Cell1DInsertExtremes(newCell1DIndex, 4, 0);
      meshDAO.Cell1DSetMarker(newCell1DIndex, 0);
      meshDAO.Cell1DSetState(newCell1DIndex, true);
      meshDAO.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);
      EXPECT_EQ(6, newCell1DIndex);
      EXPECT_TRUE(meshDAO.Cell1DIsActive(6));
      EXPECT_EQ(0, meshDAO.Cell1DMarker(6));
      EXPECT_EQ(2, meshDAO.Cell1DNumberNeighbourCell2D(6));

      vector<Eigen::MatrixXi> subCells(2);

      subCells[0].resize(2, 3);
      subCells[0].row(0)<< 0, 4, 3;
      subCells[0].row(1)<< 6, 5, 0;

      subCells[1].resize(2, 4);
      subCells[1].row(0)<< 0, 1, 2, 4;
      subCells[1].row(1)<< 1, 2, 4, 6;

      const std::vector<unsigned int> newCell2DsIndex = meshUtilities.SplitCell2D(1,
                                                                                  subCells,
                                                                                  meshDAO);
      EXPECT_EQ(std::vector<unsigned int>({ 2, 3 }), newCell2DsIndex);
      EXPECT_FALSE(meshDAO.Cell2DIsActive(1));
      EXPECT_TRUE(meshDAO.Cell2DIsActive(2));
      EXPECT_TRUE(meshDAO.Cell2DIsActive(3));
      EXPECT_EQ(meshDAO.Cell2DMarker(1), meshDAO.Cell2DMarker(2));
      EXPECT_EQ(meshDAO.Cell2DMarker(1), meshDAO.Cell2DMarker(3));
      std::list<unsigned int> updatedCell2Ds;
      EXPECT_FALSE(meshDAO.Cell2DUpdatedCell2Ds(1,
                                                updatedCell2Ds));
      updatedCell2Ds.sort();
      EXPECT_EQ(1, meshDAO.Cell2DOriginalCell2D(2));
      EXPECT_EQ(1, meshDAO.Cell2DOriginalCell2D(3));
      EXPECT_EQ(std::list<unsigned int>({2, 3}), updatedCell2Ds);
      EXPECT_EQ(2, meshDAO.Cell1DNeighbourCell2D(0, 1));
      EXPECT_EQ(3, meshDAO.Cell1DNeighbourCell2D(1, 1));
      EXPECT_EQ(3, meshDAO.Cell1DNeighbourCell2D(2, 1));
      EXPECT_EQ(3, meshDAO.Cell1DNeighbourCell2D(4, 1));
      EXPECT_EQ(2, meshDAO.Cell1DNeighbourCell2D(5, 1));

      meshDAO.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                          0,
                                          2); // right
      meshDAO.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                          1,
                                          3); // left
      EXPECT_EQ(2, meshDAO.Cell1DNeighbourCell2D(6, 0));
      EXPECT_EQ(3, meshDAO.Cell1DNeighbourCell2D(6, 1));

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R2");
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Final");

  }

  TEST(TestMeshUtilities, TestImportOVMMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    std::vector<std::vector<bool>> meshCell3DsFacesOrientation;
    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::OpenVolumeMeshInterface ovmInterface;

    const std::vector<string> lines = OVM_Mesh_Mock::FileLines();
    const Gedim::OpenVolumeMeshInterface::OVMMesh ovm_mesh = ovmInterface.StringsToOVMMesh(lines);
    ovmInterface.OVMMeshToMeshDAO(ovm_mesh,
                                  meshDao,
                                  meshCell3DsFacesOrientation);

    std::string exportFolder = "./Export/TestMeshUtilities/TestImportOVMMesh";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ImportedOVMMesh");

    ASSERT_EQ(3,
              meshDao.Dimension());
    ASSERT_EQ(64,
              meshDao.Cell0DTotalNumber());
    ASSERT_EQ(144,
              meshDao.Cell1DTotalNumber());
    ASSERT_EQ(108,
              meshDao.Cell2DTotalNumber());
    ASSERT_EQ(27,
              meshDao.Cell3DTotalNumber());

    const Gedim::OpenVolumeMeshInterface::OVMMesh reconverted_ovm_mesh = ovmInterface.MeshDAOToOVMMesh(meshDao,
                                                                                                       meshCell3DsFacesOrientation);

    ASSERT_EQ(ovm_mesh.NumCell0Ds, reconverted_ovm_mesh.NumCell0Ds);
    ASSERT_EQ(ovm_mesh.NumCell1Ds, reconverted_ovm_mesh.NumCell1Ds);
    ASSERT_EQ(ovm_mesh.NumCell2Ds, reconverted_ovm_mesh.NumCell2Ds);
    ASSERT_EQ(ovm_mesh.NumCell3Ds, reconverted_ovm_mesh.NumCell3Ds);
    ASSERT_EQ(ovm_mesh.Cell0Ds, reconverted_ovm_mesh.Cell0Ds);
    ASSERT_EQ(ovm_mesh.Cell1Ds, reconverted_ovm_mesh.Cell1Ds);
    ASSERT_EQ(ovm_mesh.Cell2Ds, reconverted_ovm_mesh.Cell2Ds);
    for (unsigned int p = 0; p < ovm_mesh.NumCell3Ds; p++)
    {
      ASSERT_EQ(ovm_mesh.Cell3Ds[p].FacesIndex, reconverted_ovm_mesh.Cell3Ds[p].FacesIndex);
      ASSERT_EQ(ovm_mesh.Cell3Ds[p].FacesOrientation, reconverted_ovm_mesh.Cell3Ds[p].FacesOrientation);
    }

    const std::vector<std::string> reconverted_lines = ovmInterface.OVMMeshToStrings(reconverted_ovm_mesh);

    ASSERT_EQ(lines.size(), reconverted_lines.size());
    for (unsigned int l = 0; l < lines.size(); l++)
      ASSERT_EQ(lines[l], reconverted_lines[l]);

    meshUtilities.ExportMeshToOpenVolume(meshDao,
                                         meshCell3DsFacesOrientation,
                                         exportFolder + "/hex_parallel_1.ovm");
  }

  TEST(TestMeshUtilities, TestFindCell2DsCommonVertices)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    std::string exportFolder = "./Export/TestMeshUtilities/TestFindCell2DsCommonVertices";
    Gedim::Output::CreateFolder(exportFolder);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConvexMesh");

    ASSERT_EQ(meshUtilities.FindCell2DsCommonVertices({0, 23, 29, 13, 30},
                                                      meshDao),
              std::vector<unsigned int>({ 11 }));
    ASSERT_EQ(meshUtilities.FindCell2DsCommonVertices({0, 23},
                                                      meshDao),
              std::vector<unsigned int>({ 11, 18 }));
    ASSERT_EQ(meshUtilities.FindCell2DsCommonVertices({ 30, 6 },
                                                      meshDao),
              std::vector<unsigned int>({ }));
  }

  TEST(TestMeshUtilities, TestFindCell2DsCommonEdges)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    std::string exportFolder = "./Export/TestMeshUtilities/TestFindCell2DsCommonEdges";
    Gedim::Output::CreateFolder(exportFolder);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConvexMesh");

    ASSERT_EQ(meshUtilities.FindCell2DsCommonEdges({0, 23, 29, 13, 30},
                                                   meshDao),
              std::vector<unsigned int>({  }));
    ASSERT_EQ(meshUtilities.FindCell2DsCommonEdges({0, 23},
                                                   meshDao),
              std::vector<unsigned int>({ 0 }));
    ASSERT_EQ(meshUtilities.FindCell2DsCommonEdges({ 30, 6 },
                                                   meshDao),
              std::vector<unsigned int>({ }));
  }

  TEST(TestMeshUtilities, TestAgglomerateTriangles)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    std::string exportFolder = "./Export/TestMeshUtilities/TestAgglomerateTriangles";
    Gedim::Output::CreateFolder(exportFolder);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConvexMesh");


    const  Gedim::MeshUtilities::AgglomerateTrianglesResult result = meshUtilities.AgglomerateTriangles({0, 23, 29, 13, 30},
                                                                                                        meshDao);

    ASSERT_EQ(vector<unsigned int>({ 11, 4, 18, 7, 22, 2, 23 }),
              result.VerticesIndex);
    ASSERT_EQ(vector<unsigned int>({ 1, 2, 51, 54, 34, 55, 33 }),
              result.EdgesIndex);
    ASSERT_EQ(vector<unsigned int>({ 0, 50, 36, 35 }),
              result.RemovedEdges);
  }

  TEST(TestMeshUtilities, TestAgglomerateMeshFromTriangularMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    std::string exportFolder = "./Export/TestMeshUtilities/TestAgglomerateMeshFromTriangularMesh";
    Gedim::Output::CreateFolder(exportFolder);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConvexMesh");

    const std::vector<std::vector<unsigned int>> trianglesToAgglomerate { {4, 31, 32},
                                                                          {0, 23, 29, 13, 30}
                                                                        };

    const Gedim::MeshUtilities::AgglomerateMeshFromTriangularMeshResult result = meshUtilities.AgglomerateMeshFromTriangularMesh(trianglesToAgglomerate,
                                                                                                                                 meshDao);

    ASSERT_EQ(vector<unsigned int>({ 0, 13, 35, 36, 50, 56 }),
              result.RemovedCell1Ds);
    ASSERT_EQ(vector<unsigned int>({ 0, 4, 13, 23, 29, 30, 31, 32 }),
              result.RemovedCell2Ds);
    ASSERT_EQ(2,
              result.ConcaveCell2Ds.size());
    ASSERT_EQ(34,
              result.ConcaveCell2Ds[0].Cell2DIndex);
    ASSERT_EQ(vector<unsigned int>({ 4, 31, 32 }),
              result.ConcaveCell2Ds[0].ConvexCell2DsIndex);
    ASSERT_EQ(35,
              result.ConcaveCell2Ds[1].Cell2DIndex);
    ASSERT_EQ(vector<unsigned int>({ 0, 23, 29, 13, 30 }),
              result.ConcaveCell2Ds[1].ConvexCell2DsIndex);

    std::vector<std::vector<unsigned int>> convexCell2DsIndex(meshDao.Cell2DTotalNumber());
    for (unsigned int c = 0; c < mesh.Mesh.NumberCell2D; c++)
    {
      if (meshDao.Cell2DIsActive(c))
        convexCell2DsIndex[c].resize(1, c);
    }
    for (unsigned int cc = 0; cc < result.ConcaveCell2Ds.size(); cc++)
    {
      const Gedim::MeshUtilities::AgglomerateMeshFromTriangularMeshResult::ConcaveCell2D& concaveCell = result.ConcaveCell2Ds[cc];
      convexCell2DsIndex[concaveCell.Cell2DIndex] = concaveCell.ConvexCell2DsIndex;
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDao,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConcaveMesh");

    std::vector<std::vector<unsigned int>> extractedConvexCell2DsIndex(meshDao.Cell2DTotalNumber());

    for (unsigned int c = 0; c < meshDao.Cell2DTotalNumber(); c++)
    {
      const unsigned int oldCell2DIndex = extractionData.NewCell2DToOldCell2D.at(c);
      extractedConvexCell2DsIndex[c] = convexCell2DsIndex[oldCell2DIndex];
    }

    const string exportMeshFolder = exportFolder + "/Mesh";
    Gedim::Output::CreateFolder(exportMeshFolder);

    meshUtilities.ExportConcaveMesh2DToCsv(meshDao,
                                           extractedConvexCell2DsIndex,
                                           ',',
                                           exportMeshFolder);
  }

  TEST(TestMeshUtilities, TEMP)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    std::string importFolder = "/home/geoscore/Dropbox/Polito/Articles/VEM_INERTIA/VEM_INERTIA/MESH/TriangularMesh/Convex";

    std::string exportFolder = "./Export/TestMeshUtilities/TEMP";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshFromCsvUtilities importerUtilities;
    Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
    meshImporterConfiguration.Folder = importFolder;
    meshImporterConfiguration.Separator = ';';
    Gedim::MeshDAOImporterFromCsv importer(importerUtilities);
    importer.Import(meshImporterConfiguration,
                    meshDao);

    meshUtilities.ComputeCell1DCell2DNeighbours(meshDao);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConvexMesh");

    const unsigned int convexCell2DNumber = meshDao.Cell2DTotalNumber();

    // generate randomly the triangles
    const std::vector<unsigned int> randomCell2Ds = Gedim::Utilities::RandomArrayNoRepetition((0.25 * convexCell2DNumber + 1),
                                                                                              convexCell2DNumber);

    const std::vector<double> percentagesToTake = { 0.25, 0.5, 0.75  };
    std::list<std::vector<unsigned int>> takenTrianglesToAgglomerate;
    std::vector<bool> used_triangles(convexCell2DNumber, 0);
    for (unsigned int t = 0; t < randomCell2Ds.size(); t++)
    {
      const unsigned int cell2DIndex = randomCell2Ds.at(t);
      const unsigned int commonCell0DIndex = meshDao.Cell2DVertex(cell2DIndex, 0);

      std::list<unsigned int> tianglesToAgglomerateList;
      int cell2DNeigh = cell2DIndex;
      do
      {
        if (used_triangles[cell2DNeigh])
          break;

        tianglesToAgglomerateList.push_back(cell2DNeigh);

        const unsigned int commonCell0DLocalPosition = meshDao.Cell2DFindVertex(cell2DNeigh,
                                                                                commonCell0DIndex);
        const unsigned int commonEdgeLocalIndex = (commonCell0DLocalPosition + 2) % 3;
        const unsigned int commonCell1DIndex = meshDao.Cell2DEdge(cell2DNeigh,
                                                                  commonEdgeLocalIndex);

        unsigned int cell2DOtherNeigh = convexCell2DNumber;
        for (unsigned int n = 0; n < meshDao.Cell1DNumberNeighbourCell2D(commonCell1DIndex); n++)
        {
          if (!meshDao.Cell1DHasNeighbourCell2D(commonCell1DIndex, n))
            continue;

          const unsigned int neigh = meshDao.Cell1DNeighbourCell2D(commonCell1DIndex, n);
          if (neigh == cell2DNeigh || neigh == cell2DIndex)
            continue;

          cell2DOtherNeigh = neigh;
        }

        if (cell2DOtherNeigh == convexCell2DNumber)
          cell2DNeigh = convexCell2DNumber;
        else
          cell2DNeigh = cell2DOtherNeigh;
      }
      while (cell2DNeigh < convexCell2DNumber);

      if (tianglesToAgglomerateList.size() < 2)
        continue;

      unsigned int sizeToTake = (percentagesToTake[t % 3] * tianglesToAgglomerateList.size() + 1);
      sizeToTake = (sizeToTake < 2 || sizeToTake > tianglesToAgglomerateList.size()) ?
                     tianglesToAgglomerateList.size() - 1 :
                     sizeToTake;

      takenTrianglesToAgglomerate.push_back(std::vector<unsigned int>(sizeToTake));

      unsigned int iterator = 0;
      for (const unsigned int triangle : tianglesToAgglomerateList)
      {
        if (iterator >= sizeToTake)
          break;

        used_triangles[triangle] = true;
        takenTrianglesToAgglomerate.back()[iterator++] = triangle;
      }
    }

    const std::vector<std::vector<unsigned int>> trianglesToAgglomerate =
        std::vector<std::vector<unsigned int>>(takenTrianglesToAgglomerate.begin(),
                                               takenTrianglesToAgglomerate.end());

    const Gedim::MeshUtilities::AgglomerateMeshFromTriangularMeshResult result = meshUtilities.AgglomerateMeshFromTriangularMesh(trianglesToAgglomerate,
                                                                                                                                 meshDao);

    std::vector<std::vector<unsigned int>> convexCell2DsIndex(meshDao.Cell2DTotalNumber());
    for (unsigned int c = 0; c < convexCell2DNumber; c++)
    {
      if (meshDao.Cell2DIsActive(c))
        convexCell2DsIndex[c].resize(1, c);
    }
    for (unsigned int cc = 0; cc < result.ConcaveCell2Ds.size(); cc++)
    {
      const Gedim::MeshUtilities::AgglomerateMeshFromTriangularMeshResult::ConcaveCell2D& concaveCell = result.ConcaveCell2Ds[cc];
      convexCell2DsIndex[concaveCell.Cell2DIndex] = concaveCell.ConvexCell2DsIndex;
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDao,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConcaveMesh");

    std::vector<std::vector<unsigned int>> extractedConvexCell2DsIndex(meshDao.Cell2DTotalNumber());

    for (unsigned int c = 0; c < meshDao.Cell2DTotalNumber(); c++)
    {
      const unsigned int oldCell2DIndex = extractionData.NewCell2DToOldCell2D.at(c);
      extractedConvexCell2DsIndex[c] = convexCell2DsIndex[oldCell2DIndex];
    }

    const string exportMeshFolder = exportFolder + "/Mesh";
    Gedim::Output::CreateFolder(exportMeshFolder);

    meshUtilities.ExportConcaveMesh2DToCsv(meshDao,
                                           extractedConvexCell2DsIndex,
                                           ',',
                                           exportMeshFolder);
  }
}

#endif // __TEST_MESH_UTILITIES2D_H