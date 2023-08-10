#include "FileTextReader.hpp"
#include "MeshDAOExporterToCsv.hpp"
#include "MeshUtilities.hpp"

#include "TriangleInterface.hpp"
#include "VTKUtilities.hpp"
#include <numeric>
#include <queue>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  void MeshUtilities::FillMesh2D(const MatrixXd& cell0Ds,
                                 const MatrixXi& cell1Ds,
                                 const vector<MatrixXi>& cell2Ds,
                                 IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(2);

    // Create Cell0Ds
    Output::Assert(cell0Ds.rows() == 3);
    const unsigned int& numCell0Ds = cell0Ds.cols();
    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell0DsInsertCoordinates(cell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
      mesh.Cell0DSetState(v, true);

    // Create Cell1Ds
    Output::Assert(cell1Ds.rows() == 2);
    unsigned int numCell1Ds = cell1Ds.cols();
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell1DsInsertExtremes(cell1Ds);
    for (int e = 0; e < cell1Ds.cols(); e++)
      mesh.Cell1DSetState(e, true);

    // Create Cell2Ds
    const unsigned int& numCell2Ds = cell2Ds.size();
    mesh.Cell2DsInitialize(numCell2Ds);
    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const MatrixXi& polygon = cell2Ds[f];
      Output::Assert(polygon.rows() == 2);
      const unsigned int& numVertices = polygon.cols();

      mesh.Cell2DInitializeVertices(f, numVertices);
      mesh.Cell2DInitializeEdges(f, numVertices);

      for (unsigned int v = 0; v < numVertices; v++)
        mesh.Cell2DInsertVertex(f, v, polygon(0, v));
      for (unsigned int e = 0; e < numVertices; e++)
        mesh.Cell2DInsertEdge(f, e, polygon(1, e));

      mesh.Cell2DSetState(f, true);
    }
  }
  // ***************************************************************************
  MeshUtilities::ComputeMesh2DCell1DsResult MeshUtilities::ComputeMesh2DCell1Ds(const Eigen::MatrixXd& cell0Ds,
                                                                                const std::vector<Eigen::VectorXi>& cell2Ds) const
  {
    ComputeMesh2DCell1DsResult result;

    const unsigned int numVertices = cell0Ds.cols();
    const unsigned int numCell2Ds = cell2Ds.size();

    Eigen::SparseMatrix<unsigned int> edges;
    edges.resize(numVertices,
                 numVertices);

    std::list<Eigen::Triplet<unsigned int>> triplets;
    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
      const Eigen::VectorXi& cell2DVertices = cell2Ds.at(c);
      const unsigned int& numCell2DVertices = cell2DVertices.size();

      for (unsigned int v = 0; v < numCell2DVertices; v++)
      {
        const unsigned int origin = cell2DVertices[v];
        const unsigned int end = cell2DVertices[(v + 1) % numCell2DVertices];
        triplets.push_back(Eigen::Triplet<unsigned int>(origin, end, 1));
        triplets.push_back(Eigen::Triplet<unsigned int>(end, origin, 1));
      }
    }

    edges.setFromTriplets(triplets.begin(), triplets.end());
    edges.makeCompressed();

    unsigned int numEdges = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
      {
        if (it.row() < it.col())
          it.valueRef() = 1 + numEdges++;
      }
    }

    result.Cell1Ds.resize(2, numEdges);

    numEdges = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
      {
        if (it.row() < it.col())
          result.Cell1Ds.col(numEdges++)<< it.row(), it.col();
      }
    }

    result.Cell2Ds.resize(numCell2Ds);

    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
      Eigen::MatrixXi& cell2D = result.Cell2Ds.at(c);
      const Eigen::VectorXi& cell2DVertices = cell2Ds.at(c);
      const unsigned int& numCell2DVertices = cell2DVertices.size();

      cell2D.resize(2, numCell2DVertices);

      for (unsigned int v = 0; v < numCell2DVertices; v++)
      {
        const unsigned int origin = cell2DVertices[v];
        const unsigned int end = cell2DVertices[(v + 1) % numCell2DVertices];

        cell2D(0, v) = origin;
        cell2D(1, v) = (origin < end) ? edges.coeff(origin, end) - 1 : edges.coeff(end, origin) - 1;
      }
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::CheckMesh2D(const CheckMesh2DConfiguration& configuration,
                                  const GeometryUtilities& geometryUtilities,
                                  const IMeshDAO& convexMesh) const
  {
    Output::Assert(convexMesh.Dimension() == 2);

    // check Cell0D are 2D
    if (configuration.Cell0D_CheckCoordinates2D)
      Output::Assert(geometryUtilities.PointsAre2D(convexMesh.Cell0DsCoordinates()));

    // check Cell0D duplications
    if (configuration.Cell0D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell0DTotalNumber(); p1++)
      {
        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell0DTotalNumber(); p2++)
        {
          Output::Assert(!geometryUtilities.PointsAreCoincident(convexMesh.Cell0DCoordinates(p1),
                                                                convexMesh.Cell0DCoordinates(p2)));
        }
      }
    }

    if (configuration.Cell1D_CheckDuplications)
    {
      for (unsigned int e1 = 0; e1 < convexMesh.Cell1DTotalNumber(); e1++)
      {
        Output::Assert(convexMesh.Cell1DByExtremes(convexMesh.Cell1DOrigin(e1),
                                                   convexMesh.Cell1DEnd(e1)) ==
                       e1);
        Output::Assert(convexMesh.Cell1DByExtremes(convexMesh.Cell1DEnd(e1),
                                                   convexMesh.Cell1DOrigin(e1)) ==
                       convexMesh.Cell1DTotalNumber());

        for (unsigned int e2 = e1 + 1; e2 < convexMesh.Cell1DTotalNumber(); e2++)
        {
          Output::Assert(!(convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DEnd(e2)));
          Output::Assert(!(convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DEnd(e2)));
        }
      }
    }

    if (configuration.Cell1D_CheckNeighbours)
    {
      for (unsigned int e = 0; e < convexMesh.Cell1DTotalNumber(); e++)
      {
        Output::Assert(convexMesh.Cell1DNumberNeighbourCell2D(e) > 0);

        for (unsigned int n = 0; n < convexMesh.Cell1DNumberNeighbourCell2D(e); n++)
        {
          if (!convexMesh.Cell1DHasNeighbourCell2D(e, n))
            continue;

          const unsigned int cell2DIndex = convexMesh.Cell1DNeighbourCell2D(e, n);
          const unsigned int cell2DNumEdges = convexMesh.Cell2DNumberEdges(cell2DIndex);

          // check edge orientation
          const unsigned int cell2DEdgeIndex = convexMesh.Cell2DFindEdge(cell2DIndex,
                                                                         e);
          const unsigned int edgeOrigin = convexMesh.Cell2DVertex(cell2DIndex,
                                                                  (cell2DEdgeIndex + 1) % cell2DNumEdges);
          const unsigned int edgeEnd = convexMesh.Cell2DVertex(cell2DIndex,
                                                               cell2DEdgeIndex);

          Output::Assert((convexMesh.Cell2DFindEdgeByExtremes(cell2DIndex,
                                                              edgeOrigin,
                                                              edgeEnd) ==
                          cell2DEdgeIndex &&
                          convexMesh.Cell2DFindEdgeByExtremes(cell2DIndex,
                                                              edgeEnd,
                                                              edgeOrigin) ==
                          cell2DNumEdges) ||
                         (convexMesh.Cell2DFindEdgeByExtremes(cell2DIndex,
                                                              edgeOrigin,
                                                              edgeEnd) ==
                          cell2DNumEdges &&
                          convexMesh.Cell2DFindEdgeByExtremes(cell2DIndex,
                                                              edgeEnd,
                                                              edgeOrigin) ==
                          cell2DEdgeIndex));
        }
      }
    }

    if (configuration.Cell1D_CheckMeasure)
    {
      for (unsigned int e = 0; e < convexMesh.Cell1DTotalNumber(); e++)
      {
        Output::Assert(geometryUtilities.IsValue1DPositive(
                         geometryUtilities.SegmentLength(convexMesh.Cell1DOriginCoordinates(e),
                                                         convexMesh.Cell1DEndCoordinates(e))));
      }
    }

    if (configuration.Cell2D_CheckEdges)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        const unsigned int cell2DNumEdges = convexMesh.Cell2DNumberEdges(p);
        for (unsigned int v = 0; v < cell2DNumEdges; v++)
        {
          const unsigned int eO = convexMesh.Cell2DVertex(p, v);
          const unsigned int eE = convexMesh.Cell2DVertex(p, (v + 1) % cell2DNumEdges);

          const unsigned int edgeFromVerticesOE = convexMesh.Cell2DFindEdgeByExtremes(p,
                                                                                      eO,
                                                                                      eE);
          const unsigned int edgeFromVerticesEO = convexMesh.Cell2DFindEdgeByExtremes(p,
                                                                                      eE,
                                                                                      eO);

          Output::Assert((edgeFromVerticesOE < cell2DNumEdges && edgeFromVerticesOE == v) ||
                         (edgeFromVerticesEO < cell2DNumEdges && edgeFromVerticesEO == v));
        }
      }
    }

    if (configuration.Cell2D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell2DTotalNumber(); p1++)
      {
        vector<unsigned int> cell2D1Vertices = convexMesh.Cell2DVertices(p1);
        sort(cell2D1Vertices.begin(), cell2D1Vertices.end());
        vector<unsigned int> cell2D1Edges = convexMesh.Cell2DEdges(p1);
        sort(cell2D1Edges.begin(), cell2D1Edges.end());

        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell2DTotalNumber(); p2++)
        {
          vector<unsigned int> cell2D2Vertices = convexMesh.Cell2DVertices(p2);
          sort(cell2D2Vertices.begin(), cell2D2Vertices.end());
          vector<unsigned int> cell2D2Edges = convexMesh.Cell2DEdges(p2);
          sort(cell2D2Edges.begin(), cell2D2Edges.end());

          Output::Assert(cell2D1Vertices.size() != cell2D2Vertices.size() || !equal(cell2D1Vertices.begin(),
                                                                                    cell2D1Vertices.end(),
                                                                                    cell2D2Vertices.begin()));
          Output::Assert(cell2D1Edges.size() != cell2D2Edges.size() || !equal(cell2D1Edges.begin(),
                                                                              cell2D1Edges.end(),
                                                                              cell2D2Edges.begin()));
        }
      }
    }

    if (configuration.Cell2D_CheckConvexity)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        const Eigen::MatrixXd cell2DVertices = convexMesh.Cell2DVerticesCoordinates(p);
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(cell2DVertices);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(cell2DVertices,
                                                                                              convexCell2DUnalignedVerticesFilter);
        const vector<unsigned int> convexHull = geometryUtilities.ConvexHull(convexCell2DUnalignedVertices);
        const Eigen::MatrixXd convexHullVertices = geometryUtilities.ExtractPoints(convexCell2DUnalignedVertices,
                                                                                   convexHull);

        Output::Assert(geometryUtilities.PolygonIsConvex(convexCell2DUnalignedVertices,
                                                         convexHullVertices));
      }
    }

    if (configuration.Cell2D_CheckMeasure)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        Output::Assert(geometryUtilities.IsValue2DPositive(
                         geometryUtilities.PolygonArea(convexMesh.Cell2DVerticesCoordinates(p))));
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::Mesh2DFromPolygon(const Eigen::MatrixXd& polygonVertices,
                                        const vector<unsigned int> vertexMarkers,
                                        const vector<unsigned int> edgeMarkers,
                                        IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(2);

    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);
    const unsigned int numPolygonVertices = polygonVertices.cols();
    Output::Assert(vertexMarkers.size() == numPolygonVertices);
    Output::Assert(edgeMarkers.size() == numPolygonVertices);

    // Create Cell0Ds
    const unsigned int& numCell0Ds = numPolygonVertices;
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      mesh.Cell0DSetState(v, true);
      mesh.Cell0DInsertCoordinates(v, polygonVertices.col(v));
      mesh.Cell0DSetMarker(v, vertexMarkers[v]);
    }

    // Create Cell1Ds
    unsigned int numCell1Ds = numPolygonVertices;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                e,
                                (e + 1) % numPolygonVertices);
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, edgeMarkers[e]);
    }

    // Create Cell2Ds
    const unsigned int& numCell2Ds = 1;
    mesh.Cell2DsInitialize(numCell2Ds);

    mesh.Cell2DInitializeVertices(0, numPolygonVertices);
    mesh.Cell2DInitializeEdges(0, numPolygonVertices);

    for (unsigned int v = 0; v < numPolygonVertices; v++)
      mesh.Cell2DInsertVertex(0, v, v);
    for (unsigned int e = 0; e < numPolygonVertices; e++)
      mesh.Cell2DInsertEdge(0, e, e);

    mesh.Cell2DSetState(0, true);
    mesh.Cell2DSetMarker(0, 0);

    // Create Cell1D neighbours
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      mesh.Cell1DInitializeNeighbourCell2Ds(e, 2);
      mesh.Cell1DInsertNeighbourCell2D(e, 1, 0);
    }
  }
  // ***************************************************************************
  vector<unsigned int> MeshUtilities::MeshCell2DRoots(const IMeshDAO& mesh) const
  {
    vector<unsigned int> rootCell2Ds(mesh.Cell2DTotalNumber());

    for (unsigned int cc = 0; cc < mesh.Cell2DTotalNumber(); cc++)
    {
      unsigned int rootCell = cc;
      while (mesh.Cell2DHasOriginalCell2D(rootCell))
        rootCell = mesh.Cell2DOriginalCell2D(rootCell);

      rootCell2Ds[cc] = rootCell;
    }

    return rootCell2Ds;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData2D MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& convexMesh) const
  {
    MeshGeometricData2D result;

    result.Cell2DsVertices.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsTriangulations.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsAreas.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsCentroids.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsDiameters.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeDirections.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeLengths.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeTangents.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeNormals.resize(convexMesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < convexMesh.Cell2DTotalNumber(); c++)
    {
      if (!convexMesh.Cell2DIsActive(c))
        continue;

      const unsigned int& domainCell2DIndex = c;

      // Extract original cell2D geometric information
      const Eigen::MatrixXd convexCell2DVertices = convexMesh.Cell2DVerticesCoordinates(domainCell2DIndex);

      // compute original cell2D triangulation
      const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
      const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
                                                                                            convexCell2DUnalignedVerticesFilter);

      const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
      vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
      for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
        convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

      const vector<Eigen::Matrix3d> convexCell2DTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(convexCell2DVertices,
                                                                                                                   convexCell2DTriangulation);

      const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();
      unsigned int cell2DTriangulationSize = numConvexCell2DTriangulation;

      // compute original cell2D area and centroids
      Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
      Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
      for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
      {
        convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
        convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
      }

      const double convexCell2DArea = convexCell2DTriangulationAreas.sum();
      const Eigen::Vector3d convexCell2DCentroid = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                                     convexCell2DTriangulationAreas,
                                                                                     convexCell2DArea);

      result.Cell2DsVertices[c] = convexCell2DVertices;

      // Compute cell2D triangulation from original cell2Ds
      result.Cell2DsTriangulations[c].resize(cell2DTriangulationSize);
      unsigned int triangulationCounter = 0;
      for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints.size(); cct++)
        result.Cell2DsTriangulations[c][triangulationCounter++] = convexCell2DTriangulationPoints[cct];

      result.Cell2DsAreas[c] = convexCell2DArea;
      result.Cell2DsCentroids[c] = convexCell2DCentroid;
      result.Cell2DsDiameters[c] = geometryUtilities.PolygonDiameter(result.Cell2DsVertices[c]);

      const unsigned int cell2DNumEdges = convexMesh.Cell2DNumberEdges(domainCell2DIndex);
      result.Cell2DsEdgeDirections[c].resize(cell2DNumEdges);
      for (unsigned int e = 0; e < cell2DNumEdges; e++)
      {
        const unsigned int origin = convexMesh.Cell2DVertex(domainCell2DIndex, e);
        const unsigned int end = convexMesh.Cell2DVertex(domainCell2DIndex,
                                                         (e + 1) % cell2DNumEdges);

        result.Cell2DsEdgeDirections[c][e] = convexMesh.Cell2DFindEdgeByExtremes(domainCell2DIndex,
                                                                                 origin,
                                                                                 end) == e;
      }

      result.Cell2DsEdgeLengths[c] = geometryUtilities.PolygonEdgeLengths(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeTangents[c] = geometryUtilities.PolygonEdgeTangents(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeNormals[c] = geometryUtilities.PolygonEdgeNormals(result.Cell2DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData2D MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& mesh,
                                                                            const std::vector<GeometryUtilities::PolygonTypes>& meshCell2DsPolygonType) const
  {
    MeshGeometricData2D result;

    result.Cell2DsVertices.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsTriangulations.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsAreas.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsCentroids.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsDiameters.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeDirections.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeLengths.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeTangents.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeNormals.resize(mesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      if (!mesh.Cell2DIsActive(c))
        continue;

      const unsigned int& domainCell2DIndex = c;

      // Extract original cell2D geometric information
      const Eigen::MatrixXd cell2DVertices = mesh.Cell2DVerticesCoordinates(domainCell2DIndex);

      // compute cell2D triangulation
      vector<Eigen::Matrix3d> cell2DTriangulationPoints;

      switch (meshCell2DsPolygonType[domainCell2DIndex]) {
        case GeometryUtilities::PolygonTypes::Triangle:
        case GeometryUtilities::PolygonTypes::Quadrilateral_Convex:
        case GeometryUtilities::PolygonTypes::Generic_Convex:
        {
          const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(cell2DVertices);
          const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(cell2DVertices,
                                                                                                convexCell2DUnalignedVerticesFilter);

          const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
          vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
          for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
            convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

          cell2DTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(cell2DVertices,
                                                                                   convexCell2DTriangulation);
        }
          break;
        case GeometryUtilities::PolygonTypes::Quadrilateral_Concave:
        case GeometryUtilities::PolygonTypes::Generic_Concave:
        {
          const vector<unsigned int> concaveCell2DTriangulation = geometryUtilities.PolygonTriangulationByEarClipping(cell2DVertices);
          cell2DTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(cell2DVertices,
                                                                                   concaveCell2DTriangulation);
        }
          break;
        default:
          throw runtime_error("Unsupported polygon type");
      }



      const unsigned int& numCell2DTriangulation = cell2DTriangulationPoints.size();
      unsigned int cell2DTriangulationSize = numCell2DTriangulation;

      // compute original cell2D area and centroids
      Eigen::VectorXd cell2DTriangulationAreas(numCell2DTriangulation);
      Eigen::MatrixXd cell2DTriangulationCentroids(3, numCell2DTriangulation);
      for (unsigned int cct = 0; cct < numCell2DTriangulation; cct++)
      {
        cell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(cell2DTriangulationPoints[cct]);
        cell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(cell2DTriangulationPoints[cct]);
      }

      const double cell2DArea = cell2DTriangulationAreas.sum();
      const Eigen::Vector3d cell2DCentroid = geometryUtilities.PolygonCentroid(cell2DTriangulationCentroids,
                                                                               cell2DTriangulationAreas,
                                                                               cell2DArea);

      result.Cell2DsVertices[c] = cell2DVertices;

      result.Cell2DsTriangulations[c].resize(cell2DTriangulationSize);
      unsigned int triangulationCounter = 0;
      for (unsigned int cct = 0; cct < cell2DTriangulationPoints.size(); cct++)
        result.Cell2DsTriangulations[c][triangulationCounter++] = cell2DTriangulationPoints[cct];

      result.Cell2DsAreas[c] = cell2DArea;
      result.Cell2DsCentroids[c] = cell2DCentroid;
      result.Cell2DsDiameters[c] = geometryUtilities.PolygonDiameter(result.Cell2DsVertices[c]);

      const unsigned int cell2DNumEdges = mesh.Cell2DNumberEdges(domainCell2DIndex);
      result.Cell2DsEdgeDirections[c].resize(cell2DNumEdges);
      for (unsigned int e = 0; e < cell2DNumEdges; e++)
      {
        const unsigned int origin = mesh.Cell2DVertex(domainCell2DIndex, e);
        const unsigned int end = mesh.Cell2DVertex(domainCell2DIndex,
                                                   (e + 1) % cell2DNumEdges);

        result.Cell2DsEdgeDirections[c][e] = mesh.Cell2DFindEdgeByExtremes(domainCell2DIndex,
                                                                           origin,
                                                                           end) == e;
      }

      result.Cell2DsEdgeLengths[c] = geometryUtilities.PolygonEdgeLengths(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeTangents[c] = geometryUtilities.PolygonEdgeTangents(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeNormals[c] = geometryUtilities.PolygonEdgeNormals(result.Cell2DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData2D MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& mesh,
                                                                            const IMeshDAO& convexMesh,
                                                                            const vector<vector<unsigned int>>& meshCell2DToConvexCell2DIndices) const
  {
    MeshGeometricData2D result;

    result.Cell2DsVertices.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsTriangulations.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsAreas.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsCentroids.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsDiameters.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeDirections.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeLengths.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeTangents.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeNormals.resize(mesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      if (!mesh.Cell2DIsActive(c))
        continue;

      const unsigned int& domainCell2DIndex = c;
      const vector<unsigned int>& domainConvexCell2DIndices = meshCell2DToConvexCell2DIndices[domainCell2DIndex];
      const unsigned int& numConvexCells = domainConvexCell2DIndices.size();

      // Get domain cell2D geometry information
      map<unsigned int, unsigned int> cell2DVerticesToPosition;
      for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(domainCell2DIndex); v++)
        cell2DVerticesToPosition.insert(pair<unsigned int, unsigned int>(mesh.Cell2DVertex(domainCell2DIndex,
                                                                                           v),
                                                                         v));

      // Extract original cell2D geometric information
      unsigned int cell2DTriangulationSize = 0;
      vector<vector<Eigen::Matrix3d>> convexCell2DTriangulationPoints(numConvexCells);
      Eigen::VectorXd convexCell2DAreas(numConvexCells);
      Eigen::MatrixXd convexCell2DCentroids(3, numConvexCells);

      for (unsigned int cc = 0; cc < numConvexCells; cc++)
      {
        const unsigned int& domainConvexCell2DIndex = domainConvexCell2DIndices[cc];
        const Eigen::MatrixXd convexCell2DVertices = convexMesh.Cell2DVerticesCoordinates(domainConvexCell2DIndex);

        // compute original cell2D triangulation
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
                                                                                              convexCell2DUnalignedVerticesFilter);

        const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
        vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
        for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
          convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

        convexCell2DTriangulationPoints[cc] = geometryUtilities.ExtractTriangulationPoints(convexCell2DVertices,
                                                                                           convexCell2DTriangulation);

        const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints[cc].size();
        cell2DTriangulationSize += numConvexCell2DTriangulation;

        // compute original cell2D area and centroids
        Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
        Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
        {
          convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cc][cct]);
          convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationPoints[cc][cct],
                                                                                          convexCell2DTriangulationAreas[cct]);
        }

        convexCell2DAreas[cc] = convexCell2DTriangulationAreas.sum();
        convexCell2DCentroids.col(cc) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                          convexCell2DTriangulationAreas,
                                                                          convexCell2DAreas[cc]);
      }

      result.Cell2DsVertices[c] = mesh.Cell2DVerticesCoordinates(domainCell2DIndex);

      // Compute cell2D triangulation from original cell2Ds
      result.Cell2DsTriangulations[c].resize(cell2DTriangulationSize);
      unsigned int triangulationCounter = 0;
      for (unsigned int cc = 0; cc < numConvexCells; cc++)
      {
        for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints[cc].size(); cct++)
          result.Cell2DsTriangulations[c][triangulationCounter++] = convexCell2DTriangulationPoints[cc][cct];
      }

      result.Cell2DsAreas[c] = convexCell2DAreas.sum();
      result.Cell2DsCentroids[c] = geometryUtilities.PolygonCentroid(convexCell2DCentroids,
                                                                     convexCell2DAreas,
                                                                     result.Cell2DsAreas[c]);
      result.Cell2DsDiameters[c] = geometryUtilities.PolygonDiameter(result.Cell2DsVertices[c]);

      const unsigned int cell2DNumEdges = mesh.Cell2DNumberEdges(domainCell2DIndex);
      result.Cell2DsEdgeDirections[c].resize(cell2DNumEdges);
      for (unsigned int e = 0; e < cell2DNumEdges; e++)
      {
        const unsigned int origin = mesh.Cell2DVertex(domainCell2DIndex, e);
        const unsigned int end = mesh.Cell2DVertex(domainCell2DIndex,
                                                   (e + 1) % cell2DNumEdges);

        result.Cell2DsEdgeDirections[c][e] = mesh.Cell2DFindEdgeByExtremes(domainCell2DIndex,
                                                                           origin,
                                                                           end) == e;
      }

      result.Cell2DsEdgeLengths[c] = geometryUtilities.PolygonEdgeLengths(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeTangents[c] = geometryUtilities.PolygonEdgeTangents(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeNormals[c] = geometryUtilities.PolygonEdgeNormals(result.Cell2DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::ComputeCell1DCell2DNeighbours(IMeshDAO& mesh) const
  {
    // Initialize cell1D neighbours
    for (unsigned int c1D = 0; c1D < mesh.Cell1DTotalNumber(); c1D++)
      mesh.Cell1DInitializeNeighbourCell2Ds(c1D, 2);

    // Compute Cell1D neighbours starting from cell2Ds
    for (unsigned int c2D = 0; c2D < mesh.Cell2DTotalNumber(); c2D++)
    {
      const unsigned int numCell2DEdges = mesh.Cell2DNumberEdges(c2D);
      for (unsigned int e = 0; e < numCell2DEdges; e++)
      {
        const unsigned int cell1D = mesh.Cell2DEdge(c2D, e);
        const unsigned int edgeOrigin =  mesh.Cell2DVertex(c2D, e);
        const unsigned int edgeEnd =  mesh.Cell2DVertex(c2D, (e + 1) % numCell2DEdges);

        if (mesh.Cell2DFindEdgeByExtremes(c2D,
                                          edgeOrigin,
                                          edgeEnd) == e) // left cell
        {
          mesh.Cell1DInsertNeighbourCell2D(cell1D,
                                           1,
                                           c2D);
        }
        else // right cell
        {
          mesh.Cell1DInsertNeighbourCell2D(cell1D,
                                           0,
                                           c2D);
        }
      }
    }
  }
  // ***************************************************************************
  std::vector<unsigned int> MeshUtilities::SplitCell2D(const unsigned int& cell2DIndex,
                                                       const std::vector<Eigen::MatrixXi> subCell2Ds,
                                                       IMeshDAO& mesh) const
  {
    const unsigned int numSubCells = subCell2Ds.size();
    unsigned int newCell2DsStartingIndex = mesh.Cell2DAppend(numSubCells);

    vector<unsigned int> newCell2DsIndex(numSubCells);

    mesh.Cell2DSetState(cell2DIndex, false);

    for (unsigned int c = 0; c < numSubCells; c++)
    {
      newCell2DsIndex[c] = newCell2DsStartingIndex + c;

      const unsigned int& newCell2DIndex = newCell2DsIndex[c];
      mesh.Cell2DAddVerticesAndEdges(newCell2DIndex,
                                     subCell2Ds[c]);

      mesh.Cell2DSetMarker(newCell2DIndex, mesh.Cell2DMarker(cell2DIndex));
      mesh.Cell2DSetState(newCell2DIndex, true);

      mesh.Cell2DInsertUpdatedCell2D(cell2DIndex, newCell2DIndex);

      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(newCell2DIndex); e++)
      {
        const unsigned int cell1DIndex = mesh.Cell2DEdge(newCell2DIndex, e);

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); n++)
        {
          if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
            continue;

          if (mesh.Cell1DNeighbourCell2D(cell1DIndex, n) == cell2DIndex)
            mesh.Cell1DInsertNeighbourCell2D(cell1DIndex,
                                             n,
                                             newCell2DIndex);
        }
      }
    }

    return newCell2DsIndex;
  }
  // ***************************************************************************
  std::vector<unsigned int> MeshUtilities::FindCell2DsCommonVertices(const std::vector<unsigned int>& cell2DsIndex,
                                                                     const IMeshDAO& mesh) const
  {
    const unsigned int numPolygons = cell2DsIndex.size();

    std::vector<unsigned int> firstPolygonVertices = mesh.Cell2DVertices(cell2DsIndex[0]);

    if (numPolygons == 1)
      return firstPolygonVertices;

    // find intersections with second triangle
    std::vector<unsigned int> secondPolygonVertices = mesh.Cell2DVertices(cell2DsIndex[1]);

    std::vector<unsigned int> intersections;

    std::sort(firstPolygonVertices.begin(),
              firstPolygonVertices.end());
    std::sort(secondPolygonVertices.begin(),
              secondPolygonVertices.end());

    std::set_intersection(firstPolygonVertices.begin(),
                          firstPolygonVertices.end(),
                          secondPolygonVertices.begin(),
                          secondPolygonVertices.end(),
                          back_inserter(intersections));

    if (numPolygons == 2)
      return intersections;

    // find intersections with third triangle
    std::vector<unsigned int> thirdPolygonVertices = mesh.Cell2DVertices(cell2DsIndex[2]);

    std::vector<unsigned int> singleIntersection;

    std::sort(thirdPolygonVertices.begin(),
              thirdPolygonVertices.end());

    std::set_intersection(intersections.begin(),
                          intersections.end(),
                          thirdPolygonVertices.begin(),
                          thirdPolygonVertices.end(),
                          back_inserter(singleIntersection));

    if (singleIntersection.size() == 0)
      return singleIntersection;

    // check other triangles
    for (unsigned int t = 3; t < cell2DsIndex.size(); t++)
    {
      std::vector<unsigned int> polygonVertices = mesh.Cell2DVertices(cell2DsIndex[t]);
      Output::Assert(find(polygonVertices.begin(),
                          polygonVertices.end(),
                          singleIntersection[0]) != polygonVertices.end());
    }

    return singleIntersection;
  }
  // ***************************************************************************
  std::vector<unsigned int> MeshUtilities::FindCell2DsCommonEdges(const std::vector<unsigned int>& cell2DsIndex,
                                                                  const IMeshDAO& mesh) const
  {
    const unsigned int numPolygons = cell2DsIndex.size();

    std::vector<unsigned int> firstPolygonEdges = mesh.Cell2DEdges(cell2DsIndex[0]);

    if (numPolygons == 1)
      return firstPolygonEdges;

    // find intersections with second polygon
    std::vector<unsigned int> secondPolygonEdges = mesh.Cell2DEdges(cell2DsIndex[1]);

    std::vector<unsigned int> intersections;

    std::sort(firstPolygonEdges.begin(),
              firstPolygonEdges.end());
    std::sort(secondPolygonEdges.begin(),
              secondPolygonEdges.end());

    std::set_intersection(firstPolygonEdges.begin(),
                          firstPolygonEdges.end(),
                          secondPolygonEdges.begin(),
                          secondPolygonEdges.end(),
                          back_inserter(intersections));

    if (numPolygons == 2)
      return intersections;

    return {};
  }
  // ***************************************************************************
  MeshUtilities::AgglomerateTrianglesResult MeshUtilities::AgglomerateTriangles(const std::vector<unsigned int>& trianglesIndexToAgglomerate,
                                                                                IMeshDAO& triangularMesh) const
  {
    const unsigned int numTriangles = trianglesIndexToAgglomerate.size();

    const std::vector<unsigned int> commonVertex = FindCell2DsCommonVertices(trianglesIndexToAgglomerate,
                                                                             triangularMesh);

    Output::Assert(commonVertex.size() == 1);

    std::list<unsigned int> commonEdges;
    std::list<unsigned int> agglomeratePolygonVertices;
    std::list<unsigned int> agglomeratePolygonEdges;

    // first triangle
    {
      const unsigned int commonVertexLocalIndex = triangularMesh.Cell2DFindVertex(trianglesIndexToAgglomerate.at(0),
                                                                                  commonVertex[0]);
      agglomeratePolygonVertices.push_back(commonVertex[0]);
      agglomeratePolygonEdges.push_back(triangularMesh.Cell2DEdge(trianglesIndexToAgglomerate.at(0),
                                                                  commonVertexLocalIndex));
    }

    // middle triangles
    for (unsigned int t = 0; t < trianglesIndexToAgglomerate.size() - 1; t++)
    {
      const std::vector<unsigned int> commonEdge = FindCell2DsCommonEdges({
                                                                            trianglesIndexToAgglomerate.at(t),
                                                                            trianglesIndexToAgglomerate.at(t + 1)
                                                                          },
                                                                          triangularMesh);

      Output::Assert(commonEdge.size() == 1);

      commonEdges.push_back(commonEdge[0]);
      const unsigned int commonEdgeLocalIndex = triangularMesh.Cell2DFindEdge(trianglesIndexToAgglomerate.at(t),
                                                                              commonEdge[0]);
      const unsigned int oppositeCommonEdgeVertexIndex = (commonEdgeLocalIndex + 2) % 3;

      agglomeratePolygonVertices.push_back(triangularMesh.Cell2DVertex(trianglesIndexToAgglomerate.at(t),
                                                                       oppositeCommonEdgeVertexIndex));
      agglomeratePolygonEdges.push_back(triangularMesh.Cell2DEdge(trianglesIndexToAgglomerate.at(t),
                                                                  oppositeCommonEdgeVertexIndex));
    }

    // last triangle
    {
      const unsigned int commonEdgeLocalIndex = triangularMesh.Cell2DFindEdge(trianglesIndexToAgglomerate.at(numTriangles - 1),
                                                                              commonEdges.back());
      const unsigned int oppositeCommonEdgeVertexIndex = (commonEdgeLocalIndex + 2) % 3;
      const unsigned int commonVertexLocalIndex = triangularMesh.Cell2DFindVertex(trianglesIndexToAgglomerate.at(numTriangles - 1),
                                                                                  commonVertex[0]);
      const unsigned int oppositeCommonVertexVertexIndex = (commonVertexLocalIndex + 1) % 3;

      agglomeratePolygonVertices.push_back(triangularMesh.Cell2DVertex(trianglesIndexToAgglomerate.at(numTriangles - 1),
                                                                       oppositeCommonVertexVertexIndex));
      agglomeratePolygonVertices.push_back(triangularMesh.Cell2DVertex(trianglesIndexToAgglomerate.at(numTriangles - 1),
                                                                       oppositeCommonEdgeVertexIndex));

      agglomeratePolygonEdges.push_back(triangularMesh.Cell2DEdge(trianglesIndexToAgglomerate.at(numTriangles - 1),
                                                                  oppositeCommonVertexVertexIndex));
      agglomeratePolygonEdges.push_back(triangularMesh.Cell2DEdge(trianglesIndexToAgglomerate.at(numTriangles - 1),
                                                                  oppositeCommonEdgeVertexIndex));
    }

    AgglomerateTrianglesResult result;
    result.VerticesIndex = std::vector<unsigned int>(agglomeratePolygonVertices.begin(),
                                                     agglomeratePolygonVertices.end());

    result.EdgesIndex = std::vector<unsigned int>(agglomeratePolygonEdges.begin(),
                                                  agglomeratePolygonEdges.end());
    result.RemovedEdges = std::vector<unsigned int>(commonEdges.begin(),
                                                    commonEdges.end());

    return result;
  }
  // ***************************************************************************
  MeshUtilities::AgglomerateMeshFromTriangularMeshResult MeshUtilities::AgglomerateMeshFromTriangularMesh(const std::vector<std::vector<unsigned int> >& trianglesIndicesToAgglomerate,
                                                                                                          IMeshDAO& triangularMesh) const
  {
    AgglomerateMeshFromTriangularMeshResult result;

    list<unsigned int> removedCell1Ds;
    list<unsigned int> removedCell2Ds;

    result.ConcaveCell2Ds.resize(trianglesIndicesToAgglomerate.size());

    for (unsigned int ac = 0; ac < trianglesIndicesToAgglomerate.size(); ac++)
    {
      const std::vector<unsigned int>& triangleToAgglomerate = trianglesIndicesToAgglomerate[ac];

      const AgglomerateTrianglesResult agglomeratedPolygon = AgglomerateTriangles(triangleToAgglomerate,
                                                                                  triangularMesh);

      // create new cell2D
      AgglomerateMeshFromTriangularMeshResult::ConcaveCell2D& concaveCell2D = result.ConcaveCell2Ds[ac];

      concaveCell2D.Cell2DIndex = triangularMesh.Cell2DAppend(1);
      concaveCell2D.ConvexCell2DsIndex = triangleToAgglomerate;

      triangularMesh.Cell2DAddVertices(concaveCell2D.Cell2DIndex,
                                       agglomeratedPolygon.VerticesIndex);
      triangularMesh.Cell2DAddEdges(concaveCell2D.Cell2DIndex,
                                    agglomeratedPolygon.EdgesIndex);
      triangularMesh.Cell2DSetMarker(concaveCell2D.Cell2DIndex,
                                     0);
      triangularMesh.Cell2DSetState(concaveCell2D.Cell2DIndex,
                                    true);

      // deactivate other cells
      for (const unsigned int cell2DIndex : triangleToAgglomerate)
      {
        Gedim::Output::Assert(triangularMesh.Cell2DIsActive(cell2DIndex));

        triangularMesh.Cell2DSetState(cell2DIndex,
                                      false);
        triangularMesh.Cell2DInsertUpdatedCell2D(cell2DIndex,
                                                 concaveCell2D.Cell2DIndex);

        removedCell2Ds.push_back(cell2DIndex);
      }

      for (const unsigned int cell1DIndex : agglomeratedPolygon.RemovedEdges)
      {
        Gedim::Output::Assert(triangularMesh.Cell1DIsActive(cell1DIndex));

        triangularMesh.Cell1DSetState(cell1DIndex,
                                      false);

        removedCell1Ds.push_back(cell1DIndex);
      }
    }

    removedCell1Ds.sort();
    removedCell2Ds.sort();

    result.RemovedCell1Ds = std::vector<unsigned int>(removedCell1Ds.begin(),
                                                      removedCell1Ds.end());
    result.RemovedCell2Ds = std::vector<unsigned int>(removedCell2Ds.begin(),
                                                      removedCell2Ds.end());

    return result;
  }
  // ***************************************************************************
  MeshUtilities::AgglomerationInformation MeshUtilities::ImportAgglomerationInformationFromCsv(const GeometryUtilities geometryUtilities,
                                                                                               const IMeshDAO& originalMesh,
                                                                                               const IMeshDAO& agglomeratedMesh,
                                                                                               const std::string& fileName,
                                                                                               const char& separator) const
  {
    AgglomerationInformation result;

    // initialize
    result.OriginalCell0DToAgglomeratedCell0Ds.resize(originalMesh.Cell0DTotalNumber(),
                                                      originalMesh.Cell0DTotalNumber());
    result.OriginalCell1DToAgglomeratedCell1Ds.resize(originalMesh.Cell1DTotalNumber(),
                                                      originalMesh.Cell1DTotalNumber());
    result.OriginalCell2DToAgglomeratedCell2Ds.resize(originalMesh.Cell2DTotalNumber(),
                                                      originalMesh.Cell2DTotalNumber());
    result.AgglomeratedCell0DToOriginalCell0Ds.resize(agglomeratedMesh.Cell0DTotalNumber(),
                                                      agglomeratedMesh.Cell0DTotalNumber());
    result.AgglomeratedCell1DToOriginalCell1Ds.resize(agglomeratedMesh.Cell1DTotalNumber());
    result.AgglomeratedCell2DToOriginalCell2Ds.resize(agglomeratedMesh.Cell2DTotalNumber());

    // read agglomeration information file
    Gedim::FileReader file(fileName);
    vector<string> lines;
    file.Open();
    file.GetAllLines(lines);
    file.Close();

    Gedim::Output::Assert(lines.size() > 1);

    unsigned int lineCounter = 0;

    // read cell2D agglomeration
    istringstream converterCell2Ds(lines[lineCounter++]);
    unsigned int numAgglomeratedCell2Ds;
    converterCell2Ds>> numAgglomeratedCell2Ds;
    Gedim::Output::Assert(numAgglomeratedCell2Ds == agglomeratedMesh.Cell2DTotalNumber());

    lineCounter++; // ingnore header
    for (unsigned int ac = 0; ac < numAgglomeratedCell2Ds; ac++)
    {
      char temp;
      istringstream converterAgglomeratedCell(lines[lineCounter++]);
      unsigned int agglomeratedCell2D, numOriginalCell2Ds;
      converterAgglomeratedCell >>agglomeratedCell2D;
      if (separator != ' ')
        converterAgglomeratedCell >> temp;
      converterAgglomeratedCell >>numOriginalCell2Ds;
      if (separator != ' ')
        converterAgglomeratedCell >> temp;

      Gedim::Output::Assert(agglomeratedCell2D < agglomeratedMesh.Cell2DTotalNumber() &&
                            numOriginalCell2Ds > 0);

      result.AgglomeratedCell2DToOriginalCell2Ds[agglomeratedCell2D].resize(numOriginalCell2Ds);
      for (unsigned int oc = 0; oc < numOriginalCell2Ds; oc++)
      {
        unsigned int originalCell2D;
        converterAgglomeratedCell >>originalCell2D;
        if (separator != ' ')
          converterAgglomeratedCell >> temp;

        Gedim::Output::Assert(originalCell2D < originalMesh.Cell2DTotalNumber());

        result.AgglomeratedCell2DToOriginalCell2Ds[agglomeratedCell2D][oc] = originalCell2D;
        result.OriginalCell2DToAgglomeratedCell2Ds[originalCell2D] = agglomeratedCell2D;
      }
    }

    lineCounter++; // ingnore empty line

    // read cell0D agglomeration
    istringstream converterCell0Ds(lines[lineCounter++]);
    unsigned int numAgglomeratedCell0Ds;
    converterCell0Ds>> numAgglomeratedCell0Ds;
    Gedim::Output::Assert(numAgglomeratedCell0Ds == agglomeratedMesh.Cell0DTotalNumber());

    lineCounter++; // ingnore header
    for (unsigned int ac = 0; ac < numAgglomeratedCell0Ds; ac++)
    {
      char temp;
      istringstream converterAgglomeratedCell(lines[lineCounter++]);
      unsigned int agglomeratedCell0D, originalCell0D;
      converterAgglomeratedCell >>agglomeratedCell0D;
      if (separator != ' ')
        converterAgglomeratedCell >> temp;
      converterAgglomeratedCell >>originalCell0D;
      if (separator != ' ')
        converterAgglomeratedCell >> temp;

      Gedim::Output::Assert(agglomeratedCell0D < agglomeratedMesh.Cell0DTotalNumber());

      if (originalCell0D < originalMesh.Cell0DTotalNumber())
        result.OriginalCell0DToAgglomeratedCell0Ds[originalCell0D] = agglomeratedCell0D;

      result.AgglomeratedCell0DToOriginalCell0Ds[agglomeratedCell0D] = originalCell0D;
    }

    // create cell0D-cell1D relation
    struct Edge
    {
        unsigned int Cell1DIndex;
        unsigned int Cell0DEnd;
    };
    std::vector<std::list<Edge>> cell0DsCell1Ds(originalMesh.Cell0DTotalNumber());
    vector<list<unsigned int>> originalCell0DsCell1Ds(originalMesh.Cell0DTotalNumber());

    for (unsigned int c1D = 0; c1D < originalMesh.Cell1DTotalNumber(); c1D++)
    {
      const unsigned int originalCell1DOrigin = originalMesh.Cell1DOrigin(c1D);
      const unsigned int originalCell1DEnd = originalMesh.Cell1DEnd(c1D);

      originalCell0DsCell1Ds[originalCell1DOrigin].push_back(c1D);
      originalCell0DsCell1Ds[originalCell1DEnd].push_back(c1D);
      cell0DsCell1Ds[originalMesh.Cell1DOrigin(c1D)].push_back({c1D, originalMesh.Cell1DEnd(c1D)});
    }

    // compute agglomerate cell1D and original cell1D relation
    for (unsigned int ac1D = 0; ac1D < agglomeratedMesh.Cell1DTotalNumber(); ac1D++)
    {
      const unsigned int cell1DAgglomeratedOrigin = agglomeratedMesh.Cell1DOrigin(ac1D);
      const unsigned int cell1DAgglomeratedEnd = agglomeratedMesh.Cell1DEnd(ac1D);

      const unsigned int cell1DOriginalOrigin = result.AgglomeratedCell0DToOriginalCell0Ds[cell1DAgglomeratedOrigin];
      const unsigned int cell1DOriginalEnd = result.AgglomeratedCell0DToOriginalCell0Ds[cell1DAgglomeratedEnd];

      if (cell1DOriginalOrigin >=  originalMesh.Cell0DTotalNumber() ||
          cell1DOriginalEnd >= originalMesh.Cell0DTotalNumber())
        continue;

      const std::list<Edge>& originCell1Ds = cell0DsCell1Ds.at(cell1DOriginalOrigin);
      std::list<Edge>::const_iterator findOriginEdge = std::find_if(originCell1Ds.begin(),
                                                                    originCell1Ds.end(),
                                                                    [&] (const Edge& edge)
      { return edge.Cell0DEnd == cell1DOriginalEnd; });

      const std::list<Edge>& endCell1Ds = cell0DsCell1Ds.at(cell1DOriginalEnd);
      std::list<Edge>::const_iterator findEndEdge = std::find_if(endCell1Ds.begin(),
                                                                 endCell1Ds.end(),
                                                                 [&] (const Edge& edge)
      { return edge.Cell0DEnd == cell1DOriginalOrigin; });

      if (findOriginEdge != originCell1Ds.end())
      {
        const unsigned int originalCell1D = findOriginEdge->Cell1DIndex;
        result.OriginalCell1DToAgglomeratedCell1Ds[originalCell1D] = ac1D;
        result.AgglomeratedCell1DToOriginalCell1Ds[ac1D].resize(1, originalCell1D);
      }
      else if (findEndEdge != endCell1Ds.end())
      {
        const unsigned int originalCell1D = findEndEdge->Cell1DIndex;
        result.OriginalCell1DToAgglomeratedCell1Ds[originalCell1D] = ac1D;
        result.AgglomeratedCell1DToOriginalCell1Ds[ac1D].resize(1, originalCell1D);
      }
      else
      {
        // search original cell1Ds aligned and inside the agglomerated cell1D
        const Eigen::Vector3d& endCoordinates = originalMesh.Cell0DCoordinates(cell1DOriginalEnd);

        list<unsigned int> originalCell1Ds;
        queue<unsigned int> origins;
        origins.push(cell1DOriginalOrigin);

        while (!origins.empty())
        {
          const unsigned int cell0D = origins.front();
          origins.pop();

          if (originalCell0DsCell1Ds[cell0D].size() == 0)
            continue;

          const Vector3d cell0DCoordinates = originalMesh.Cell0DCoordinates(cell0D);
          const list<unsigned int>& cell1Ds = originalCell0DsCell1Ds[cell0D];

          unsigned int cell1DAligned = originalMesh.Cell1DTotalNumber();
          for (unsigned int cell1D : cell1Ds)
          {
            if (result.OriginalCell1DToAgglomeratedCell1Ds[cell1D] !=
                originalMesh.Cell1DTotalNumber())
              continue;

            const unsigned int cell1DOtherCell0DIndex = (originalMesh.Cell1DOrigin(cell1D) == cell0D) ?
                                                          originalMesh.Cell1DEnd(cell1D) :
                                                          originalMesh.Cell1DOrigin(cell1D);

            if (cell1DOtherCell0DIndex == cell1DOriginalEnd)
            {
              cell1DAligned = cell1D;
              break;
            }

            const Vector3d cell1DOtherCell0DCoordinates = originalMesh.Cell0DCoordinates(cell1DOtherCell0DIndex);

            if (!geometryUtilities.PointIsAligned(cell0DCoordinates,
                                                  endCoordinates,
                                                  cell1DOtherCell0DCoordinates))
              continue;

            if (geometryUtilities.PointSegmentPosition(cell1DOtherCell0DCoordinates,
                                                       cell0DCoordinates,
                                                       endCoordinates) !=
                Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment)
              continue;

            cell1DAligned = cell1D;
            if (cell1DOtherCell0DIndex != cell1DOriginalEnd)
              origins.push(cell1DOtherCell0DIndex);

            break;
          }

          Gedim::Output::Assert(cell1DAligned != originalMesh.Cell1DTotalNumber());

          result.OriginalCell1DToAgglomeratedCell1Ds[cell1DAligned] = ac1D;
          originalCell1Ds.push_back(cell1DAligned);
        }

        result.AgglomeratedCell1DToOriginalCell1Ds[ac1D] = vector<unsigned int>(originalCell1Ds.begin(),
                                                                                originalCell1Ds.end());
      }
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::ExportConcaveMesh2DToCsv(const IMeshDAO& mesh,
                                               const std::vector<std::vector<unsigned int> >& convexCell2DsIndex,
                                               const char& separator,
                                               const std::string& exportFolderPath) const
  {
    // export mesh
    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
    Gedim::MeshFromCsvUtilities::Configuration exportConfiguration;
    exportConfiguration.Separator = separator;
    exportConfiguration.Folder = exportFolderPath;
    Gedim::MeshDAOExporterToCsv exporter(meshFromCsvUtilities);
    exporter.Export(exportConfiguration,
                    mesh);

    // export concave to convex mesh file
    {
      ofstream mapFile;

      mapFile.open(exportFolderPath +
                   "/hierarchy_map.txt");
      mapFile.precision(16);

      if (mapFile.fail())
        throw runtime_error("Error on hierarchy_map file");

      mapFile<< mesh.Cell2DTotalNumber()<< endl;
      mapFile<< "# newCellId, sizeOldCellIdsContainer, oldCellIds"<< endl;
      for (unsigned int v = 0; v < mesh.Cell2DTotalNumber(); v++)
      {
        mapFile<< scientific<< v<< separator;
        mapFile<< scientific<< convexCell2DsIndex[v].size();
        for (unsigned int cc = 0; cc < convexCell2DsIndex[v].size(); cc++)
          mapFile<< scientific<< separator<< convexCell2DsIndex[v][cc];
        mapFile<< endl;
      }

      mapFile<< endl;
      mapFile<< mesh.Cell0DTotalNumber()<< endl;
      mapFile<< "# newVertId, oldVertId"<< endl;
      for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
      {
        mapFile<< scientific<< v<< separator;
        mapFile<< scientific<< v<< endl;
      }

      mapFile.close();
    }
  }
  // ***************************************************************************
  void MeshUtilities::ChangePolygonMeshMarkers(const Eigen::MatrixXd& polygonVertices,
                                               const vector<unsigned int>& cell0DMarkers,
                                               const vector<unsigned int>& cell1DMarkers,
                                               IMeshDAO& mesh) const
  {
    Output::Assert(mesh.Dimension() == 2);

    const unsigned int numVertices = polygonVertices.cols();

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      if (mesh.Cell0DMarker(v) == 0)
        continue;

      const unsigned int newMarker = (mesh.Cell0DMarker(v) <= numVertices) ?
                                       cell0DMarkers.at(mesh.Cell0DMarker(v) - 1) :
                                       cell1DMarkers.at(mesh.Cell0DMarker(v) - numVertices - 1);

      mesh.Cell0DSetMarker(v,
                           newMarker);
    }

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
    {
      if (mesh.Cell1DMarker(e) == 0)
        continue;

      const unsigned int newMarker = cell1DMarkers.at(mesh.Cell1DMarker(e) - numVertices - 1);

      mesh.Cell1DSetMarker(e,
                           newMarker);
    }
  }
  // ***************************************************************************
}
