#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::Polyhedron GeometryUtilities::CreateTetrahedronWithOrigin(const Eigen::Vector3d& origin,
                                                                               const Eigen::Vector3d& lengthVector,
                                                                               const Eigen::Vector3d& heightVector,
                                                                               const Eigen::Vector3d& widthVector) const
  {
    GeometryUtilities::Polyhedron tetrahedron;

    // create vertices
    tetrahedron.Vertices.resize(3, 4);
    tetrahedron.Vertices.col(0) << origin;
    tetrahedron.Vertices.col(1) << origin + lengthVector;
    tetrahedron.Vertices.col(2) << origin + widthVector;
    tetrahedron.Vertices.col(3) << origin + heightVector;

    // create edges
    tetrahedron.Edges.resize(2, 6);
    tetrahedron.Edges.col(0) << 0, 1;
    tetrahedron.Edges.col(1) << 0, 2;
    tetrahedron.Edges.col(2) << 1, 2;
    tetrahedron.Edges.col(3) << 0, 3;
    tetrahedron.Edges.col(4) << 1, 3;
    tetrahedron.Edges.col(5) << 2, 3;

    // create faces
    tetrahedron.Faces.resize(4, MatrixXi::Zero(2, 3));
    tetrahedron.Faces[0].row(0)<< 0, 1, 2;
    tetrahedron.Faces[1].row(0)<< 0, 1, 3;
    tetrahedron.Faces[2].row(0)<< 0, 2, 3;
    tetrahedron.Faces[3].row(0)<< 1, 2, 3;

    tetrahedron.Faces[0].row(1)<< 0, 2, 1;
    tetrahedron.Faces[1].row(1)<< 0, 4, 3;
    tetrahedron.Faces[2].row(1)<< 1, 5, 3;
    tetrahedron.Faces[3].row(1)<< 2, 5, 4;

    return tetrahedron;
  }
  // ***************************************************************************
  GeometryUtilities::Polyhedron GeometryUtilities::CreateTetrahedronWithVertices(const Eigen::Vector3d& v1,
                                                                                 const Eigen::Vector3d& v2,
                                                                                 const Eigen::Vector3d& v3,
                                                                                 const Eigen::Vector3d& v4) const
  {
    GeometryUtilities::Polyhedron tetrahedron;

    // create vertices
    tetrahedron.Vertices.resize(3, 4);
    tetrahedron.Vertices.col(0) << v1;
    tetrahedron.Vertices.col(1) << v2;
    tetrahedron.Vertices.col(2) << v3;
    tetrahedron.Vertices.col(3) << v4;

    // create edges
    tetrahedron.Edges.resize(2, 6);
    tetrahedron.Edges.col(0)<< 0, 1;
    tetrahedron.Edges.col(1)<< 0, 2;
    tetrahedron.Edges.col(2)<< 1, 2;
    tetrahedron.Edges.col(3)<< 0, 3;
    tetrahedron.Edges.col(4)<< 1, 3;
    tetrahedron.Edges.col(5)<< 2, 3;

    // create faces
    tetrahedron.Faces.resize(4, MatrixXi::Zero(2, 3));
    tetrahedron.Faces[0].row(0)<< 0, 1, 2;
    tetrahedron.Faces[1].row(0)<< 0, 1, 3;
    tetrahedron.Faces[2].row(0)<< 0, 2, 3;
    tetrahedron.Faces[3].row(0)<< 1, 2, 3;

    tetrahedron.Faces[0].row(1)<< 0, 2, 1;
    tetrahedron.Faces[1].row(1)<< 0, 4, 3;
    tetrahedron.Faces[2].row(1)<< 1, 5, 3;
    tetrahedron.Faces[3].row(1)<< 2, 5, 4;

    return tetrahedron;
  }
  // ***************************************************************************
  GeometryUtilities::Polyhedron GeometryUtilities::CreateParallelepipedWithOrigin(const Eigen::Vector3d& origin,
                                                                                  const Eigen::Vector3d& lengthVector,
                                                                                  const Eigen::Vector3d& heightVector,
                                                                                  const Eigen::Vector3d& widthVector) const
  {
    Gedim::GeometryUtilities::Polyhedron parallelpiped;

    // create vertices
    parallelpiped.Vertices.setZero(3, 8);
    parallelpiped.Vertices.col(0)<< origin;
    parallelpiped.Vertices.col(1)<< origin + lengthVector;
    parallelpiped.Vertices.col(2)<< origin + lengthVector + widthVector;
    parallelpiped.Vertices.col(3)<< origin + widthVector;
    parallelpiped.Vertices.col(4)<< origin + heightVector;
    parallelpiped.Vertices.col(5)<< origin + heightVector + lengthVector;
    parallelpiped.Vertices.col(6)<< origin + heightVector + lengthVector + widthVector;
    parallelpiped.Vertices.col(7)<< origin + heightVector + widthVector;

    // create edges
    parallelpiped.Edges.setZero(2, 12);
    parallelpiped.Edges.col(0)<< 0, 1;
    parallelpiped.Edges.col(1)<< 1, 2;
    parallelpiped.Edges.col(2)<< 2, 3;
    parallelpiped.Edges.col(3)<< 3, 0;
    parallelpiped.Edges.col(4)<< 4, 5;
    parallelpiped.Edges.col(5)<< 5, 6;
    parallelpiped.Edges.col(6)<< 6, 7;
    parallelpiped.Edges.col(7)<< 7, 4;
    parallelpiped.Edges.col(8)<< 0, 4;
    parallelpiped.Edges.col(9)<< 1, 5;
    parallelpiped.Edges.col(10)<< 2, 6;
    parallelpiped.Edges.col(11)<< 3, 7;

    // create faces
    parallelpiped.Faces.resize(6, MatrixXi::Zero(2, 4));
    parallelpiped.Faces[0].row(0)<< 0, 1, 2, 3;
    parallelpiped.Faces[0].row(1)<< 0, 1, 2, 3;

    parallelpiped.Faces[1].row(0)<< 4, 5, 6, 7;
    parallelpiped.Faces[1].row(1)<< 4, 5, 6, 7;

    parallelpiped.Faces[2].row(0)<< 0, 3, 7, 4;
    parallelpiped.Faces[2].row(1)<< 3, 11, 7, 8;

    parallelpiped.Faces[3].row(0)<< 1, 2, 6, 5;
    parallelpiped.Faces[3].row(1)<< 1, 10, 5, 9;

    parallelpiped.Faces[4].row(0)<< 0, 1, 5, 4;
    parallelpiped.Faces[4].row(1)<< 0, 9, 4, 8;

    parallelpiped.Faces[5].row(0)<< 3, 2, 6, 7;
    parallelpiped.Faces[5].row(1)<< 2, 10, 6, 11;

    return parallelpiped;
  }

  MatrixXd GeometryUtilities::PolyhedronEdgeTangents(const Eigen::MatrixXd& polyhedronVertices,
                                                     const Eigen::MatrixXi& polyhedronEdges) const
  {
    MatrixXd edgeTangents(3, polyhedronEdges.cols());

    for (unsigned int e = 0; e < polyhedronEdges.cols(); e++)
    {
      edgeTangents.col(e) = SegmentTangent(polyhedronVertices.col(polyhedronEdges(0, e)),
                                           polyhedronVertices.col(polyhedronEdges(1, e)));
    }

    return edgeTangents;
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::PolyhedronFaceVertices(const Eigen::MatrixXd& polyhedronVertices,
                                                             const vector<Eigen::MatrixXi>& polyhedronFaces) const
  {
    vector<Eigen::MatrixXd> facesVertices(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const Eigen::VectorXi& faceVerticesIndices = polyhedronFaces[f].row(0);
      Eigen::MatrixXd& faceVertices = facesVertices[f];
      faceVertices.setZero(3, faceVerticesIndices.size());
      for (unsigned int v = 0; v < faceVerticesIndices.size(); v++)
        faceVertices.col(v)<< polyhedronVertices.col(faceVerticesIndices[v]);
    }

    return facesVertices;
  }
  // ***************************************************************************
  vector<vector<bool>> GeometryUtilities::PolyhedronFaceEdgeDirections(const Eigen::MatrixXd& polyhedronVertices,
                                                                       const Eigen::MatrixXi& polyhedronEdges,
                                                                       const vector<Eigen::MatrixXi>& polyhedronFaces) const
  {
    vector<vector<bool>> facesEdgeDirections(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const Eigen::VectorXi& faceEdgesIndices = polyhedronFaces[f].row(1);
      facesEdgeDirections[f].resize(faceEdgesIndices.size(), true);
      for (unsigned int e = 0; e < faceEdgesIndices.size(); e++)
        facesEdgeDirections[f][e] = (polyhedronEdges(0, faceEdgesIndices[e]) == polyhedronFaces[f](0, e));
    }

    return facesEdgeDirections;
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::PolyhedronFaceEdgeTangents(const Eigen::MatrixXd& polyhedronVertices,
                                                                 const Eigen::MatrixXi& polyhedronEdges,
                                                                 const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                 const vector<vector<bool>>& polyhedronFaceEdgeDirections,
                                                                 const Eigen::MatrixXd& polyhedronEdgeTangents) const
  {
    vector<Eigen::MatrixXd> facesEdgeTangents(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const Eigen::VectorXi& faceEdgesIndices = polyhedronFaces[f].row(1);
      Eigen::MatrixXd& faceEdgeTangents = facesEdgeTangents[f];
      faceEdgeTangents.setZero(3, faceEdgesIndices.size());
      for (unsigned int e = 0; e < faceEdgesIndices.size(); e++)
        faceEdgeTangents.col(e)<< (polyhedronFaceEdgeDirections[f][e] ? 1.0 : -1.0) *
                                  polyhedronEdgeTangents.col(faceEdgesIndices[e]);
    }

    return facesEdgeTangents;
  }
  // ***************************************************************************
  vector<Matrix3d> GeometryUtilities::PolyhedronFaceRotationMatrices(const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                     const vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                                     const vector<Eigen::Vector3d>& polyhedronFaceTranslations) const
  {
    vector<Matrix3d> rotationMatrices;
    rotationMatrices.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
    {
      rotationMatrices.push_back(PolygonRotationMatrix(polyhedronFaceVertices[f],
                                                       polyhedronFaceNormals[f],
                                                       polyhedronFaceTranslations[f]));
    }

    return rotationMatrices;
  }
  // ***************************************************************************
  vector<Vector3d> GeometryUtilities::PolyhedronFaceTranslations(const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<Vector3d> translations;
    translations.reserve(polyhedronFaceVertices.size());
    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
      translations.push_back(PolygonTranslation(polyhedronFaceVertices[f]));

    return translations;
  }
  // ***************************************************************************
  vector<Vector3d> GeometryUtilities::PolyhedronFaceNormals(const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<Vector3d> faceNormals;
    faceNormals.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
      faceNormals.push_back(PolygonNormal(polyhedronFaceVertices[f]));

    return faceNormals;
  }
  // ***************************************************************************
  vector<Vector3d> GeometryUtilities::PolyhedronFaceBarycenter(const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<Vector3d> faceBarycenters;
    faceBarycenters.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
      faceBarycenters.push_back(PolygonBarycenter(polyhedronFaceVertices[f]));

    return faceBarycenters;
  }
  // ***************************************************************************
  vector<bool> GeometryUtilities::PolyhedronFaceNormalDirections(const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                 const Eigen::Vector3d& pointInsidePolyhedron,
                                                                 const vector<Vector3d>& polyhedronFaceNormals) const
  {
    vector<bool> faceDirections;
    faceDirections.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
    {
      const Eigen::Vector3d& normal = polyhedronFaceNormals[f];
      const PointPlanePositionTypes pointFacePosition = PointPlanePosition(pointInsidePolyhedron,
                                                                           normal,
                                                                           polyhedronFaceVertices[f].col(0));
      Output::Assert(pointFacePosition == PointPlanePositionTypes::Negative ||
                     pointFacePosition == PointPlanePositionTypes::Positive);

      faceDirections.push_back(!(pointFacePosition == PointPlanePositionTypes::Positive));
    }

    return faceDirections;
  }
  // ***************************************************************************
  vector<vector<unsigned int> > GeometryUtilities::PolyhedronFaceTriangulations(const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                const vector<vector<unsigned int>>& localFaceTriangulations) const
  {
    vector<vector<unsigned int>> facesTriangulation(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      vector<unsigned int>& faceTriangulation = facesTriangulation[f];

      vector<unsigned int> faceLocalTriangulation = localFaceTriangulations[f];
      faceTriangulation.resize(faceLocalTriangulation.size());

      for (unsigned int v = 0; v < faceLocalTriangulation.size(); v++)
        faceTriangulation[v] = polyhedronFaces[f](0, faceLocalTriangulation[v]);
    }

    return facesTriangulation;
  }
  // ***************************************************************************
  vector<vector<unsigned int>> GeometryUtilities::PolyhedronFaceTriangulationsByFirstVertex(const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                            const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<vector<unsigned int>> localFacesTriangulations(polyhedronFaces.size());
    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
      localFacesTriangulations[f] = PolygonTriangulationByFirstVertex(polyhedronFaceVertices[f]);

    return PolyhedronFaceTriangulations(polyhedronFaces,
                                        localFacesTriangulations);
  }
  // ***************************************************************************
  vector<vector<unsigned int>> GeometryUtilities::PolyhedronFaceTriangulationsByInternalPoint(const Eigen::MatrixXd& polyhedronVertices,
                                                                                              const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                              const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                              const vector<Eigen::Vector3d>& polyhedronFaceInternalPoints) const
  {
    vector<vector<unsigned int>> facesTriangulation(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      vector<unsigned int>& faceTriangulation = facesTriangulation[f];

      vector<unsigned int> faceLocalTriangulation = PolygonTriangulationByInternalPoint(polyhedronFaceVertices[f],
                                                                                        polyhedronFaceInternalPoints[f]);
      faceTriangulation.resize(faceLocalTriangulation.size());

      for (unsigned int v = 0; v < faceLocalTriangulation.size(); v++)
        faceTriangulation[v] = (v % 3 == 0) ? polyhedronVertices.cols() + f :
                                              polyhedronFaces[f](0, faceLocalTriangulation[v]);
    }

    return facesTriangulation;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolyhedronTetrahedronsByFaceTriangulations(const Eigen::MatrixXd& polyhedronVertices,
                                                                                     const vector<vector<unsigned int>>& faceTriangulations,
                                                                                     const Eigen::Vector3d& polyhedronInternalPoint)
  {
    list<unsigned int> tetrahedronList;

    const unsigned int numPolyhedronVertices = polyhedronVertices.cols();
    for (unsigned int f = 0; f < faceTriangulations.size(); f++)
    {
      const unsigned int numFaceTriangulation = faceTriangulations[f].size() / 3;
      for (unsigned int ft = 0; ft < numFaceTriangulation; ft++)
      {
        tetrahedronList.push_back(faceTriangulations[f][3 * ft]);
        tetrahedronList.push_back(faceTriangulations[f][3 * ft + 1]);
        tetrahedronList.push_back(faceTriangulations[f][3 * ft + 2]);
        tetrahedronList.push_back(numPolyhedronVertices);
      }
    }

    Output::Assert(tetrahedronList.size() % 4 == 0);

    return vector<unsigned int>(tetrahedronList.begin(), tetrahedronList.end());
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolyhedronTetrahedronsByFaceTriangulations(const Eigen::MatrixXd& polyhedronVertices,
                                                                                     const vector<vector<unsigned int> >& faceTriangulations,
                                                                                     const vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                                                     const Eigen::Vector3d& polyhedronInternalPoint)
  {
    list<unsigned int> tetrahedronList;

    const unsigned int numPolyhedronVertices = polyhedronVertices.cols();
    for (unsigned int f = 0; f < faceTriangulations.size(); f++)
    {
      const unsigned int numFaceTriangulation = faceTriangulations[f].size() / 3;
      for (unsigned int ft = 0; ft < numFaceTriangulation; ft++)
      {
        tetrahedronList.push_back(faceTriangulations[f][3 * ft]);
        tetrahedronList.push_back(faceTriangulations[f][3 * ft + 1]);
        tetrahedronList.push_back(faceTriangulations[f][3 * ft + 2]);
        tetrahedronList.push_back(numPolyhedronVertices + polyhedronFaceInternalPoints.size());
      }
    }

    Output::Assert(tetrahedronList.size() % 4 == 0);

    return vector<unsigned int>(tetrahedronList.begin(), tetrahedronList.end());
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::ExtractTetrahedronPoints(const Eigen::MatrixXd& polyhedronVertices,
                                                               const Eigen::Vector3d& polyhedronInternalPoint,
                                                               const vector<unsigned int>& pointTetrahedrons) const
  {
    const unsigned int numTetra = pointTetrahedrons.size() / 4;
    vector<MatrixXd> tetra(numTetra, Eigen::MatrixXd::Zero(3, 4));

    for (unsigned int t = 0; t < numTetra; t++)
    {
      Eigen::MatrixXd& tetraVertices = tetra[t];
      tetraVertices.col(0)<< polyhedronVertices.col(pointTetrahedrons[4 * t]);
      tetraVertices.col(1)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 1]);
      tetraVertices.col(2)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 2]);
      tetraVertices.col(3)<< polyhedronInternalPoint;
    }

    return tetra;
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::ExtractTetrahedronPoints(const Eigen::MatrixXd& polyhedronVertices,
                                                               const Eigen::Vector3d& polyhedronInternalPoint,
                                                               const vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                               const vector<unsigned int>& pointTetrahedrons) const
  {
    const unsigned int numTetra = pointTetrahedrons.size() / 4;
    vector<MatrixXd> tetra(numTetra, Eigen::MatrixXd::Zero(3, 4));

    for (unsigned int t = 0; t < numTetra; t++)
    {
      Eigen::MatrixXd& tetraVertices = tetra[t];
      tetraVertices.col(0)<< polyhedronFaceInternalPoints[pointTetrahedrons[4 * t] - polyhedronVertices.size()];
      tetraVertices.col(1)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 1]);
      tetraVertices.col(2)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 2]);
      tetraVertices.col(3)<< polyhedronInternalPoint;
    }

    return tetra;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolyhedronCoordinateSystem(const Eigen::MatrixXd& polyhedronVertices,
                                                                     const Eigen::MatrixXi& polyhedronEdges)
  {
    vector<unsigned int> coordinateSystem(4);

    coordinateSystem[0] = 0;
    unsigned int vertexIndex = 1;
    for(unsigned int e = 0; e < polyhedronEdges.cols(); e++)
    {
      const unsigned int& originId = polyhedronEdges(0, e);
      const unsigned int& endId = polyhedronEdges(1, e);

      bool first = (originId == 0);
      bool second = (endId == 0);
      if (first || second)
        coordinateSystem[vertexIndex++] = first ? endId : originId;

      if (vertexIndex == 4)
        break;
    }

    Output::Assert(vertexIndex == 4);

    return coordinateSystem;
  }
  // ***************************************************************************
}
