#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::SplitPolygonWithPlaneResult GeometryUtilities::SplitPolygonWithPlane(const Eigen::MatrixXd& polygonVertices,
                                                                                          const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                          const Eigen::Vector3d& planeNormal,
                                                                                          const Eigen::Vector3d& planeOrigin) const
  {
    SplitPolygonWithPlaneResult result;

    bool positiveUsed = false;
    bool negativeUsed = false;

    std::list<unsigned int> positiveVertices;
    std::list<unsigned int> negativeVertices;
    std::list<unsigned int> pointsOnPlane;
    std::list<Eigen::Vector3d> newVertices;

    // note: for convex polygons only
    for (unsigned int i = 0; i < polygonVertices.cols(); i++)
    {
      // Place the vertex in the front or back polygon (or both)
      const unsigned int& vertexIndex = i;
      PointPlanePositionTypes vertexPlanePosition = PointPlanePosition(polygonVertices.col(vertexIndex),
                                                                       planeNormal,
                                                                       planeOrigin);
      switch (vertexPlanePosition)
      {
        case PointPlanePositionTypes::Positive:
          positiveVertices.push_back(vertexIndex);
          positiveUsed = true;
          break;

        case PointPlanePositionTypes::OnPlane:
          positiveVertices.push_back(vertexIndex);
          negativeVertices.push_back(vertexIndex);
          pointsOnPlane.push_back(vertexIndex);
          break;

        case PointPlanePositionTypes::Negative:
          negativeVertices.push_back(vertexIndex);
          negativeUsed = true;
          break;
        default:
          throw runtime_error("Unsupported plane position");
      }

      // Add a vertex to both polygons where edges intersect the plane
      unsigned int nextVertexIndex = (i + 1) % polygonVertices.cols();

      PointPlanePositionTypes nextVertexPlanePosition = PointPlanePosition(polygonVertices.col(nextVertexIndex),
                                                                           planeNormal,
                                                                           planeOrigin);
      if ((vertexPlanePosition == PointPlanePositionTypes::Positive &&
           nextVertexPlanePosition == PointPlanePositionTypes::Negative) ||
          (vertexPlanePosition == PointPlanePositionTypes::Negative &&
           nextVertexPlanePosition == PointPlanePositionTypes::Positive))
      {
        IntersectionSegmentPlaneResult edgePlaneIntersection = IntersectionSegmentPlane(polygonVertices.col(vertexIndex),
                                                                                        polygonVertices.col(nextVertexIndex),
                                                                                        planeNormal,
                                                                                        planeOrigin);
        Output::Assert(edgePlaneIntersection.Type == IntersectionSegmentPlaneResult::Types::SingleIntersection &&
                       edgePlaneIntersection.SingleIntersection.Type == PointSegmentPositionTypes::InsideSegment);

        const unsigned int newVertexIndex = polygonVertices.cols() +
                                            pointsOnPlane.size();

        newVertices.push_back(polygonVertices.col(vertexIndex) +
                              edgePlaneIntersection.SingleIntersection.CurvilinearCoordinate *
                              polygonEdgeTangents.col(vertexIndex));
        pointsOnPlane.push_back(newVertexIndex);
        positiveVertices.push_back(newVertexIndex);
        negativeVertices.push_back(newVertexIndex);
      }
    }

    result.PositiveVertices = vector<unsigned int>(positiveVertices.begin(), positiveVertices.end());
    result.NegativeVertices = vector<unsigned int>(negativeVertices.begin(), negativeVertices.end());
    result.PointsOnPlane = vector<unsigned int>(pointsOnPlane.begin(), pointsOnPlane.end());
    result.NewVertices = vector<Eigen::Vector3d>(newVertices.begin(), newVertices.end());
    result.Type = (positiveUsed) ? (negativeUsed ? SplitPolygonWithPlaneResult::Types::Split :
                                                   SplitPolygonWithPlaneResult::Types::Positive) :
                                   negativeUsed ? SplitPolygonWithPlaneResult::Types::Negative :
                                                  SplitPolygonWithPlaneResult::Types::OnPlane;

    return result;
  }
  // ***************************************************************************
  void GeometryUtilities::SplitPolyhedronByPlane(const Eigen::MatrixXd& polyhedronVertices,
                                                 const Eigen::MatrixXi& polyhedronEdges,
                                                 const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                 const Eigen::MatrixXd& polyhedronEdgeTangents,
                                                 const Eigen::Vector3d& planeNormal,
                                                 const Eigen::Vector3d& planeOrigin) const
  {
    std::list<Eigen::Vector3d> newVertices;
    std::list<unsigned int> pointsOnPlane;
    std::list<unsigned int> positivePolyhedronVerticesIndices;
    std::list<unsigned int> negativePolyhedronVerticesIndices;
    std::list<Eigen::MatrixXi> positivePolyhedronFaces;
    std::list<Eigen::MatrixXi> negativePolyhedronFaces;
    bool positiveUsed = false, negativeUsed = false;

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const VectorXi faceVerticesIds = polyhedronFaces[f].row(0);
      Eigen::MatrixXd faceVertices = ExtractPoints(polyhedronVertices,
                                                   vector<unsigned int>(faceVerticesIds.data(),
                                                                        faceVerticesIds.data() +
                                                                        faceVerticesIds.size()));
      Eigen::MatrixXd faceEdgeTangents = ExtractPoints(polyhedronEdgeTangents,
                                                       vector<unsigned int>(faceVerticesIds.data(),
                                                                            faceVerticesIds.data() +
                                                                            faceVerticesIds.size()));
      SplitPolygonWithPlaneResult splitFaceByPlane = SplitPolygonWithPlane(faceVertices,
                                                                           faceEdgeTangents,
                                                                           planeNormal,
                                                                           planeOrigin);

      switch (splitFaceByPlane.Type)
      {
        case SplitPolygonWithPlaneResult::Types::Positive:
        {
          const unsigned int& positiveFaceVertices = splitFaceByPlane.PositiveVertices.size();

          positivePolyhedronFaces.push_back(Eigen::MatrixXi(2, positiveFaceVertices));
          Eigen::MatrixXi& positiveFace = positivePolyhedronFaces.back();
          for (unsigned int v = 0; v < positiveFaceVertices; v++)
          {
            positivePolyhedronVerticesIndices.push_back(polyhedronFaces[f](0, splitFaceByPlane.PositiveVertices[v]));
                positiveFace(0, v) = polyhedronFaces[f](0, splitFaceByPlane.PositiveVertices[v]);
          }

          positiveUsed = true;
        }
          break;

        case SplitPolygonWithPlaneResult::Types::OnPlane:
        case SplitPolygonWithPlaneResult::Types::Split:
        {
          // TODO: create new split faces
          //          frontPolyhedron->m_polygons.push_back(front);
          //          backPolyhedron->m_polygons.push_back(back);
        }
          break;

        case SplitPolygonWithPlaneResult::Types::Negative:
        {
          const unsigned int& negativeFaceVertices = splitFaceByPlane.NegativeVertices.size();

          negativePolyhedronFaces.push_back(Eigen::MatrixXi(2, negativeFaceVertices));
          Eigen::MatrixXi& negativeFace = negativePolyhedronFaces.back();
          for (unsigned int v = 0; v < negativeFaceVertices; v++)
          {
            negativePolyhedronVerticesIndices.push_back(polyhedronFaces[f](0, splitFaceByPlane.NegativeVertices[v]));
                negativeFace(0, v) = polyhedronFaces[f](0, splitFaceByPlane.NegativeVertices[v]);
          }

          negativeUsed = true;
        }
          break;
        default:
          throw runtime_error("Unsupported split polygon with plane type");
      }
    }

    // TODO: create new polyhedrons
    //    result.clear();
    //    if (frontUsed == backUsed)
    //    {
    //      // Create polygons to fill in the holes inside the new polyhedra
    //      frontPolyhedron->m_polygons.push_back(Polygon::fromPoints(pointsOnPlane, plane.normal, frontPolyhedron->m_vertices));
    //      backPolyhedron->m_polygons.push_back(Polygon::fromPoints(pointsOnPlane, -plane.normal, backPolyhedron->m_vertices));
    //      result.push_back(frontPolyhedron);
    //      result.push_back(backPolyhedron);
    //      return true;
    //    }
    //    else
    //    {
    //      delete frontPolyhedron;
    //      delete backPolyhedron;
    //      return false;
    //    }
  }
  // ***************************************************************************
}
