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
                                                                                          const Eigen::Vector3d& planeOrigin,
                                                                                          const Eigen::Vector3d& polygonTranslation,
                                                                                          const Eigen::Matrix3d& polygonRotationMatrix) const
  {
    SplitPolygonWithPlaneResult result;
    
    bool positiveUsed = false;
    bool negativeUsed = false;
    
    std::list<unsigned int> positiveVerticesIndex;
    std::list<unsigned int> negativeVerticesIndex;
    std::list<unsigned int> pointsOnPlaneIndex;
    std::list<Eigen::Vector3d> newVertices;
    std::list<unsigned int> newVerticesEdgeIndex;
    
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
          positiveVerticesIndex.push_back(vertexIndex);
          positiveUsed = true;
          break;
          
        case PointPlanePositionTypes::OnPlane:
          positiveVerticesIndex.push_back(vertexIndex);
          negativeVerticesIndex.push_back(vertexIndex);
          pointsOnPlaneIndex.push_back(vertexIndex);
          break;
          
        case PointPlanePositionTypes::Negative:
          negativeVerticesIndex.push_back(vertexIndex);
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
                                            pointsOnPlaneIndex.size();
        
        newVertices.push_back(polygonVertices.col(vertexIndex) +
                              edgePlaneIntersection.SingleIntersection.CurvilinearCoordinate *
                              polygonEdgeTangents.col(vertexIndex));
        newVerticesEdgeIndex.push_back(vertexIndex);
        pointsOnPlaneIndex.push_back(newVertexIndex);
        positiveVerticesIndex.push_back(newVertexIndex);
        negativeVerticesIndex.push_back(newVertexIndex);
      }
    }
    
    vector<unsigned int> positiveVertices(positiveVerticesIndex.begin(), positiveVerticesIndex.end());
    vector<unsigned int> negativeVertices(negativeVerticesIndex.begin(), negativeVerticesIndex.end());
    result.PointsOnPlane = vector<unsigned int>(pointsOnPlaneIndex.begin(), pointsOnPlaneIndex.end());
    result.NewVertices = vector<Eigen::Vector3d>(newVertices.begin(), newVertices.end());
    result.NewVerticesEdgeIndex = vector<unsigned int>(newVerticesEdgeIndex.begin(), newVerticesEdgeIndex.end());
    result.Type = (positiveUsed) ? (negativeUsed ? SplitPolygonWithPlaneResult::Types::Split :
                                                   SplitPolygonWithPlaneResult::Types::Positive) :
                                   negativeUsed ? SplitPolygonWithPlaneResult::Types::Negative :
                                                  SplitPolygonWithPlaneResult::Types::OnPlane;

    cerr<< "FACE SPLIT"<< endl;
    cerr<< "**> "<< "positiveVertices:\n"<< positiveVertices<< endl;
    cerr<< "**> "<< "negativeVertices:\n"<< negativeVertices<< endl;
    cerr<< "**> "<< "result.PointsOnPlane:\n"<< result.PointsOnPlane<< endl;
    cerr<< "**> "<< "result.NewVertices:\n"<< result.NewVertices<< endl;
    cerr<< "**> "<< "result.NewVerticesEdgeIndex:\n"<< result.NewVerticesEdgeIndex<< endl;
    cerr<< "**> "<< "result.Type:\n"<< (unsigned int)result.Type<< endl;

    // sort points with correct order
    Eigen::MatrixXd globalVertices(3, polygonVertices.cols() + result.NewVertices.size());
    globalVertices.block(0, 0, 3, polygonVertices.cols()) = polygonVertices;
    for (unsigned int v = 0; v < result.NewVertices.size(); v++)
      globalVertices.col(polygonVertices.cols() + v)<< result.NewVertices[v];

    if (positiveVertices.size() > 0)
    {
      const Eigen::MatrixXd positive3DVertices = ExtractPoints(globalVertices,
                                                               positiveVertices);
      const Eigen::MatrixXd positive2DVertices = RotatePointsFrom3DTo2D(positive3DVertices,
                                                                        polygonRotationMatrix.transpose(),
                                                                        polygonTranslation);
      const vector<unsigned int> convexHull = ConvexHull(positive2DVertices);
      Output::Assert(convexHull.size() == positiveVertices.size());

      result.PositiveVertices.resize(convexHull.size());
      for (unsigned int c = 0; c < convexHull.size(); c++)
        result.PositiveVertices[c] = positiveVertices[convexHull[c]];
    }

    if (negativeVertices.size() > 0)
    {
      const Eigen::MatrixXd negative3DVertices = ExtractPoints(globalVertices,
                                                               negativeVertices);
      const Eigen::MatrixXd negative2DVertices = RotatePointsFrom3DTo2D(negative3DVertices,
                                                                        polygonRotationMatrix.transpose(),
                                                                        polygonTranslation);

      const vector<unsigned int> convexHull = ConvexHull(negative2DVertices);
      Output::Assert(convexHull.size() == negativeVertices.size());

      result.NegativeVertices.resize(convexHull.size());
      for (unsigned int c = 0; c < convexHull.size(); c++)
        result.NegativeVertices[c] = negativeVertices[convexHull[c]];
    }

    cerr<< "***> "<< "result.PositiveVertices:\n"<< result.PositiveVertices<< endl;
    cerr<< "***> "<< "result.NegativeVertices:\n"<< result.NegativeVertices<< endl;

    return result;
  }
  // ***************************************************************************
  void GeometryUtilities::SplitPolyhedronWithPlane(const Eigen::MatrixXd& polyhedronVertices,
                                                   const Eigen::MatrixXi& polyhedronEdges,
                                                   const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                   const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                   const vector<Eigen::MatrixXd>& polyhedronFaceEdgeTangents,
                                                   const vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                   const vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices,
                                                   const Eigen::Vector3d& planeNormal,
                                                   const Eigen::Vector3d& planeOrigin,
                                                   const Eigen::Matrix3d& planeRotationMatrix,
                                                   const Eigen::Vector3d& planeTranslation) const
  {
    std::list<Eigen::Vector3d> newVertices;
    std::unordered_map<unsigned int, unsigned int> newVerticesByEdgeIndex;
    std::list<unsigned int> pointsOnPlane;
    std::list<unsigned int> positivePolyhedronVerticesIndices;
    std::list<unsigned int> negativePolyhedronVerticesIndices;
    std::list<Eigen::MatrixXi> positivePolyhedronFaces;
    std::list<Eigen::MatrixXi> negativePolyhedronFaces;
    bool positiveUsed = false, negativeUsed = false;
    
    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const unsigned int numFaceVertices = polyhedronFaces[f].cols();
      const Eigen::MatrixXd& faceVertices = polyhedronFaceVertices[f];
      Eigen::MatrixXd faceEdgeTangents = polyhedronFaceEdgeTangents[f];
      SplitPolygonWithPlaneResult splitFaceByPlane = SplitPolygonWithPlane(faceVertices,
                                                                           faceEdgeTangents,
                                                                           planeNormal,
                                                                           planeOrigin,
                                                                           polyhedronFaceTranslations[f],
                                                                           polyhedronFaceRotationMatrices[f]);
      
      switch (splitFaceByPlane.Type)
      {
        case SplitPolygonWithPlaneResult::Types::Positive:
        {
          const unsigned int& positiveFaceVertices = splitFaceByPlane.PositiveVertices.size();
          
          positivePolyhedronFaces.push_back(Eigen::MatrixXi(2, positiveFaceVertices));
          Eigen::MatrixXi& positiveFace = positivePolyhedronFaces.back();
          for (unsigned int v = 0; v < positiveFaceVertices; v++)
          {
            const unsigned int polyhedronVertexIndex = polyhedronFaces[f](0, splitFaceByPlane.PositiveVertices[v]);
            positivePolyhedronVerticesIndices.push_back(polyhedronVertexIndex);
            positiveFace(0, v) = polyhedronVertexIndex;
          }
          
          positiveUsed = true;
        }
          break;
          
        case SplitPolygonWithPlaneResult::Types::OnPlane:
        case SplitPolygonWithPlaneResult::Types::Split:
        {
          vector<unsigned int> newVerticesIndices(splitFaceByPlane.NewVertices.size());
          for (unsigned int nv = 0; nv < splitFaceByPlane.NewVertices.size(); nv++)
          {
            const unsigned int polyhedronEdgeIndex = polyhedronFaces[f](1, splitFaceByPlane.NewVerticesEdgeIndex[nv]);
            if (newVerticesByEdgeIndex.find(polyhedronEdgeIndex) == newVerticesByEdgeIndex.end())
            {
              // Add new vertex on the list
              const unsigned int polyhedronNewVertexIndex = polyhedronVertices.cols() + newVertices.size();
              newVertices.push_back(splitFaceByPlane.NewVertices[nv]);
              newVerticesByEdgeIndex.insert(make_pair(polyhedronEdgeIndex, polyhedronNewVertexIndex));
            }

            newVerticesIndices[nv] = newVerticesByEdgeIndex.at(polyhedronEdgeIndex);
          }

          const unsigned int& positiveFaceVertices = splitFaceByPlane.PositiveVertices.size();
          positivePolyhedronFaces.push_back(Eigen::MatrixXi(2, positiveFaceVertices));
          Eigen::MatrixXi& positiveFace = positivePolyhedronFaces.back();
          for (unsigned int v = 0; v < positiveFaceVertices; v++)
          {
            const unsigned int polyhedronVertexIndex = splitFaceByPlane.PositiveVertices[v] < numFaceVertices ?
                                                         polyhedronFaces[f](0, splitFaceByPlane.PositiveVertices[v]) :
              newVerticesIndices[splitFaceByPlane.PositiveVertices[v] - numFaceVertices];
              positivePolyhedronVerticesIndices.push_back(polyhedronVertexIndex);
              positiveFace(0, v) = polyhedronVertexIndex;
            }

            const unsigned int& negativeFaceVertices = splitFaceByPlane.NegativeVertices.size();
            negativePolyhedronFaces.push_back(Eigen::MatrixXi(2, negativeFaceVertices));
            Eigen::MatrixXi& negativeFace = negativePolyhedronFaces.back();
            for (unsigned int v = 0; v < negativeFaceVertices; v++)
            {
              const unsigned int polyhedronVertexIndex = splitFaceByPlane.NegativeVertices[v] < numFaceVertices ?
                                                           polyhedronFaces[f](0, splitFaceByPlane.NegativeVertices[v]) :
                newVerticesIndices[splitFaceByPlane.NegativeVertices[v] - numFaceVertices];
                negativePolyhedronVerticesIndices.push_back(polyhedronVertexIndex);
                negativeFace(0, v) = polyhedronVertexIndex;
              }
            }
              break;

        case SplitPolygonWithPlaneResult::Types::Negative:
              {
                const unsigned int& negativeFaceVertices = splitFaceByPlane.NegativeVertices.size();

                negativePolyhedronFaces.push_back(Eigen::MatrixXi(2, negativeFaceVertices));
                Eigen::MatrixXi& negativeFace = negativePolyhedronFaces.back();
                for (unsigned int v = 0; v < negativeFaceVertices; v++)
                {
                  const unsigned int polyhedronVertexIndex = polyhedronFaces[f](0, splitFaceByPlane.NegativeVertices[v]);
                  negativePolyhedronVerticesIndices.push_back(polyhedronVertexIndex);
                  negativeFace(0, v) = polyhedronVertexIndex;
                }

                negativeUsed = true;
              }
                break;
              default:
                throw runtime_error("Unsupported split polygon with plane type");
            }
          }

          cerr<< "SPLIT POLYHEDRON:\n"<< endl;
          cerr<< "*> newVertices:\n"<< newVertices<< endl;
          cerr<< "*> newVerticesByEdgeIndex:\n"<< newVerticesByEdgeIndex<< endl;
          cerr<< "*> pointsOnPlane:\n"<< pointsOnPlane<< endl;
          cerr<< "*> positivePolyhedronVerticesIndices:\n"<< positivePolyhedronVerticesIndices<< endl;
          cerr<< "*> negativePolyhedronVerticesIndices:\n"<< negativePolyhedronVerticesIndices<< endl;
          cerr<< "*> positivePolyhedronFaces:\n"<< positivePolyhedronFaces<< endl;
          cerr<< "*> negativePolyhedronFaces:\n"<< negativePolyhedronFaces<< endl;
          cerr<< "*> positiveUsed:\n"<< positiveUsed<< endl;
          cerr<< "*> negativeUsed:\n"<< negativeUsed<< endl;

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
