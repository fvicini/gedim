#include "RefinementUtilities.hpp"
#include "VTKUtilities.hpp"


namespace Gedim
{
  // ***************************************************************************
  RefinementUtilities::TetrahedronMaxEdgeDirection RefinementUtilities::ComputeTetrahedronMaxEdgeDirection(const Eigen::MatrixXi& polyhedronEdges,
                                                                                                           const Eigen::VectorXd& edgesLength) const
  {
    Gedim::Output::Assert(polyhedronEdges.cols() == 6);
    Gedim::Output::Assert(edgesLength.size() == 6);

    TetrahedronMaxEdgeDirection result;

    Eigen::VectorXd::Index maxEdgeLocalIndex;
    edgesLength.maxCoeff(&maxEdgeLocalIndex);

    result.MaxEdgeIndex = maxEdgeLocalIndex;

    const unsigned int edgeOrigin = polyhedronEdges(0, result.MaxEdgeIndex);
    const unsigned int edgeEnd = polyhedronEdges(1, result.MaxEdgeIndex);

    unsigned int vertexNumber = 0, vertexFound = 0;
    while (vertexFound < 2 &&
           vertexNumber < 4)
    {
      if (vertexNumber != edgeOrigin &&
          vertexNumber != edgeEnd)
      {
        result.OppositeVerticesIndex[vertexFound++] = vertexNumber;
      }

      vertexNumber++;
    }

    Gedim::Output::Assert(vertexFound == 2);

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolyhedron_Result RefinementUtilities::RefinePolyhedronCell_ByPlane(const unsigned int& cell3DIndex,
                                                                                                 const Eigen::MatrixXd& cell3DVertices,
                                                                                                 const Eigen::MatrixXi& cell3DEdges,
                                                                                                 const std::vector<bool>& cell3DEdgesDirection,
                                                                                                 const Eigen::VectorXd& cell3DEdgesLength,
                                                                                                 const std::vector<Eigen::MatrixXi>& cell3DFaces,
                                                                                                 const std::vector<Eigen::MatrixXd>& cell3DFaces3DVertices,
                                                                                                 const std::vector<Eigen::MatrixXd>& cell3DFacesEdges3DTangent,
                                                                                                 const std::vector<Eigen::Vector3d>& cell3DFacesTranslation,
                                                                                                 const std::vector<Eigen::Matrix3d>& cell3DFacesRotationMatrix,
                                                                                                 const double& cell3DVolume,
                                                                                                 const Eigen::Vector3d& planeNormal,
                                                                                                 const Eigen::Vector3d& planeOrigin,
                                                                                                 const Eigen::Matrix3d& planeRotationMatrix,
                                                                                                 const Eigen::Vector3d& planeTranslation,
                                                                                                 IMeshDAO& mesh) const
  {
    RefinePolyhedron_Result result;

    if (mesh.Cell3DHasUpdatedCell3Ds(cell3DIndex))
    {
      result.ResultType = RefinePolyhedron_Result::ResultTypes::Cell3DAlreadySplitted;
      return result;
    }

    // Check if the cell refinement creates null cell3Ds
    if (geometryUtilities.IsValueZero(0.5 * cell3DVolume,
                                      geometryUtilities.Tolerance3D()))
    {
      result.ResultType = RefinePolyhedron_Result::ResultTypes::Cell3DSplitUnderTolerance;
      return result;
    }

    const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult split_cell = geometryUtilities.SplitPolyhedronWithPlane(cell3DVertices,
                                                                                                                           cell3DEdges,
                                                                                                                           cell3DFaces,
                                                                                                                           cell3DFaces3DVertices,
                                                                                                                           cell3DFacesEdges3DTangent,
                                                                                                                           cell3DFacesTranslation,
                                                                                                                           cell3DFacesRotationMatrix,
                                                                                                                           planeNormal,
                                                                                                                           planeOrigin,
                                                                                                                           planeRotationMatrix,
                                                                                                                           planeTranslation);

    // Create new mesh elements
    const unsigned int numOriginalVertices = cell3DVertices.cols();
    const unsigned int numNewVertices = split_cell.Vertices.NewVerticesOriginalEdge.size();

    std::vector<unsigned int> splitCell0DsIndex(split_cell.Vertices.Vertices.cols());
    std::vector<unsigned int> splitCell1DsIndex(split_cell.Edges.NewEdgesOriginalEdges.size());
    std::vector<unsigned int> splitCell2DsIndex(split_cell.Faces.NewFacesOriginalFaces.size());

    std::list<RefinePolyhedron_Result::RefinedCell1D> newCell1Ds;

    for (unsigned int v = 0; v < numOriginalVertices; v++)
      splitCell0DsIndex[v] = mesh.Cell3DVertex(cell3DIndex, v);

    // Create new vertices
    result.NewCell0DsIndex.resize(numNewVertices);
    for (unsigned int nv = 0; nv < numNewVertices; nv++)
    {
      const unsigned int originalEdgeIndex = split_cell.Vertices.NewVerticesOriginalEdge.at(nv);
      const unsigned int originalCell1DIndex = mesh.Cell3DEdge(cell3DIndex,
                                                               originalEdgeIndex);

      const SplitCell1D_Result splitResult = SplitCell1D(originalCell1DIndex,
                                                         split_cell.Vertices.Vertices.col(numOriginalVertices + nv),
                                                         mesh);

      result.NewCell0DsIndex[nv] = splitResult.NewCell0DIndex;
      splitCell0DsIndex[numOriginalVertices + nv] = splitResult.NewCell0DIndex;

      RefinePolyhedron_Result::RefinedCell1D refinedCell1D;
      refinedCell1D.Type = RefinePolyhedron_Result::RefinedCell1D::Types::Updated;
      refinedCell1D.NewCell1DsIndex = splitResult.NewCell1DsIndex;
      refinedCell1D.OriginalCell1DIndex = originalCell1DIndex;
      refinedCell1D.NewCell0DIndex = splitResult.NewCell0DIndex;
      refinedCell1D.OriginalCell3DEdgeIndex = originalEdgeIndex;
      newCell1Ds.push_back(refinedCell1D);
    }

    // Create new edges
    for (unsigned int ne = 0; ne < split_cell.Edges.NewEdgesOriginalEdges.size(); ne++)
    {
      if (split_cell.Edges.NewEdgesOriginalEdges[ne] >= 0)
      {
        splitCell1DsIndex[ne] = mesh.Cell3DEdge(cell3DIndex,
                                                split_cell.Edges.NewEdgesOriginalEdges[ne]);
        continue;
      }

      const unsigned int originalFaceIndex = split_cell.Edges.NewEdgesOriginalFace.at(ne);
      const unsigned int originalCell2DIndex = mesh.Cell3DFace(cell3DIndex,
                                                               originalFaceIndex);

      const unsigned int newEdgeOrigin = split_cell.Edges.Edges(0, ne);
      const unsigned int newEdgeEnd = split_cell.Edges.Edges(1, ne);

      const unsigned int cell0DIndexOrigin = splitCell0DsIndex.at(newEdgeOrigin);
      const unsigned int cell0DIndexEnd = splitCell0DsIndex.at(newEdgeEnd);

      RefinePolyhedron_Result::RefinedCell1D newCell1D;
      newCell1D.Type = RefinePolyhedron_Result::RefinedCell1D::Types::New;
      newCell1D.NewCell1DsIndex.resize(1);

      newCell1D.NewCell1DsIndex[0] = mesh.Cell1DAppend(1);
      const unsigned int newCell1DIndex = newCell1D.NewCell1DsIndex[0];
      splitCell1DsIndex[ne] = newCell1DIndex;

      mesh.Cell1DInsertExtremes(newCell1DIndex,
                                cell0DIndexEnd,
                                cell0DIndexOrigin);
      mesh.Cell1DSetMarker(newCell1DIndex,
                           mesh.Cell2DMarker(originalCell2DIndex));
      mesh.Cell1DSetState(newCell1DIndex, true);
      newCell1Ds.push_back(newCell1D);
    }

    // Create new Faces
    for (unsigned int nf = 0; nf < split_cell.Faces.NewFacesOriginalFaces.size(); nf++)
    {
      if (split_cell.Faces.NewFacesOriginalFaces[nf] >= 0)
        continue;

      const Eigen::MatrixXi newFace = split_cell.Faces.Faces[nf];

      Eigen::MatrixXi cell2DExtremes(2, newFace.cols());

      for (unsigned int nfv = 0; nfv < newFace.cols(); nfv++)
      {
        cell2DExtremes(0, nfv) = splitCell0DsIndex.at(newFace(0, nfv));
        cell2DExtremes(1, nfv) = splitCell1DsIndex.at(newFace(1, nfv));
      }

      RefinePolyhedron_Result::RefinedCell2D newCell2D;
      newCell2D.Type = RefinePolyhedron_Result::RefinedCell2D::Types::New;
      newCell2D.NewCell2DsIndex.resize(1);
      newCell2D.NewCell2DsIndex[0] = mesh.Cell2DAppend(1);
      const unsigned int newCell2DIndex = newCell2D.NewCell2DsIndex[0];
      splitCell2DsIndex[nf] = newCell2DIndex;
      mesh.Cell2DAddVerticesAndEdges(newCell2DIndex,
                                     cell2DExtremes);

      mesh.Cell2DSetMarker(newCell2DIndex, 0);
      mesh.Cell2DSetState(newCell2DIndex, true);
      mesh.Cell2DInitializeNeighbourCell3Ds(newCell2DIndex,
                                            2);
      // TODO: add neighbours
    }

    result.NewCell1DsIndex = std::vector<RefinePolyhedron_Result::RefinedCell1D>(newCell1Ds.begin(),
                                                                                 newCell1Ds.end());

    return result;
  }
  // ***************************************************************************
}
