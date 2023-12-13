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

    return result;
  }
  // ***************************************************************************
}
