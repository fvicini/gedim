#include "RefinementUtilities.hpp"
#include "VTKUtilities.hpp"


namespace Gedim
{
  // ***************************************************************************
  RefinementUtilities::RefinementUtilities(const GeometryUtilities& geometryUtilities,
                                           const MeshUtilities& meshUtilities) :
    geometryUtilities(geometryUtilities),
    meshUtilities(meshUtilities)
  {
  }
  RefinementUtilities::~RefinementUtilities()
  {
  }
  // ***************************************************************************
  RefinementUtilities::SplitCell1D_Result RefinementUtilities::SplitCell1D(const unsigned int& cell1DIndex,
                                                                           const Eigen::Vector3d& newVertexCoordinate,
                                                                           IMeshDAO& mesh) const
  {
    SplitCell1D_Result result;

    const unsigned int cell1DOriginIndex = mesh.Cell1DOrigin(cell1DIndex);
    const unsigned int cell1DEndIndex = mesh.Cell1DEnd(cell1DIndex);

    result.NewCell0DIndex = mesh.Cell0DAppend(1);
    const unsigned int& newCell0DIndex = result.NewCell0DIndex;
    mesh.Cell0DInsertCoordinates(newCell0DIndex,
                                 newVertexCoordinate);
    mesh.Cell0DSetMarker(newCell0DIndex,
                         mesh.Cell1DMarker(cell1DIndex));
    mesh.Cell0DSetState(newCell0DIndex,
                        true);

    // Split max edge into sub-edges
    Eigen::MatrixXi newCell1DsExtreme(2, 2);
    newCell1DsExtreme.col(0)<< cell1DOriginIndex, newCell0DIndex;
    newCell1DsExtreme.col(1)<< newCell0DIndex, cell1DEndIndex;
    result.NewCell1DsIndex = meshUtilities.SplitCell1D(cell1DIndex,
                                                       newCell1DsExtreme,
                                                       mesh);

    return result;
  }
  // ***************************************************************************
  bool RefinementUtilities::AreVerticesAligned(const Eigen::MatrixXd& cell2DVertices,
                                               const unsigned int fromVertex,
                                               const unsigned int toVertex) const
  {
    const unsigned int cell2DNumVertices = cell2DVertices.cols();

    // check contigous vertices
    if ((fromVertex + 1) % cell2DNumVertices == toVertex ||
        (toVertex + 1) % cell2DNumVertices == fromVertex)
    {
      return true;
    }

    // check aligned vertices
    if (geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                         cell2DVertices.col(toVertex),
                                         cell2DVertices.col((fromVertex + 1) % cell2DNumVertices)) ||
        geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                         cell2DVertices.col(toVertex),
                                         cell2DVertices.col((toVertex + 1) % cell2DNumVertices)))
    {
      return true;
    }

    return false;
  }
  // ***************************************************************************
}
