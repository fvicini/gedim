#include "MeshUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  void MeshUtilities::FillMesh1D(const GeometryUtilities&,
                                 const Vector3d& segmentOrigin,
                                 const Vector3d& segmentTangent,
                                 const vector<double>& coordinates,
                                 IMeshDAO& mesh) const
  {
    if (coordinates.size() == 0)
      return;

    mesh.InitializeDimension(1);

    const unsigned int& numCell0Ds = coordinates.size();
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int c = 0; c < numCell0Ds; c++)
    {
      mesh.Cell0DSetState(c, true);
      mesh.Cell0DInsertCoordinates(c, segmentOrigin + coordinates[c] * segmentTangent);
    }
    mesh.Cell0DSetMarker(0, 1);
    mesh.Cell0DSetMarker(numCell0Ds - 1, 2);

    const unsigned int numCell1Ds = numCell0Ds - 1;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                e,
                                e + 1);
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, 0);
    }
  }
  // ***************************************************************************
  MeshUtilities::FilterMeshData MeshUtilities::FilterMesh1D(const std::vector<unsigned int>& cell1DsFilter,
                                                            const IMeshDAO& mesh) const
  {
    list<unsigned int> cell1Ds;
    std::set<unsigned int> cell0Ds;

    for (const unsigned int cell1DIndex : cell1DsFilter)
    {
      if (!mesh.Cell1DIsActive(cell1DIndex))
        continue;

      cell1Ds.push_back(cell1DIndex);

      for (unsigned int v = 0; v < 2; v++)
      {
        const unsigned int cell0DIndex = mesh.Cell1DVertex(cell1DIndex, v);

        if (cell0Ds.find(cell0DIndex) == cell0Ds.end())
          cell0Ds.insert(cell0DIndex);
      }
    }

    MeshUtilities::FilterMeshData result;

    result.Cell0Ds = std::vector<unsigned int>(cell0Ds.begin(),
                                               cell0Ds.end());
    result.Cell1Ds = std::vector<unsigned int>(cell1Ds.begin(),
                                               cell1Ds.end());

    return result;
  }
  // ***************************************************************************
  MeshUtilities::ExtractMeshData MeshUtilities::ExtractMesh1D(const std::vector<unsigned int>& cell0DsFilter,
                                                              const std::vector<unsigned int>& cell1DsFilter,
                                                              const IMeshDAO& originalMesh,
                                                              IMeshDAO& mesh) const
  {
    ExtractMeshData result;
    result.NewCell0DToOldCell0D.resize(cell0DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.NewCell1DToOldCell1D.resize(cell1DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.OldCell0DToNewCell0D.resize(originalMesh.Cell0DTotalNumber(), std::numeric_limits<unsigned int>::max());
    result.OldCell1DToNewCell1D.resize(originalMesh.Cell1DTotalNumber(), std::numeric_limits<unsigned int>::max());

    Eigen::MatrixXd newCell0Ds(3, cell0DsFilter.size());
    for (unsigned int v = 0; v < cell0DsFilter.size(); v++)
    {
      const unsigned int oldCell0DIndex = cell0DsFilter[v];
      result.NewCell0DToOldCell0D[v] = oldCell0DIndex;
      result.OldCell0DToNewCell0D[oldCell0DIndex] = v;

      newCell0Ds.col(v) = originalMesh.Cell0DCoordinates(oldCell0DIndex);
    }

    Eigen::MatrixXi newCell1Ds(2, cell1DsFilter.size());
    for (unsigned int e = 0; e < cell1DsFilter.size(); e++)
    {
      const unsigned int oldCell1DIndex = cell1DsFilter[e];
      result.NewCell1DToOldCell1D[e] = oldCell1DIndex;
      result.OldCell1DToNewCell1D[oldCell1DIndex] = e;

      const Eigen::VectorXi cell1DExtremes = originalMesh.Cell1DExtremes(oldCell1DIndex);

      newCell1Ds(0, e) = result.OldCell0DToNewCell0D.at(cell1DExtremes[0]);
      newCell1Ds(1, e) = result.OldCell0DToNewCell0D.at(cell1DExtremes[1]);
    }

    mesh.InitializeDimension(1);

    const unsigned int& numCell0Ds = newCell0Ds.cols();
    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell0DsInsertCoordinates(newCell0Ds);
    for (unsigned int c = 0; c < numCell0Ds; c++)
      mesh.Cell0DSetState(c, true);

    const unsigned int numCell1Ds = newCell1Ds.cols();
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell1DsInsertExtremes(newCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
      mesh.Cell1DSetState(e, true);

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::Mesh1DFromSegment(const GeometryUtilities& geometryUtilities,
                                        const Eigen::MatrixXd& segmentVertices,
                                        const vector<unsigned int> vertexMarkers,
                                        IMeshDAO& mesh) const
  {
    FillMesh1D(geometryUtilities,
               segmentVertices.col(0),
               segmentVertices.col(1),
               { 0.0, 1.0 },
               mesh);

    mesh.Cell0DSetMarker(0, vertexMarkers[0]);
    mesh.Cell0DSetMarker(mesh.Cell0DTotalNumber() - 1, vertexMarkers[1]);
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData1D MeshUtilities::FillMesh1DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& convexMesh) const
  {
    MeshGeometricData1D result;

    result.Cell1DsVertices.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsTangents.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsLengths.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsSquaredLengths.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsCentroids.resize(convexMesh.Cell1DTotalNumber());

    for (unsigned int c = 0; c < convexMesh.Cell1DTotalNumber(); c++)
    {
      if (!convexMesh.Cell1DIsActive(c))
        continue;

      // Extract original cell1D geometric information

      result.Cell1DsVertices[c] = convexMesh.Cell1DCoordinates(c);
      result.Cell1DsTangents[c] = geometryUtilities.SegmentTangent(result.Cell1DsVertices[c].col(0),
                                                                   result.Cell1DsVertices[c].col(1));
      result.Cell1DsLengths[c] = geometryUtilities.SegmentLength(result.Cell1DsVertices[c].col(0),
                                                                 result.Cell1DsVertices[c].col(1));
      result.Cell1DsSquaredLengths[c] = result.Cell1DsLengths[c] * result.Cell1DsLengths[c];
      result.Cell1DsCentroids[c] = geometryUtilities.SimplexBarycenter(result.Cell1DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  std::vector<unsigned int> MeshUtilities::SplitCell1D(const unsigned int& cell1DIndex,
                                                       const Eigen::MatrixXi subCell1Ds,
                                                       IMeshDAO& mesh) const
  {
    const unsigned int numSubCells = subCell1Ds.cols();
    unsigned int newCell1DsStartingIndex = mesh.Cell1DAppend(numSubCells);

    mesh.Cell1DSetState(cell1DIndex, false);

    vector<unsigned int> newCell1DsIndex(numSubCells);

    for (unsigned int c = 0; c < numSubCells; c++)
    {
      newCell1DsIndex[c] = newCell1DsStartingIndex + c;

      const unsigned int& newCell1DIndex = newCell1DsIndex[c];

      mesh.Cell1DInsertExtremes(newCell1DIndex,
                                subCell1Ds(0, c),
                                subCell1Ds(1, c));

      mesh.Cell1DSetMarker(newCell1DIndex, mesh.Cell1DMarker(cell1DIndex));
      mesh.Cell1DSetState(newCell1DIndex, true);

      mesh.Cell1DInsertUpdatedCell1D(cell1DIndex,
                                     newCell1DIndex);

      const unsigned int numCell1DNumberNeighbourCell2D = mesh.Cell1DNumberNeighbourCell2D(cell1DIndex);

      if (numCell1DNumberNeighbourCell2D > 0)
      {

        mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex,
                                              numCell1DNumberNeighbourCell2D);

        for (unsigned int n = 0; n < numCell1DNumberNeighbourCell2D; n++)
        {
          if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
            continue;

          mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                           n,
                                           mesh.Cell1DNeighbourCell2D(cell1DIndex,
                                                                      n));
        }
      }

      const unsigned int numCell1DNumberNeighbourCell3D = mesh.Cell1DNumberNeighbourCell3D(cell1DIndex);

      if (numCell1DNumberNeighbourCell3D > 0)
      {
        mesh.Cell1DInitializeNeighbourCell3Ds(newCell1DIndex,
                                              numCell1DNumberNeighbourCell3D);

        for (unsigned int n = 0; n < numCell1DNumberNeighbourCell3D; n++)
        {
          if (!mesh.Cell1DHasNeighbourCell3D(cell1DIndex, n))
            continue;

          mesh.Cell1DInsertNeighbourCell3D(newCell1DIndex,
                                           n,
                                           mesh.Cell1DNeighbourCell3D(cell1DIndex,
                                                                      n));
        }
      }
    }

    return newCell1DsIndex;
  }
  // ***************************************************************************
}
