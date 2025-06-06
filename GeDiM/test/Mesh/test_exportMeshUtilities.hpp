#ifndef __test_exportMeshUtilities_H
#define __test_exportMeshUtilities_H

#include "Eigen/Eigen"

namespace GedimUnitTesting
{
struct ExportMeshData final
{
    Eigen::MatrixXd Cell0Ds;
    Eigen::MatrixXi Cell1Ds;
    std::vector<Eigen::MatrixXi> Cell2Ds;
    std::vector<std::vector<unsigned int>> Cell3DsVertices;
    std::vector<std::vector<unsigned int>> Cell3DsEdges;
    std::vector<std::vector<unsigned int>> Cell3DsFaces;
    std::vector<unsigned int> Cell0DsMarker;
    std::vector<unsigned int> Cell1DsMarker;
    std::vector<unsigned int> Cell2DsMarker;
    std::vector<unsigned int> Cell3DsMarker;
};

struct ExportMeshGeometricData3D final
{
    std::vector<Eigen::MatrixXd> PolyhedronsVertices;
    std::vector<Eigen::MatrixXi> PolyhedronsEdges;
    std::vector<std::vector<Eigen::MatrixXi>> PolyhedronsFaces;
    std::vector<double> PolyhedronsVolume;
    std::vector<double> PolyhedronsDiameter;
    std::vector<Eigen::Vector3d> PolyhedronsCentroid;
    std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsTetrahedronsVertices;
    std::vector<std::vector<Eigen::Vector3d>> PolyhedronsFacesTranslation;
    std::vector<std::vector<Eigen::Matrix3d>> PolyhedronsFacesRotationMatrix;
    std::vector<Eigen::MatrixXd> PolyhedronsFacesNormal;
    std::vector<std::vector<bool>> PolyhedronsFacesNormalDirection;
    std::vector<std::vector<std::vector<bool>>> PolyhedronsFacesEdgesDirection;
    std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFaces2DVertices;
    std::vector<std::vector<std::vector<Eigen::Matrix3d>>> PolyhedronsFacesTriangulations2DVertices;
    std::vector<Eigen::VectorXd> PolyhedronsFacesArea;
    std::vector<std::vector<Eigen::Vector3d>> PolyhedronsFaces2DCentroid;
    std::vector<Eigen::VectorXd> PolyhedronsFacesDiameter;
    std::vector<std::vector<Eigen::VectorXd>> PolyhedronsFacesEdgesLength;
    std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DCentroid;
    std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DTangent;
    std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DTangentNormalized;
    std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DNormal;
};

struct ExportMeshGeometricData2D final
{
    std::vector<Eigen::MatrixXd> PolygonsVertices;
    std::vector<std::vector<Eigen::Matrix3d>> PolygonsTriangulations;
    std::vector<double> PolygonsArea;
    std::vector<Eigen::Vector3d> PolygonsCentroid;
    std::vector<double> PolygonsDiameter;
    std::vector<std::vector<bool>> PolygonsEdgesDirection;
    std::vector<Eigen::VectorXd> PolygonsEdgesLength;
    std::vector<Eigen::MatrixXd> PolygonsEdgesTangent;
    std::vector<Eigen::MatrixXd> PolygonsEdgesTangentNormalized;
    std::vector<Eigen::MatrixXd> PolygonsEdgesCentroid;
    std::vector<Eigen::MatrixXd> PolygonsEdgesNormal;
};

class ExportMeshUtilities
{
  public:
  public:
    ExportMeshUtilities(ExportMeshUtilities &other) = delete;
    void operator=(const ExportMeshUtilities &) = delete;

  public:
    ExportMeshUtilities()
    {
    }

    static ExportMeshData ImportMesh2DFromText(const std::string &file_path);
    static void ExportMesh2DToText(const ExportMeshData &mesh_data, const std::string &file_path);
    static ExportMeshGeometricData2D ImportMeshGeometricData2DFromText(const std::string &file_path);
    static void ExportMeshGeometricData2DToText(const ExportMeshGeometricData2D &mesh_geometric_data, const std::string &file_path);

    static ExportMeshData ImportMesh3DFromText(const std::string &file_path);
    static void ExportMesh3DToText(const ExportMeshData &mesh_data, const std::string &file_path);
    static ExportMeshGeometricData3D ImportMeshGeometricData3DFromText(const std::string &file_path);
    static void ExportMeshGeometricData3DToText(const ExportMeshGeometricData3D &mesh_geometric_data, const std::string &file_path);
};
} // namespace GedimUnitTesting

#endif // __IOUtilities_H
