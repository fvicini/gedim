#ifndef __VTPExtension_H
#define __VTPExtension_H

#include "Macro.hpp"

#include <vector>
#include <map>

#include "Eigen/Eigen"

#if ENABLE_VTK == 1
#include "vtkCellType.h"
#include "vtkSmartPointer.h"
#include "vtkWeakPointer.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIdList.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkCellTypes.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkAppendFilter.h"
#endif // ENABLE_VTK


namespace Gedim
{
  struct VTPGeometry
  {
      enum Types
      {
        Point = 0,
        Segment = 1,
        Polygon = 2,
        Polyhedron = 3
      };

      Types Type;

      VTPGeometry(Types type) :
        Type(type)
      {}

      virtual ~VTPGeometry() {}
  };

  struct VTPPoint final : VTPGeometry
  {
      const Eigen::Vector3d& Point;

      VTPPoint(const Eigen::Vector3d& point) :
        VTPGeometry(Types::Point),
        Point(point)
      { }
  };

  struct VTPSegment final : VTPGeometry
  {
      const Eigen::MatrixXd& Vertices;
      const Eigen::VectorXi& Edge;

      VTPSegment(const Eigen::MatrixXd& vertices,
                 const Eigen::VectorXi& edge) :
        VTPGeometry(Types::Segment),
        Vertices(vertices),
        Edge(edge)
      { }
  };

  struct VTPPolygon final : VTPGeometry
  {
      const Eigen::MatrixXd& Vertices;
      const Eigen::MatrixXi& Edges;

      VTPPolygon(const Eigen::MatrixXd& vertices,
                 const Eigen::MatrixXi& edges) :
        VTPGeometry(Types::Polygon),
        Vertices(vertices),
        Edges(edges)
      { }
  };

  struct VTPPolyhedron final : VTPGeometry
  {
      const Eigen::MatrixXd& Vertices; /// size 3xnumVertices
      const Eigen::MatrixXi& Edges; /// size 2xnumEdges
      const std::vector<Eigen::MatrixXi>& Faces; /// size numFacesx2xnumFaceVertices

      VTPPolyhedron(const Eigen::MatrixXd& vertices,
                    const Eigen::MatrixXi& edges,
                    const std::vector<Eigen::MatrixXi>& faces) :
        VTPGeometry(Types::Polyhedron),
        Vertices(vertices),
        Edges(edges),
        Faces(faces)
      { }
  };

  struct VTPProperty final
  {
      enum Formats
      {
        Points = 0,
        Cells = 1
      };

      std::string Label;
      Formats Format;

      std::list<unsigned int> Size;
      std::list<const double*> Data;
  };

  /// \brief The VTPExtension class
  /// \copyright See top level LICENSE file for details.
  class VTPExtension final
  {
    public:
      enum ExportFormat
      {
        Binary = 0,
        Ascii = 1,
        Appended = 2
      };

    private:
      ExportFormat _exportFormat;
      std::list<VTPGeometry> geometries;

      std::list<VTPProperty> properties;
    public:
      VTPExtension();
      ~VTPExtension();

      void SetExportFormat(const ExportFormat& exportFormat) { _exportFormat = exportFormat; }

      void AddProperty(const std::string& solutionLabel = "Property",
                       const VTPProperty::Formats& format = VTPProperty::Formats::Points);

      inline void AddPoint(const Eigen::Vector3d& point)
      {
        geometries.push_back(VTPPoint(point));
      }

      inline void AddSegment(const Eigen::MatrixXd& vertices,
                             const Eigen::VectorXi& edge)
      {
        geometries.push_back(VTPSegment(vertices,
                                        edge));
      }

      inline void AddPolygon(const Eigen::MatrixXd& vertices,
                             const Eigen::MatrixXi& edges)
      {
        geometries.push_back(VTPPolygon(vertices,
                                        edges));
      }

      inline void AddPolyhedron(const Eigen::MatrixXd& vertices,
                                const Eigen::MatrixXi& edges,
                                const std::vector<Eigen::MatrixXi>& faces)
      {
        geometries.push_back(VTPPolyhedron(vertices,
                                           edges,
                                           faces));
      }

      void AddGeometrySolution(const unsigned int& solutionSize,
                               const double* solution);

      void Export(const std::string& filePath) const;
  };

  class GeometryToPolyData final
  {
    protected:
      const VTPGeometry& geometry;

      unsigned int numberSolutions;
      std::vector<VTPProperty::Formats> solutionFormats;
      std::vector<std::string> solutionLabels;
      std::vector<unsigned int> solutionSizes;
      std::vector<const double*> solutions;

      void* vtkPointsPointer;
      void* vtkVerticesPointer;
      void* vtkLinesPointer;
      void* vtkFacesPointer;
      void* vtkPolyDataPointer;
      std::vector<void*> vtkDoubleArrayPointers;

#if ENABLE_VTK == 1
      void AddPoint(const Eigen::Vector3d& point,
                    vtkSmartPointer<vtkPoints>& points,
                    vtkSmartPointer<vtkCellArray>& vertices) const;
      void AddLine(const unsigned int& originId,
                   const unsigned int& endId,
                   vtkSmartPointer<vtkIdList>& pointIds,
                   vtkSmartPointer<vtkCellArray>& lines) const;
      void AddFace(const Eigen::VectorXi& faceVerticesIds,
                   vtkSmartPointer<vtkIdList>& pointIds,
                   vtkSmartPointer<vtkCellArray>& faces) const;
#endif

      void InitializePoints(void* vtkPointsPointer, void* vtkVerticesPointer) const;
      void AddPoints(void* vtkPointsPointer, void* vtkVerticesPointer) const;
      void InitializeLines(void* vtkLinesPointer) const;
      void AddLines(void* vtkLinesPointer) const;
      void InitializeFaces(void* vtkFacesPointer) const;
      void AddFaces(void* vtkFacesPointer) const;
      void AddSolution(const unsigned int& solutionPosition, void* vtkDoubleArrayPointer) const;

    public:
      GeometryToPolyData(const VTPGeometry& geometry);
      ~GeometryToPolyData();

      const void* PolyDataPointer() const { return vtkPolyDataPointer; }

      void InitializeSolutions(const unsigned int& _numberSolutions);
      void SetSolutionOptions(const unsigned int& solutionPosition,
                              const std::string& solutionLabel,
                              const VTPProperty::Formats& format = VTPProperty::Formats::Points);
      void AddSolution(const unsigned int& solutionPosition,
                       const unsigned int& solutionSize,
                       const double* solution);

      void Convert();
  };
}

#endif // __VTPExtension_H
