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
  struct IVTPGeometry
  {
      enum Types
      {
        Point = 0,
        Segment = 1,
        Polygon = 2,
        Polyhedron = 3
      };

      virtual ~IVTPGeometry() {}

      virtual Types Type() const = 0;
  };

  struct VTPPoint final : public IVTPGeometry
  {
      const Eigen::Vector3d& Point;

      VTPPoint(const Eigen::Vector3d& point) :
        Point(point)
      { }

      Types Type() const { return Types::Point; }
  };

  struct VTPSegment final : public IVTPGeometry
  {
      const Eigen::MatrixXd& Vertices;
      const Eigen::VectorXi& Edge;

      VTPSegment(const Eigen::MatrixXd& vertices,
                 const Eigen::VectorXi& edge) :
        Vertices(vertices),
        Edge(edge)
      { }

      Types Type() const { return Types::Segment; }
  };

  struct VTPPolygon final : public IVTPGeometry
  {
      const Eigen::MatrixXd& Vertices;
      const Eigen::MatrixXi& Edges;

      VTPPolygon(const Eigen::MatrixXd& vertices,
                 const Eigen::MatrixXi& edges) :
        Vertices(vertices),
        Edges(edges)
      { }

      Types Type() const { return Types::Polygon; }
  };

  struct VTPPolyhedron final : public IVTPGeometry
  {
      const Eigen::MatrixXd& Vertices; /// size 3xnumVertices
      const Eigen::MatrixXi& Edges; /// size 2xnumEdges
      const std::vector<Eigen::MatrixXi>& Faces; /// size numFacesx2xnumFaceVertices

      VTPPolyhedron(const Eigen::MatrixXd& vertices,
                    const Eigen::MatrixXi& edges,
                    const std::vector<Eigen::MatrixXi>& faces) :
        Vertices(vertices),
        Edges(edges),
        Faces(faces)
      { }

      Types Type() const { return Types::Polyhedron; }
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
  class VTPUtilities final
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
      std::list<IVTPGeometry*> geometries;

      std::unordered_map<std::string, VTPProperty> properties;
    public:
      VTPUtilities();
      ~VTPUtilities();

      void SetExportFormat(const ExportFormat& exportFormat) { _exportFormat = exportFormat; }

      void AddProperty(const std::string& solutionLabel = "Property",
                       const VTPProperty::Formats& format = VTPProperty::Formats::Points);

      inline void AddPoint(const Eigen::Vector3d& point)
      {
        geometries.push_back(new VTPPoint(point));
      }

      inline void AddSegment(const Eigen::MatrixXd& vertices,
                             const Eigen::VectorXi& edge)
      {
        geometries.push_back(new VTPSegment(vertices,
                                            edge));
      }

      inline void AddPolygon(const Eigen::MatrixXd& vertices,
                             const Eigen::MatrixXi& edges)
      {
        geometries.push_back(new VTPPolygon(vertices,
                                            edges));
      }

      inline void AddPolyhedron(const Eigen::MatrixXd& vertices,
                                const Eigen::MatrixXi& edges,
                                const std::vector<Eigen::MatrixXi>& faces)
      {
        geometries.push_back(new VTPPolyhedron(vertices,
                                               edges,
                                               faces));
      }

      void AddGeometryProperty(const std::string& solutionLabel,
                               const unsigned int& solutionSize,
                               const double* solution);

      void Export(const std::string& filePath) const;
  };

  class IGeometryToPolyData
  {
    public:
      virtual ~IGeometryToPolyData() { }

      virtual void InitializeSolutions(const unsigned int& _numberSolutions) = 0;
      virtual void SetSolutionOptions(const unsigned int& solutionPosition,
                                      const std::string& solutionLabel,
                                      const VTPProperty::Formats& format = VTPProperty::Formats::Points) = 0;
      virtual void AddSolution(const unsigned int& solutionPosition,
                               const unsigned int& solutionSize,
                               const double* solution) = 0;

      virtual void Convert() = 0;

      virtual const void* PolyDataPointer() const = 0;
  };

  template <typename T>
  class GeometryToPolyData final : public IGeometryToPolyData
  {
    protected:
      const T& geometry;

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
      GeometryToPolyData(const T& geometry);
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
