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
      const Eigen::MatrixXi& Face;

      VTPPolygon(const Eigen::MatrixXd& vertices,
                 const Eigen::MatrixXi& edges,
                 const Eigen::MatrixXi& face) :
        Vertices(vertices),
        Edges(edges),
        Face(face)
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

      unsigned int Size;
      const double* Data;
  };

  class IGeometryToPolyData
  {
    public:
      virtual ~IGeometryToPolyData() { }

      virtual void InitializeSolutions(const unsigned int& _numberSolutions) = 0;
      virtual void SetSolutionOptions(const unsigned int& solutionPosition,
                                      const std::string& solutionLabel,
                                      const VTPProperty::Formats& format) = 0;
      virtual void AddSolution(const unsigned int& solutionPosition,
                               const unsigned int& solutionSize,
                               const double* solution) = 0;

#if ENABLE_VTK == 1
      virtual vtkSmartPointer<vtkPolyData> Convert() const = 0;
#endif
  };

  class VTKUtilities final
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

      std::list<IGeometryToPolyData*> geometries;

#if ENABLE_VTK == 1
      vtkNew<vtkAppendFilter> exportData;
#endif

      void AddToExportData(IGeometryToPolyData& polyData,
                           const std::vector<VTPProperty>& properties);

    public:
      VTKUtilities();
      ~VTKUtilities();

      void SetExportFormat(const ExportFormat& exportFormat) { _exportFormat = exportFormat; }

      void AddPoint(const Eigen::Vector3d& point,
                    const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::MatrixXd& vertices,
                      const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::MatrixXd& vertices,
                      const Eigen::VectorXi& edge,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolygon(const Eigen::MatrixXd& vertices,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolygon(const Eigen::MatrixXd& vertices,
                      const Eigen::MatrixXi& edges,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolyhedron(const Eigen::MatrixXd& vertices,
                         const Eigen::MatrixXi& edges,
                         const std::vector<Eigen::MatrixXi>& faces,
                         const std::vector<VTPProperty>& properties = {});

      void Export(const std::string& filePath) const;
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

#if ENABLE_VTK == 1
      void AddPoint(const Eigen::Vector3d& point,
                    vtkSmartPointer<vtkPoints>& points,
                    vtkSmartPointer<vtkCellArray>& vertices) const;
      void AddLine(const unsigned int& originId,
                   const unsigned int& endId,
                   vtkSmartPointer<vtkCellArray>& lines) const;
      void AddFace(const Eigen::VectorXi& faceVerticesIds,
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

      void InitializeSolutions(const unsigned int& _numberSolutions);
      void SetSolutionOptions(const unsigned int& solutionPosition,
                              const std::string& solutionLabel,
                              const VTPProperty::Formats& format = VTPProperty::Formats::Points);
      void AddSolution(const unsigned int& solutionPosition,
                       const unsigned int& solutionSize,
                       const double* solution);

#if ENABLE_VTK == 1
      vtkSmartPointer<vtkPolyData> Convert() const;
#endif
  };
}

#endif // __VTPExtension_H
