#include "VTPUtilities.hpp"

#include "IOUtilities.hpp"
#include "CommonUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  VTKUtilities::VTKUtilities()
  {
  }
  VTKUtilities::~VTKUtilities()
  {
  }
  // ***************************************************************************
  void VTKUtilities::AddToExportData(IGeometryToPolyData& polyData,
                                     const vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    if (properties.size() > 0)
    {
      polyData.InitializeSolutions(properties.size());
      unsigned int p = 0;
      for (const VTPProperty& property : properties)
      {
        polyData.SetSolutionOptions(p,
                                    property.Label,
                                    property.Format);
        polyData.AddSolution(p,
                             property.Size,
                             property.Data);
        p++;
      }
    }

    exportData->AddInputData(polyData.Convert());
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddPoint(const Eigen::Vector3d& point,
                              const vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    VTPPoint vtpPoint(point);
    GeometryToPolyData<VTPPoint> polyData(vtpPoint);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddSegment(const Eigen::MatrixXd& vertices,
                                const std::vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    Eigen::MatrixXi edge =  (Eigen::MatrixXi(2,1)<< 0, 1).finished();
    VTPSegment vtpSegment(vertices,
                          edge);
    GeometryToPolyData<VTPSegment> polyData(vtpSegment);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddPolygon(const Eigen::MatrixXd& vertices,
                                const std::vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    Eigen::MatrixXi edges(2, vertices.cols());
    for (unsigned int v = 0; v < vertices.cols(); v++)
      edges.col(v)<< v, (v + 1) % vertices.cols();

    Eigen::MatrixXi face(2, vertices.cols());
    face.row(0)<< Eigen::VectorXi::LinSpaced(vertices.cols(),
                                             0.0,
                                             vertices.cols() - 1).transpose();
    face.row(1)<< Eigen::VectorXi::LinSpaced(vertices.cols(),
                                             0.0,
                                             vertices.cols() - 1).transpose();

    VTPPolygon vtpPolygon(vertices,
                          edges,
                          face);
    GeometryToPolyData<VTPPolygon> polyData(vtpPolygon);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddPolyhedron(const Eigen::MatrixXd& vertices,
                                   const Eigen::MatrixXi& edges,
                                   const std::vector<Eigen::MatrixXi>& faces,
                                   const std::vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    VTPPolyhedron vtpPolyhedron(vertices,
                                edges,
                                faces);
    GeometryToPolyData<VTPPolyhedron> polyData(vtpPolyhedron);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::Export(const std::string& filePath) const
  {
#if ENABLE_VTK == 1
    exportData->Update();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filePath.c_str());
    writer->SetInputData(exportData->GetOutput());

    switch (_exportFormat)
    {
      case VTKUtilities::Binary:
        writer->SetDataModeToBinary();
        break;
      case VTKUtilities::Ascii:
        writer->SetDataModeToAscii();
        break;
      case VTKUtilities::Appended:
        writer->SetDataModeToAppended();
        break;
      default:
        throw runtime_error("Export Format not supported");
    }

    writer->Write();
#endif
  }
  // ***************************************************************************
  template class GeometryToPolyData<VTPPoint>;
  template class GeometryToPolyData<VTPSegment>;
  template class GeometryToPolyData<VTPPolygon>;
  template class GeometryToPolyData<VTPPolyhedron>;
  // ***************************************************************************
  template <typename T>
  GeometryToPolyData<T>::GeometryToPolyData(const T& geometry) :
    geometry(geometry)
  {
    numberSolutions = 0;
  }
  template <typename T>
  GeometryToPolyData<T>::~GeometryToPolyData()
  {
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::InitializeSolutions(const unsigned int& _numberSolutions)
  {
    Output::Assert(_numberSolutions > 0);
    numberSolutions = _numberSolutions;

    solutionLabels.resize(numberSolutions, "Solution");
    solutionFormats.resize(numberSolutions, VTPProperty::Formats::Points);
    solutionSizes.resize(numberSolutions, 0);
    solutions.resize(numberSolutions, NULL);
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::SetSolutionOptions(const unsigned int& solutionPosition,
                                                 const string& solutionLabel,
                                                 const VTPProperty::Formats& format)
  {
    Output::Assert(solutionPosition < numberSolutions);

    solutionLabels[solutionPosition] = solutionLabel;
    solutionFormats[solutionPosition] = format;
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddSolution(const unsigned int& solutionPosition,
                                          const unsigned int& solutionSize,
                                          const double* solution)
  {
    Output::Assert(solutionPosition < solutionLabels.size());

    if (solutions[solutionPosition] != NULL)
      throw runtime_error("Solution " + to_string(solutionPosition) + " already inserted");

    solutions[solutionPosition] = solution;
    solutionSizes[solutionPosition] = solutionSize;
  }
  // ***************************************************************************
#if ENABLE_VTK == 1
  template <typename T>
  void GeometryToPolyData<T>::AddPoint(const Eigen::Vector3d& point,
                                       vtkSmartPointer<vtkPoints>& points,
                                       vtkSmartPointer<vtkCellArray>& vertices) const
  {
    vtkIdType pointId[1];
    pointId[0] = points->InsertNextPoint(point(0), point(1), point(2));

    vertices->InsertNextCell(1, pointId);
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddLine(const unsigned int& originId,
                                      const unsigned int& endId,
                                      vtkSmartPointer<vtkCellArray>& lines) const
  {
    vtkNew<vtkIdList> pointIds;

    pointIds->Initialize();
    pointIds->Allocate(2);

    pointIds->InsertNextId(originId);
    pointIds->InsertNextId(endId);

    pointIds->Squeeze();

    lines->InsertNextCell(pointIds);

    pointIds->Reset();
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddFace(const Eigen::VectorXi& faceVerticesIds,
                                      vtkSmartPointer<vtkIdList>& pointIds,
                                      vtkSmartPointer<vtkCellArray>& faces) const
  {
    pointIds->Initialize();
    pointIds->Allocate(faceVerticesIds.size());

    for (unsigned int fp = 0 ; fp < faceVerticesIds.size(); fp++)
      pointIds->InsertNextId(faceVerticesIds[fp]);

    pointIds->Squeeze();

    faces->InsertNextCell(pointIds);

    pointIds->Reset();
  }
#endif
  // ***************************************************************************
  template<>
  void GeometryToPolyData<VTPPoint>::InitializePoints(void* vtkPointsPointer,
                                                      void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    const unsigned int numberTotalPoints = 1;
    points->Allocate(numberTotalPoints);
    vertices->Allocate(numberTotalPoints);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template<>
  void GeometryToPolyData<VTPSegment>::InitializePoints(void* vtkPointsPointer,
                                                        void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    const unsigned int numberTotalPoints = 2;
    points->Allocate(numberTotalPoints);
    vertices->Allocate(numberTotalPoints);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template<>
  void GeometryToPolyData<VTPPolygon>::InitializePoints(void* vtkPointsPointer,
                                                        void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    const unsigned int numberTotalPoints = geometry.Vertices.cols();
    points->Allocate(numberTotalPoints);
    vertices->Allocate(numberTotalPoints);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template<>
  void GeometryToPolyData<VTPPolyhedron>::InitializePoints(void* vtkPointsPointer,
                                                           void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    const unsigned int numberTotalPoints = geometry.Vertices.cols();
    points->Allocate(numberTotalPoints);
    vertices->Allocate(numberTotalPoints);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPoint>::AddPoints(void* vtkPointsPointer,
                                               void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    AddPoint(geometry.Point,
             points,
             vertices);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPSegment>::AddPoints(void* vtkPointsPointer,
                                                 void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    AddPoint(geometry.Vertices.col(0),
             points,
             vertices);
    AddPoint(geometry.Vertices.col(1),
             points,
             vertices);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolygon>::AddPoints(void* vtkPointsPointer,
                                                 void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    for (unsigned int i = 0 ; i < geometry.Vertices.cols(); i++)
    {
      AddPoint(geometry.Vertices.col(i),
               points,
               vertices);
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolyhedron>::AddPoints(void* vtkPointsPointer,
                                                    void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    for (unsigned int i = 0 ; i < geometry.Vertices.cols(); i++)
    {
      AddPoint(geometry.Vertices.col(i),
               points,
               vertices);
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPoint>::InitializeLines(void* vtkLinesPointer) const
  {
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPSegment>::InitializeLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;

    unsigned int numberTotalLines = 1;
    lines->Allocate(numberTotalLines);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolygon>::InitializeLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;

    unsigned int numberTotalLines = geometry.Vertices.cols();
    lines->Allocate(numberTotalLines);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolyhedron>::InitializeLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;

    unsigned int numberTotalLines = geometry.Vertices.cols();
    lines->Allocate(numberTotalLines);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPoint>::AddLines(void* vtkLinesPointer) const
  {
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPSegment>::AddLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;

    AddLine(geometry.Edge[0],
        geometry.Edge[1],
        lines);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolygon>::AddLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;

    for (unsigned int i = 0 ; i < geometry.Edges.cols(); i++)
    {
      AddLine(geometry.Edges(0, i),
              geometry.Edges(1, i),
              lines);
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolyhedron>::AddLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;
    for (unsigned int i = 0 ; i < geometry.Edges.cols(); i++)
    {
      AddLine(geometry.Edges(0, i),
              geometry.Edges(1, i),
              lines);
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPoint>::InitializeFaces(void* vtkFacesPointer) const
  {
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPSegment>::InitializeFaces(void* vtkFacesPointer) const
  {
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolygon>::InitializeFaces(void* vtkFacesPointer) const
  {
    if (vtkFacesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> faces = *(vtkSmartPointer<vtkCellArray>*)vtkFacesPointer;

    unsigned int numberTotalFaces = 1;
    faces->Allocate(numberTotalFaces);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolyhedron>::InitializeFaces(void* vtkFacesPointer) const
  {
    if (vtkFacesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> faces = *(vtkSmartPointer<vtkCellArray>*)vtkFacesPointer;

    unsigned int numberTotalFaces = geometry.Faces.size();
    faces->Allocate(numberTotalFaces);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPoint>::AddFaces(void* vtkFacesPointer) const
  {
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPSegment>::AddFaces(void* vtkFacesPointer) const
  {
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolygon>::AddFaces(void* vtkFacesPointer) const
  {
    if (vtkFacesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> faces = *(vtkSmartPointer<vtkCellArray>*)vtkFacesPointer;
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();

    AddFace(geometry.Face.row(0).transpose(),
            pointIds,
            faces);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolyhedron>::AddFaces(void* vtkFacesPointer) const
  {
    if (vtkFacesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> faces = *(vtkSmartPointer<vtkCellArray>*)vtkFacesPointer;
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();

    for (unsigned int i = 0 ; i < geometry.Faces.size(); i++)
    {
      AddFace(geometry.Faces[i].row(i),
              pointIds,
              faces);
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddSolution(const unsigned int& solutionPosition,
                                          void* vtkDoubleArrayPointer) const
  {
    if (solutionPosition >= numberSolutions)
      throw runtime_error("");

    const unsigned int& solutionSize = solutionSizes[solutionPosition];
    const double* solution = solutions[solutionPosition];

    if (solutionSize == 0 || solution == NULL || vtkDoubleArrayPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkDoubleArray& vtkSolution = *(vtkDoubleArray*)vtkDoubleArrayPointer;
    vtkSolution.SetNumberOfValues(solutionSize);
    for (unsigned int p = 0; p < solutionSize; p++)
      vtkSolution.SetValue(p, solution[p]);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
#if ENABLE_VTK == 1
  template <typename T>
  vtkSmartPointer<vtkPolyData> GeometryToPolyData<T>::Convert() const
  {
    vtkNew<vtkPolyData> polyData;

    std::vector<vtkSmartPointer<vtkDoubleArray>> vtkDoubleArrayPointers;
    vtkDoubleArrayPointers.reserve(numberSolutions);

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> vertices;
    vtkNew<vtkCellArray> lines;
    vtkNew<vtkCellArray> faces;

    InitializePoints(&points, &vertices);
    InitializeLines(&lines);
    InitializeFaces(&faces);

    polyData->SetPoints(points);
    polyData->SetVerts(vertices);
    polyData->SetLines(lines);
    polyData->SetPolys(faces);

    AddPoints(&points, &vertices);
    AddLines(&lines);
    AddFaces(&faces);

    if (numberSolutions > 0)
    {
      for (unsigned int s = 0; s < numberSolutions; s++)
      {
        vtkDoubleArrayPointers.push_back(vtkNew<vtkDoubleArray>());

        vtkDoubleArray* vtkSolution = (vtkDoubleArray*)vtkDoubleArrayPointers.back();

        const string& solutionLabel = solutionLabels[s];
        VTPProperty::Formats solutionFormat = solutionFormats[s];

        vtkSolution->SetName(solutionLabel.c_str());

        switch (solutionFormat)
        {
          case VTPProperty::Formats::Points:
            polyData->GetPointData()->AddArray(vtkSolution);

            if (s == 0)
              polyData->GetPointData()->SetActiveScalars(solutionLabel.c_str());
            break;
          default:
            throw runtime_error("Solution Format not supported");
        }

        AddSolution(s, vtkSolution);
      }
    }

    return polyData;
  }
#endif // ENABLE_VTK
  // ***************************************************************************
}
