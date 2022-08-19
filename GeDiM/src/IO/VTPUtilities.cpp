#include "VTPUtilities.hpp"

#include "IOUtilities.hpp"
#include "CommonUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  VTPUtilities::VTPUtilities()
  {
    _exportFormat = VTPUtilities::Binary;
  }
  VTPUtilities::~VTPUtilities()
  {
    properties.clear();
    geometries.clear();
  }
  // ***************************************************************************
  void VTPUtilities::AddProperty(const string& solutionLabel,
                                 const VTPProperty::Formats& format)
  {
    properties.push_back(VTPProperty());

    VTPProperty& property = properties.back();
    property.Label = solutionLabel;
    property.Format = format;
  }
  // ***************************************************************************

  // ***************************************************************************
  void VTPUtilities::AddGeometrySolution(const unsigned int& solutionSize,
                                         const double* solution)
  {
    VTPProperty& property = properties.back();

    property.Data.push_back(solution);
    property.Size.push_back(solutionSize);
  }
  // ***************************************************************************
  void VTPUtilities::Export(const string& filePath) const
  {
    Utilities::Unused(filePath); //< Is unused with VTK off

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkAppendFilter> append =	vtkSmartPointer<vtkAppendFilter>::New();

    vector<GeometryToPolyData> domainToPolyDatas;
    domainToPolyDatas.reserve(geometries.size());

    vector<vector<unsigned int>> geometryPropertySizes(properties.size());
    vector<vector<const double*>> geometryPropertyDatas(properties.size());

    unsigned int p = 0;
    for (const VTPProperty& property : properties)
    {
      geometryPropertySizes[p] = vector<unsigned int>(property.Size.begin(), property.Size.end());
      geometryPropertyDatas[p] = vector<const double*>(property.Data.begin(), property.Data.end());
    }

    unsigned int g = 0;
    for (const VTPGeometry& geometry : geometries)
    {
      domainToPolyDatas.push_back(GeometryToPolyData(geometry));

      GeometryToPolyData& domainToPolyData = domainToPolyDatas.back();

      domainToPolyData.InitializeSolutions(properties.size());

      unsigned int p = 0;
      for (const VTPProperty& property : properties)
      {
        domainToPolyData.SetSolutionOptions(p, property.Label, property.Format);
        domainToPolyData.AddSolution(p, geometryPropertySizes[p][g], geometryPropertyDatas[p][g]);
        p++;
      }

      domainToPolyData.Convert();

      vtkWeakPointer<vtkPolyData> polyData = (vtkPolyData*)domainToPolyData.PolyDataPointer();
      append->AddInputData(polyData);

      g++;
    }

    append->Update();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filePath.c_str());
    writer->SetInputData(append->GetOutput());

    switch (_exportFormat)
    {
      case VTPUtilities::Binary:
        writer->SetDataModeToBinary();
        break;
      case VTPUtilities::Ascii:
        writer->SetDataModeToAscii();
        break;
      case VTPUtilities::Appended:
        writer->SetDataModeToAppended();
        break;
      default:
        throw runtime_error("Export Format not supported");
        break;
    }

    writer->Write();
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  GeometryToPolyData::GeometryToPolyData(const VTPGeometry& geometry) :
    geometry(geometry)
  {
    numberSolutions = 0;

    vtkPointsPointer = NULL;
    vtkVerticesPointer = NULL;
    vtkLinesPointer = NULL;
    vtkFacesPointer = NULL;
    vtkPolyDataPointer = NULL;
  }
  GeometryToPolyData::~GeometryToPolyData()
  {
#if ENABLE_VTK == 1
    if (vtkPolyDataPointer != NULL)
    {
      vtkWeakPointer<vtkPolyData> polyData = (vtkPolyData*)vtkPolyDataPointer;
      polyData.GetPointer()->Delete();
      vtkPolyDataPointer = NULL;
    }

    if (vtkPointsPointer != NULL)
    {
      vtkWeakPointer<vtkPoints> points = (vtkPoints*)vtkPointsPointer;
      points.GetPointer()->Delete();
      vtkPointsPointer = NULL;
    }

    if (vtkVerticesPointer != NULL)
    {
      vtkWeakPointer<vtkCellArray> vertices = (vtkCellArray*)vtkVerticesPointer;
      vertices.GetPointer()->Delete();
      vtkVerticesPointer = NULL;
    }

    if (vtkLinesPointer != NULL)
    {
      vtkWeakPointer<vtkCellArray> lines = (vtkCellArray*)vtkLinesPointer;
      lines.GetPointer()->Delete();
      vtkLinesPointer = NULL;
    }

    if (vtkFacesPointer != NULL)
    {
      vtkWeakPointer<vtkCellArray> faces = (vtkCellArray*)vtkFacesPointer;
      faces.GetPointer()->Delete();
      vtkFacesPointer = NULL;
    }

    for (unsigned int i = 0; i < vtkDoubleArrayPointers.size(); i++)
    {
      void* vtkDoubleArrayPointer = vtkDoubleArrayPointers[i];
      if (vtkDoubleArrayPointer != NULL)
      {
        vtkWeakPointer<vtkDoubleArray> vtkSolution = (vtkDoubleArray*)vtkDoubleArrayPointer;
        vtkSolution.GetPointer()->Delete();
      }
    }
    vtkDoubleArrayPointers.clear();
#endif
  }
  // ***************************************************************************
  void GeometryToPolyData::InitializeSolutions(const unsigned int& _numberSolutions)
  {
    Output::Assert(_numberSolutions > 0);
    numberSolutions = _numberSolutions;

    solutionLabels.resize(numberSolutions, "Solution");
    solutionFormats.resize(numberSolutions, VTPProperty::Formats::Points);
    solutionSizes.resize(numberSolutions, 0);
    solutions.resize(numberSolutions, NULL);
    vtkDoubleArrayPointers.reserve(numberSolutions);
  }
  // ***************************************************************************
  void GeometryToPolyData::SetSolutionOptions(const unsigned int& solutionPosition,
                                              const string& solutionLabel,
                                              const VTPProperty::Formats& format)
  {
    Output::Assert(solutionPosition < numberSolutions);

    solutionLabels[solutionPosition] = solutionLabel;
    solutionFormats[solutionPosition] = format;
  }
  // ***************************************************************************
  void GeometryToPolyData::AddSolution(const unsigned int& solutionPosition,
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
  void GeometryToPolyData::AddPoint(const Eigen::Vector3d& point,
                                    vtkSmartPointer<vtkPoints>& points,
                                    vtkSmartPointer<vtkCellArray>& vertices) const
  {
    vtkIdType pointId[1];
    pointId[0] = points->InsertNextPoint(point(0), point(1), point(2));

    vertices->InsertNextCell(1, pointId);
  }
  // ***************************************************************************
  void GeometryToPolyData::AddLine(const unsigned int& originId,
                                   const unsigned int& endId,
                                   vtkSmartPointer<vtkIdList>& pointIds,
                                   vtkSmartPointer<vtkCellArray>& lines) const
  {
    pointIds->Initialize();
    pointIds->Allocate(2);

    pointIds->InsertNextId(originId);
    pointIds->InsertNextId(endId);

    pointIds->Squeeze();

    lines->InsertNextCell(pointIds);

    pointIds->Reset();
  }
  // ***************************************************************************
  void GeometryToPolyData::AddFace(const Eigen::VectorXi& faceVerticesIds,
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
  void GeometryToPolyData::InitializePoints(void* vtkPointsPointer,
                                            void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    unsigned int numberTotalPoints = 0;

    switch (geometry.Type)
    {
      case VTPGeometry::Point:
        numberTotalPoints += 1;
        break;
      case VTPGeometry::Segment:
        numberTotalPoints += 2;
        break;
      case VTPGeometry::Polygon:
        numberTotalPoints += static_cast<const VTPPolygon&>(geometry).Vertices.cols();
        break;
      case VTPGeometry::Polyhedron:
        numberTotalPoints += static_cast<const VTPPolyhedron&>(geometry).Vertices.cols();
        break;
      default:
        throw runtime_error("");
    }

    if (numberTotalPoints > 0)
    {
      points->Allocate(numberTotalPoints);
      vertices->Allocate(numberTotalPoints);
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  void GeometryToPolyData::AddPoints(void* vtkPointsPointer,
                                     void* vtkVerticesPointer) const
  {
    if (vtkPointsPointer == NULL || vtkVerticesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkPoints>& points = *(vtkSmartPointer<vtkPoints>*)vtkPointsPointer;
    vtkSmartPointer<vtkCellArray>& vertices = *(vtkSmartPointer<vtkCellArray>*)vtkVerticesPointer;

    switch (geometry.Type)
    {
      case VTPGeometry::Point:
      {
        AddPoint(static_cast<const VTPPoint&>(geometry).Point,
                 points,
                 vertices);
      }
        break;
      case VTPGeometry::Segment:
      {
        const VTPSegment& segment = static_cast<const VTPSegment&>(geometry);

        AddPoint(segment.Vertices.col(0),
                 points,
                 vertices);
        AddPoint(segment.Vertices.col(1),
                 points,
                 vertices);
      }
        break;
      case VTPGeometry::Polygon:
      {
        const VTPPolygon& polygon = static_cast<const VTPPolygon&>(geometry);
        for (unsigned int i = 0 ; i < polygon.Vertices.cols(); i++)
        {
          AddPoint(polygon.Vertices.col(i),
                   points,
                   vertices);
        }
      }
        break;
      case VTPGeometry::Polyhedron:
      {
        const VTPPolyhedron& polyhedron = static_cast<const VTPPolyhedron&>(geometry);
        for (unsigned int i = 0 ; i < polyhedron.Vertices.cols(); i++)
        {
          AddPoint(polyhedron.Vertices.col(i),
                   points,
                   vertices);
        }
      }
        break;
      default:
        throw runtime_error("");
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  void GeometryToPolyData::InitializeLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;

    unsigned int numberTotalLines = 0;

    switch (geometry.Type)
    {
      case VTPGeometry::Point:
        break;
      case VTPGeometry::Segment:
        numberTotalLines += 1;
        break;
      case VTPGeometry::Polygon:
        numberTotalLines += static_cast<const VTPPolygon&>(geometry).Vertices.cols();
        break;
      case VTPGeometry::Polyhedron:
        numberTotalLines += static_cast<const VTPPolyhedron&>(geometry).Edges.cols();
        break;
      default:
        throw runtime_error("");
    }

    if (numberTotalLines > 0)
      lines->Allocate(numberTotalLines);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  void GeometryToPolyData::AddLines(void* vtkLinesPointer) const
  {
    if (vtkLinesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> lines = *(vtkSmartPointer<vtkCellArray>*)vtkLinesPointer;
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();

    switch (geometry.Type)
    {
      case VTPGeometry::Point:
        break;
      case VTPGeometry::Segment:
      {
        const VTPSegment& segment = static_cast<const VTPSegment&>(geometry);
        AddLine(segment.Edge[0],
            segment.Edge[1],
            pointIds,
            lines);
      }
        break;
      case VTPGeometry::Polygon:
      {
        const VTPPolygon& polygon = static_cast<const VTPPolygon&>(geometry);
        for (unsigned int i = 0 ; i < polygon.Edges.cols(); i++)
        {
          AddLine(polygon.Edges(0, i),
                  polygon.Edges(1, i),
                  pointIds,
                  lines);
        }
      }
        break;
      case VTPGeometry::Polyhedron:
      {
        const VTPPolyhedron& polyhedron = static_cast<const VTPPolyhedron&>(geometry);
        for (unsigned int i = 0 ; i < polyhedron.Edges.cols(); i++)
        {
          AddLine(polyhedron.Edges(0, i),
                  polyhedron.Edges(1, i),
                  pointIds,
                  lines);
        }
      }
        break;
      default:
        throw runtime_error("");
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  void GeometryToPolyData::InitializeFaces(void* vtkFacesPointer) const
  {
    if (vtkFacesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> faces = *(vtkSmartPointer<vtkCellArray>*)vtkFacesPointer;

    unsigned int numberTotalFaces = 0;

    switch (geometry.Type)
    {
      case VTPGeometry::Point:
        break;
      case VTPGeometry::Segment:
        break;
      case VTPGeometry::Polygon:
        numberTotalFaces += 1;
        break;
      case VTPGeometry::Polyhedron:
        numberTotalFaces += static_cast<const VTPPolyhedron&>(geometry).Faces.size();
        break;
      default:
        throw runtime_error("");
    }

    if (numberTotalFaces > 0)
      faces->Allocate(numberTotalFaces);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  void GeometryToPolyData::AddFaces(void* vtkFacesPointer) const
  {
    if (vtkFacesPointer == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSmartPointer<vtkCellArray> faces = *(vtkSmartPointer<vtkCellArray>*)vtkFacesPointer;
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();

    switch (geometry.Type)
    {
      case VTPGeometry::Point:
        break;
      case VTPGeometry::Segment:
        break;
      case VTPGeometry::Polygon:
      {
        const VTPPolygon& polygon = static_cast<const VTPPolygon&>(geometry);
        Eigen::VectorXi faceVerticesId = Eigen::VectorXi::LinSpaced(polygon.Vertices.cols(),
                                                                    0.0,
                                                                    polygon.Vertices.cols() - 1);
        AddFace(faceVerticesId,
                pointIds,
                faces);
      }
        break;
      case VTPGeometry::Polyhedron:
      {
        const VTPPolyhedron& polyhedron = static_cast<const VTPPolyhedron&>(geometry);
        for (unsigned int i = 0 ; i < polyhedron.Faces.size(); i++)
        {
          AddFace(polyhedron.Faces[i].row(i),
                  pointIds,
                  faces);
        }
      }
        break;
      default:
        throw runtime_error("");
    }
#endif // ENABLE_VTK
  }
  // ***************************************************************************
  void GeometryToPolyData::AddSolution(const unsigned int& solutionPosition,
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
  void GeometryToPolyData::Convert()
  {
#if ENABLE_VTK == 1
    vtkPolyDataPointer = vtkPolyData::New();
    vtkPointsPointer = vtkPoints::New();
    vtkVerticesPointer = vtkCellArray::New();
    vtkLinesPointer = vtkCellArray::New();
    vtkFacesPointer = vtkCellArray::New();

    vtkWeakPointer<vtkPolyData> polyData = (vtkPolyData*)vtkPolyDataPointer;
    vtkWeakPointer<vtkPoints> points = (vtkPoints*)vtkPointsPointer;
    vtkWeakPointer<vtkCellArray> vertices = (vtkCellArray*)vtkVerticesPointer;
    vtkWeakPointer<vtkCellArray> lines = (vtkCellArray*)vtkLinesPointer;
    vtkWeakPointer<vtkCellArray> faces = (vtkCellArray*)vtkFacesPointer;

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
        vtkDoubleArrayPointers.push_back(vtkDoubleArray::New());

        vtkDoubleArray* vtkSolution = (vtkDoubleArray*)vtkDoubleArrayPointers.back();

        string& solutionLabel = solutionLabels[s];
        VTPProperty::Formats solutionFormat = solutionFormats[s];

        vtkSolution->SetName(solutionLabel.c_str());

        switch (solutionFormat) {
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
#endif // ENABLE_VTK
  }
  // ***************************************************************************
}
