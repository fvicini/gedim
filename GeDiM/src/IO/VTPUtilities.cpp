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
#if ENABLE_VTK == 1
    exportData = vtkSmartPointer<vtkAppendFilter>::New();
#endif
  }
  VTKUtilities::~VTKUtilities()
  {
#if ENABLE_VTK == 1
    exportData->RemoveAllInputs();
    exportData->Delete();
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddPoint(const Eigen::Vector3d& point)
  {
    VTPPoint vtpPoint(point);
    GeometryToPolyData<VTPPoint> polyData(vtpPoint);

#if ENABLE_VTK == 1
    exportData->AddInputData((vtkPolyData*)polyData.Convert());
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
  //  // ***************************************************************************
  //  VTPUtilities::VTPUtilities()
  //  {
  //    _exportFormat = VTPUtilities::Binary;
  //  }
  //  VTPUtilities::~VTPUtilities()
  //  {
  //    for (IVTPGeometry* geometry : geometries)
  //      delete geometry;
  //    geometries.clear();

  //    properties.clear();
  //  }
  //  // ***************************************************************************
  //  void VTPUtilities::AddProperty(const string& solutionLabel,
  //                                 const VTPProperty::Formats& format)
  //  {
  //    properties.insert(make_pair(solutionLabel, VTPProperty()));

  //    VTPProperty& property = properties.at(solutionLabel);
  //    property.Label = solutionLabel;
  //    property.Format = format;
  //  }
  //  // ***************************************************************************
  //  void VTPUtilities::AddGeometryProperty(const string& solutionLabel,
  //                                         const unsigned int& solutionSize,
  //                                         const double* solution)
  //  {
  //    VTPProperty& property = properties.at(solutionLabel);

  //    property.Data.push_back(solution);
  //    property.Size.push_back(solutionSize);
  //  }
  //  // ***************************************************************************
  //  void VTPUtilities::Export(const string& filePath) const
  //  {
  //    Utilities::Unused(filePath); //< Is unused with VTK off

  //#if ENABLE_VTK == 1
  //    vtkSmartPointer<vtkAppendFilter> append =	vtkSmartPointer<vtkAppendFilter>::New();

  //    vector<IGeometryToPolyData*> domainToPolyDatas;
  //    domainToPolyDatas.reserve(geometries.size());

  //    vector<vector<unsigned int>> geometryPropertySizes(properties.size());
  //    vector<vector<const double*>> geometryPropertyDatas(properties.size());

  //    unsigned int p = 0;
  //    for (const pair<string, VTPProperty>& propertyPair : properties)
  //    {
  //      const VTPProperty& property = propertyPair.second;
  //      geometryPropertySizes[p] = vector<unsigned int>(property.Size.begin(), property.Size.end());
  //      geometryPropertyDatas[p] = vector<const double*>(property.Data.begin(), property.Data.end());
  //      p++;
  //    }

  //    unsigned int g = 0;
  //    for (const IVTPGeometry* geometry : geometries)
  //    {
  //      switch (geometry->Type())
  //      {
  //        case IVTPGeometry::Types::Point:
  //          domainToPolyDatas.push_back(new GeometryToPolyData<VTPPoint>(*static_cast<const VTPPoint*>(geometry)));
  //          break;
  //        case IVTPGeometry::Types::Segment:
  //          domainToPolyDatas.push_back(new GeometryToPolyData<VTPSegment>(*static_cast<const VTPSegment*>(geometry)));
  //          break;
  //        case IVTPGeometry::Types::Polygon:
  //          domainToPolyDatas.push_back(new GeometryToPolyData<VTPPolygon>(*static_cast<const VTPPolygon*>(geometry)));
  //          break;
  //        case IVTPGeometry::Types::Polyhedron:
  //          domainToPolyDatas.push_back(new GeometryToPolyData<VTPPolyhedron>(*static_cast<const VTPPolyhedron*>(geometry)));
  //          break;
  //        default:
  //          throw runtime_error("");
  //      }

  //      IGeometryToPolyData& domainToPolyData = *domainToPolyDatas.back();

  //      domainToPolyData.InitializeSolutions(properties.size());

  //      unsigned int p = 0;
  //      for (const pair<string, VTPProperty>& propertyPair : properties)
  //      {
  //        const VTPProperty& property = propertyPair.second;
  //        domainToPolyData.SetSolutionOptions(p,
  //                                            property.Label,
  //                                            property.Format);
  //        domainToPolyData.AddSolution(p,
  //                                     geometryPropertySizes[p][g],
  //                                     geometryPropertyDatas[p][g]);
  //        p++;
  //      }

  //      domainToPolyData.Convert();

  //      vtkWeakPointer<vtkPolyData> polyData = (vtkPolyData*)domainToPolyData.PolyDataPointer();
  //      append->AddInputData(polyData);

  //      g++;
  //    }

  //    append->Update();

  //    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  //    writer->SetFileName(filePath.c_str());
  //    writer->SetInputData(append->GetOutput());

  //    switch (_exportFormat)
  //    {
  //      case VTPUtilities::Binary:
  //        writer->SetDataModeToBinary();
  //        break;
  //      case VTPUtilities::Ascii:
  //        writer->SetDataModeToAscii();
  //        break;
  //      case VTPUtilities::Appended:
  //        writer->SetDataModeToAppended();
  //        break;
  //      default:
  //        throw runtime_error("Export Format not supported");
  //        break;
  //    }

  //    writer->Write();
  //#endif // ENABLE_VTK
  //  }
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
    //#if ENABLE_VTK == 1
    //    if (vtkPolyDataPointer != NULL)
    //    {
    //      vtkWeakPointer<vtkPolyData> polyData = (vtkPolyData*)vtkPolyDataPointer;
    //      polyData.GetPointer()->Delete();
    //      vtkPolyDataPointer = NULL;
    //    }

    //    if (vtkPointsPointer != NULL)
    //    {
    //      vtkWeakPointer<vtkPoints> points = (vtkPoints*)vtkPointsPointer;
    //      points.GetPointer()->Delete();
    //      vtkPointsPointer = NULL;
    //    }

    //    if (vtkVerticesPointer != NULL)
    //    {
    //      vtkWeakPointer<vtkCellArray> vertices = (vtkCellArray*)vtkVerticesPointer;
    //      vertices.GetPointer()->Delete();
    //      vtkVerticesPointer = NULL;
    //    }

    //    if (vtkLinesPointer != NULL)
    //    {
    //      vtkWeakPointer<vtkCellArray> lines = (vtkCellArray*)vtkLinesPointer;
    //      lines.GetPointer()->Delete();
    //      vtkLinesPointer = NULL;
    //    }

    //    if (vtkFacesPointer != NULL)
    //    {
    //      vtkWeakPointer<vtkCellArray> faces = (vtkCellArray*)vtkFacesPointer;
    //      faces.GetPointer()->Delete();
    //      vtkFacesPointer = NULL;
    //    }

    //    for (unsigned int i = 0; i < vtkDoubleArrayPointers.size(); i++)
    //    {
    //      void* vtkDoubleArrayPointer = vtkDoubleArrayPointers[i];
    //      if (vtkDoubleArrayPointer != NULL)
    //      {
    //        vtkWeakPointer<vtkDoubleArray> vtkSolution = (vtkDoubleArray*)vtkDoubleArrayPointer;
    //        vtkSolution.GetPointer()->Delete();
    //      }
    //    }
    //    vtkDoubleArrayPointers.clear();
    //#endif
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
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();

    AddLine(geometry.Edge[0],
        geometry.Edge[1],
        pointIds,
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
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();

    for (unsigned int i = 0 ; i < geometry.Edges.cols(); i++)
    {
      AddLine(geometry.Edges(0, i),
              geometry.Edges(1, i),
              pointIds,
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
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();

    for (unsigned int i = 0 ; i < geometry.Edges.cols(); i++)
    {
      AddLine(geometry.Edges(0, i),
              geometry.Edges(1, i),
              pointIds,
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

    Eigen::VectorXi faceVerticesId = Eigen::VectorXi::LinSpaced(geometry.Vertices.cols(),
                                                                0.0,
                                                                geometry.Vertices.cols() - 1);
    AddFace(faceVerticesId,
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
  template <typename T>
  void* GeometryToPolyData<T>::Convert() const
  {
    void* vtkPolyDataPointer = nullptr;

#if ENABLE_VTK == 1
    vtkPolyDataPointer = vtkPolyData::New();
    void* vtkPointsPointer = vtkPoints::New();
    void* vtkVerticesPointer = vtkCellArray::New();
    void* vtkLinesPointer = vtkCellArray::New();
    void* vtkFacesPointer = vtkCellArray::New();
    std::vector<void*> vtkDoubleArrayPointers;
    vtkDoubleArrayPointers.reserve(numberSolutions);

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
#endif // ENABLE_VTK

    return vtkPolyDataPointer;
  }
  // ***************************************************************************
}
