#include "VoroInterface.hpp"

#include "Macro.hpp"
#include "Eigen/Eigen"
#include "CommonUtilities.hpp"


namespace Gedim
{
// ************************************************************************* //
VoroInterface::VoroInterface(const Gedim::GeometryUtilities& geometryUtilities):
    geometryUtilities(geometryUtilities)
{

}
// ************************************************************************* //
#if ENABLE_VORO == 0
void VoroInterface::GenerateVoronoiTassellations(const Eigen::MatrixXd& polyhedronVertices,
                                                 const Eigen::MatrixXi& polyhedronEdges,
                                                 const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                 const unsigned int& numPoints,
                                                 Gedim::IMeshDAO& mesh)
{
    Gedim::Utilities::Unused(polyhedronVertices);
    Gedim::Utilities::Unused(polyhedronEdges);
    Gedim::Utilities::Unused(polyhedronFaces);
    Gedim::Utilities::Unused(numPoints);
    Gedim::Utilities::Unused(mesh);
    throw runtime_error("Not active module VORO");
}
#else
// ************************************************************************* //
bool VoroInterface::InsertNewPoints(Cell0D& cell0D,
                                    list<Cell0D>& cell0Ds)
{
    if(cell0Ds.size() == 0)
    {
        cell0D.id = 0;
        cell0Ds.push_back(cell0D);
        return true;
    }
    else{
        for(auto it = cell0Ds.begin(); it != cell0Ds.end(); it++)
        {
            if( geometryUtilities.IsValue1DGreater((*it).x, cell0D.x))
            {
                cell0D.id = cell0Ds.size();
                cell0Ds.insert(it, cell0D);
                return true;
            }
            else if(geometryUtilities.IsValue1DZero((*it).x - cell0D.x))
            {
                if( geometryUtilities.IsValue1DGreater((*it).y, cell0D.y))
                {
                    cell0D.id = cell0Ds.size();
                    cell0Ds.insert(it, cell0D);
                    return true;
                }
                else if(geometryUtilities.IsValue1DZero((*it).y - cell0D.y))
                {
                    if( geometryUtilities.IsValue1DGreater((*it).z, cell0D.z))
                    {
                        cell0D.id = cell0Ds.size();
                        cell0Ds.insert(it, cell0D);
                        return true;
                    }
                    else if(geometryUtilities.IsValue1DZero((*it).z - cell0D.z))
                    {
                        cell0D.id = (*it).id;
                        return false;
                    }
                }
            }
        }

        cell0D.id = cell0Ds.size();
        cell0Ds.push_back(cell0D);
        return true;
    }
}
// ************************************************************************* //
void VoroInterface::GenerateVoronoiTassellations(const Eigen::MatrixXd& polyhedronVertices,
                                                 const Eigen::MatrixXi& polyhedronEdges,
                                                 const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                 const unsigned int& numPoints,
                                                 Gedim::IMeshDAO& mesh)
{
    const unsigned int numVerticesDomain = polyhedronVertices.cols();
    const unsigned int numEdgesDomain = polyhedronEdges.cols();
    const unsigned int numFacesDomain = polyhedronFaces.size();


    vector<Eigen::Vector3d> edgeDomainTangent(numEdgesDomain);
    vector<double> edgeDomainTangentSquaredNorm(numEdgesDomain);
    vector<Eigen::Vector3d> edgeDomainLineOrigin(numEdgesDomain);
    for(unsigned int e = 0; e < numEdgesDomain; e++)
    {
        edgeDomainLineOrigin[e] = polyhedronVertices.col(polyhedronEdges(0, e));
        const Eigen::Vector3d end = polyhedronVertices.col(polyhedronEdges(1, e));

        Eigen::Vector3d tangent = end - edgeDomainLineOrigin[e];

        edgeDomainTangent[e] = tangent;
        edgeDomainTangentSquaredNorm[e] = tangent.squaredNorm();
    }


    vector<Eigen::Vector3d> faceDomainNormal(numFacesDomain);
    for(unsigned int f = 0; f < numFacesDomain; f++)
    {
        const unsigned int numVerticesFaceDomain = polyhedronFaces[f].cols();
        Eigen::MatrixXd vert = Eigen::MatrixXd::Zero(3, numVerticesFaceDomain);

        for(unsigned int v = 0; v < numVerticesFaceDomain; v++)
            vert.col(v) = polyhedronVertices.col(polyhedronFaces[f](0,v));

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(vert);

        faceDomainNormal[f] = normal;
    }

    vector<unsigned int> map_marker(numFacesDomain, 0);

    const double x_min = polyhedronVertices.row(0).minCoeff(), x_max = polyhedronVertices.row(0).maxCoeff();
    const double y_min = polyhedronVertices.row(1).minCoeff(), y_max = polyhedronVertices.row(1).maxCoeff();
    const double z_min = polyhedronVertices.row(2).minCoeff(), z_max = polyhedronVertices.row(2).maxCoeff();
    const double cvol = (x_max - x_min) * (y_max - y_min) * (x_max - x_min);
    double vvol = 0.0;

    // Set up the number of blocks that the container is divided into
    const int n_x = 6, n_y = 6, n_z = 6;

    // Initialize Voronoi mesh
    voro::container con
        (
            x_min, x_max, y_min, y_max, z_min, z_max,
            n_x, n_y, n_z,
            false, false, false, 8
            );

    GenerateRandomPoints(polyhedronVertices,
                         numPoints,
                         con);

    // Get cell number
    unsigned int numberOfCells3D = con.total_particles();
    unsigned int numberOfCellsMesh = 0;

    vector<vector<unsigned int>> cell3DVertices(numberOfCells3D);


    vector<Cell2D> cell2Ds;
    vector<vector<unsigned int>> cellEdges(numberOfCells3D);
    vector<vector<unsigned int>> cellFaces(numberOfCells3D);
    vector<vector<int>> cellFaceNeighbour(numberOfCells3D);
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> cell1Ds;


    list<Cell0D> cell0Ds;

    // Loop on Voronoi cells
    voro::c_loop_all vl(con);
    int ijk, q; double *pp;
    voro::voronoicell_neighbor c;
    if (vl.start())
    {
        do{
            if (con.compute_cell(c, vl))
            {
                ijk = vl.ijk;
                q = vl.q;
                pp = con.p[ijk] + con.ps*q;

                // Cell barycenter
                const double cx = *pp;
                const double cy = pp[1];
                const double cz = pp[2];

                // Get cell3D ID
                const int ID = con.id[ijk][q];
                numberOfCellsMesh++;

                // Get cell vertices
                std::vector<double> vertices;
                c.vertices(cx, cy, cz, vertices);

                // Get cell face cardinality
                std::vector<int> faceOrders;
                c.face_orders(faceOrders);

                // Get face-vertices connectivity
                std::vector<int> localFaceVertices;
                c.face_vertices(localFaceVertices);

                // Get cell neighbors IDs
                std::vector<int> faceNeighs;
                c.neighbors(faceNeighs);


                // conservative resize
                unsigned int ID_max = numberOfCells3D - 1;
                for(unsigned int id_max = 0; id_max < faceNeighs.size(); id_max++)
                    if(faceNeighs[id_max] >= 0)
                        if(static_cast<unsigned int>(faceNeighs[id_max]) > ID_max)
                            ID_max = faceNeighs[id_max];
                ID_max = max(ID_max, static_cast<unsigned int>(ID));
                if(ID_max > numberOfCells3D - 1)
                {
                    numberOfCells3D = ID_max + 1;
                    cell3DVertices.resize(numberOfCells3D);
                    cellEdges.resize(numberOfCells3D);
                    cellFaces.resize(numberOfCells3D);
                    cellFaceNeighbour.resize(numberOfCells3D);
                }

                // Fill mesh connectivity structures
                // cell vertices
                const unsigned int nCoords = vertices.size();
                const unsigned int numVertices = nCoords/3;
                cell3DVertices[ID].resize(numVertices);

                Eigen::MatrixXd cellVertCoordinates = Eigen::MatrixXd::Zero(3, numVertices);

                unsigned int v = 0;
                unsigned int i = 0;
                while (i < nCoords)
                {
                    double x = vertices[i++];
                    double y = vertices[i++];
                    double z = vertices[i++];

                    Cell0D cell0D(x,y,z);
                    InsertNewPoints(cell0D, cell0Ds);

                    cell3DVertices[ID][v] = cell0D.id;
                    cellVertCoordinates.col(v) << x, y, z;
                    v++;
                }

                // cell faces
                unsigned int nFaces = faceOrders.size();
                cellFaceNeighbour[ID].resize(nFaces);
                cellFaces[ID].resize(nFaces);

                int count = 0;
                for (unsigned int i = 0; i < nFaces; i++)
                {
                    const unsigned numFaceVertices = faceOrders[i];
                    vector<unsigned int> faceVertices(numFaceVertices);
                    Eigen::MatrixXd faceVertCoordinates = Eigen::MatrixXd::Zero(3, numFaceVertices);
                    cellFaceNeighbour[ID][i] = faceNeighs[i];

                    count++;
                    for (unsigned int j = 0; j < numFaceVertices; j++)
                    {
                        faceVertices[j] = cell3DVertices[ID][localFaceVertices[count]];
                        faceVertCoordinates.col(j) = cellVertCoordinates.col(localFaceVertices[count]);
                        count++;
                    }

                    if(cellFaceNeighbour[ID][i] < 0) // if it is a boundary face
                    {
                        // prima volta che incontro la faccia  -  nuovi lati possibili - faccia di bordo

                        Cell2D cell2D;
                        cell2D.id = cell2Ds.size();;
                        cell2D.vertices = faceVertices;

                        if( map_marker[abs(faceNeighs[i]) - 1] == 0)
                        {
                            const Eigen::Vector3d barycenter = faceVertCoordinates.rowwise().mean();

                            for(unsigned int f = 0; f < numFacesDomain; f++)
                            {

                                if(geometryUtilities.IsPointOnPlane(barycenter,
                                                                     faceDomainNormal[f],
                                                                     polyhedronVertices.col(polyhedronFaces[f](0,0))))
                                {
                                    map_marker[abs(faceNeighs[i]) - 1] = numVerticesDomain + numEdgesDomain + f + 1;
                                    cell2D.marker = map_marker[abs(faceNeighs[i]) - 1];
                                    break;
                                }

                            }

                            if( map_marker[abs(faceNeighs[i]) - 1] == 0)
                                throw runtime_error("Marker face wrong");

                        }
                        else
                            cell2D.marker = map_marker[abs(faceNeighs[i]) - 1];

                        cell2D.edges.resize(numFaceVertices);
                        for(unsigned int v = 0; v < numFaceVertices; v++)
                        {
                            const unsigned int id_min = min(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);
                            const unsigned int id_max = max(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);

                            auto it = cell1Ds.find({id_min,id_max});

                            if( it != cell1Ds.end())
                            {
                                cell2D.edges[v] = it -> second[0];

                                if((it -> second[1]) == 0)
                                {
                                    for(unsigned int e = 0; e < numEdgesDomain; e++)
                                    {

                                        if(geometryUtilities.IsPointOnLine(faceVertCoordinates.col(v),
                                                                            edgeDomainLineOrigin[e],
                                                                            edgeDomainTangent[e],
                                                                            edgeDomainTangentSquaredNorm[e]) &&
                                            geometryUtilities.IsPointOnLine(faceVertCoordinates.col((v+1)%numFaceVertices),
                                                                            edgeDomainLineOrigin[e],
                                                                            edgeDomainTangent[e],
                                                                            edgeDomainTangentSquaredNorm[e]))
                                        {
                                            (it -> second[1]) = numVerticesDomain + e + 1;
                                            break;
                                        }

                                        if((it -> second[1]) == 0)
                                            (it -> second[1]) = map_marker[abs(faceNeighs[i]) - 1];

                                    }
                                }

                            }
                            else
                            {
                                unsigned int marker = 0;

                                for(unsigned int e = 0; e < numEdgesDomain; e++)
                                {

                                    if(geometryUtilities.IsPointOnLine(faceVertCoordinates.col(v),
                                                                        edgeDomainLineOrigin[e],
                                                                        edgeDomainTangent[e],
                                                                        edgeDomainTangentSquaredNorm[e]) &&
                                        geometryUtilities.IsPointOnLine(faceVertCoordinates.col((v+1)%numFaceVertices),
                                                                        edgeDomainLineOrigin[e],
                                                                        edgeDomainTangent[e],
                                                                        edgeDomainTangentSquaredNorm[e]))
                                    {
                                        marker = numVerticesDomain + e + 1;
                                        break;
                                    }

                                    if(marker == 0)
                                        marker = map_marker[abs(faceNeighs[i]) - 1];

                                }
                                vector<unsigned int> cell1D(2);
                                cell1D[0] = cell1Ds.size();
                                cell1D[1] = marker;
                                cell1Ds.insert({{id_min, id_max}, cell1D});

                                cell2D.edges[v] = cell1Ds.size() - 1;
                            }
                        }

                        cell2Ds.push_back(cell2D);

                        cellFaces[ID][i] = cell2D.id;
                        for(unsigned int e = 0; e < cell2D.edges.size(); e++)
                            if(find(cellEdges[ID].begin(), cellEdges[ID].end(), cell2D.edges[e]) == cellEdges[ID].end())
                                cellEdges[ID].push_back(cell2D.edges[e]);

                    }
                    else
                    {
                        if(cellFaceNeighbour[faceNeighs[i]].size()==0)
                        {
                            // prima volta che incontro la faccia -  nuovi lati possibili

                            Cell2D cell2D;
                            cell2D.id = cell2Ds.size();
                            cell2D.marker = 0;
                            cell2D.vertices = faceVertices;

                            cell2D.edges.resize(numFaceVertices);
                            for(unsigned int v = 0; v < numFaceVertices; v++)
                            {
                                const unsigned int id_min = min(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);
                                const unsigned int id_max = max(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);

                                auto it = cell1Ds.find({id_min,id_max});

                                if( it != cell1Ds.end())
                                    cell2D.edges[v] = it -> second[0];
                                else{
                                    vector<unsigned int> cell1D(2);
                                    cell1D[0] = cell1Ds.size();
                                    cell1D[1] = 0;
                                    cell1Ds.insert({{id_min, id_max}, cell1D});
                                    cell2D.edges[v] = cell1Ds.size() - 1;
                                }
                            }

                            cell2Ds.push_back(cell2D);

                            cellFaces[ID][i] = cell2D.id;
                            for(unsigned int e = 0; e < cell2D.edges.size(); e++)
                                if(find(cellEdges[ID].begin(), cellEdges[ID].end(), cell2D.edges[e]) == cellEdges[ID].end())
                                    cellEdges[ID].push_back(cell2D.edges[e]);

                        }
                        else{
                            for(unsigned int f = 0; f < cellFaceNeighbour[faceNeighs[i]].size(); f++)
                            {

                                // faccia già esistente -  nessun nuovo lato

                                if(cellFaceNeighbour[faceNeighs[i]][f] == ID)
                                {
                                    cellFaces[ID][i] = cellFaces[faceNeighs[i]][f];

                                    for(unsigned int e = 0; e < cell2Ds[cellFaces[ID][i]].edges.size(); e++)
                                        if(find(cellEdges[ID].begin(), cellEdges[ID].end(), cell2Ds[cellFaces[ID][i]].edges[e]) == cellEdges[ID].end())
                                            cellEdges[ID].push_back(cell2Ds[cellFaces[ID][i]].edges[e]);

                                    break;
                                }

                                if(f == cellFaceNeighbour[faceNeighs[i]].size() - 1)
                                    throw runtime_error("ID face not found");
                            }

                        }
                    }
                }
                vvol += c.volume();

            }
        }
        while(vl.inc());
    }

    //    if(!geometryUtilities.IsValue3DZero(cvol - vvol))
    //        throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");

    if(abs(cvol - vvol) > 1.0e-12)
        throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");

    /// <li> Set Cell0Ds
    const unsigned int numberOfPointsMesh = cell0Ds.size();
    const unsigned int numberOfEdgesMesh = cell1Ds.size();
    const unsigned int numberOfFacesMesh = cell2Ds.size();


    mesh.InitializeDimension(3);
    mesh.Cell0DsInitialize(numberOfPointsMesh);
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    mesh.Cell2DsInitialize(numberOfFacesMesh);
    mesh.Cell3DsInitialize(numberOfCellsMesh);

    for (const Cell0D& cell0D : cell0Ds)
    {
        const unsigned int id = cell0D.id;

        mesh.Cell0DSetState(id, true);
        mesh.Cell0DInsertCoordinates(id,
                                     Eigen::Vector3d(cell0D.x,
                                                     cell0D.y,
                                                     cell0D.z));


        const bool is_xmin = (geometryUtilities.IsValue1DZero(cell0D.x - x_min)); // l == 0
        const bool is_xmax = (geometryUtilities.IsValue1DZero(cell0D.x - x_max)); // l == (numLengthPoints - 1)
        const bool is_ymin = (geometryUtilities.IsValue1DZero(cell0D.y - y_min)); // h == 0
        const bool is_ymax = (geometryUtilities.IsValue1DZero(cell0D.y - y_max)); // h == (numHeightPoints - 1)
        const bool is_zmin = (geometryUtilities.IsValue1DZero(cell0D.z - z_min)); // w == 0
        const bool is_zmax = (geometryUtilities.IsValue1DZero(cell0D.z - z_max)); // w == (numWidthPoints - 1)

        const unsigned int marker = 1 * (is_xmin && is_ymin && is_zmin) +
                                    2 * (is_xmax && is_ymin && is_zmin) +
                                    3 * (is_xmax && is_ymax && is_zmin) +
                                    4 * (is_xmin && is_ymax && is_zmin) +
                                    5 * (is_xmin && is_ymin && is_zmax) +
                                    6 * (is_xmax && is_ymin && is_zmax) +
                                    7 * (is_xmax && is_ymax && is_zmax) +
                                    8 * (is_xmin && is_ymax && is_zmax) +
                                    9 * (!is_xmin && !is_xmax && is_ymin && is_zmin) +
                                    10 * (is_xmax && !is_ymin && !is_ymax && is_zmin) +
                                    11 * (!is_xmin && !is_xmax && is_ymax && is_zmin) +
                                    12 * (is_xmin && !is_ymin && !is_ymax && is_zmin) +
                                    13 * (!is_xmin && !is_xmax && is_ymin && is_zmax) +
                                    14 * (is_xmax && !is_ymin && !is_ymax && is_zmax) +
                                    15 * (!is_xmin && !is_xmax && is_ymax && is_zmax) +
                                    16 * (is_xmin && !is_ymin && !is_ymax && is_zmax) +
                                    17 * (is_xmin && is_ymin && !is_zmin && !is_zmax) +
                                    18 * (is_xmax && is_ymin && !is_zmin && !is_zmax) +
                                    19 * (is_xmax && is_ymax && !is_zmin && !is_zmax) +
                                    20 * (is_xmin && is_ymax && !is_zmin && !is_zmax) +
                                    21 * (!is_xmin && !is_xmax && !is_ymin && !is_ymax && is_zmin) +
                                    22 * (!is_xmin && !is_xmax && !is_ymin && !is_ymax && is_zmax) +
                                    23 * (is_xmin && !is_ymin && !is_ymax && !is_zmin && !is_zmax) +
                                    24 * (is_xmax && !is_ymin && !is_ymax && !is_zmin && !is_zmax) +
                                    25 * (!is_xmin && !is_xmax && is_ymin  && !is_zmin && !is_zmax) +
                                    26 * (!is_xmin && !is_xmax &&  is_ymax && !is_zmin && !is_zmax);

        mesh.Cell0DSetMarker(id, marker);
    }

    /// <li> Set Edges
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    for (const auto& cell1D : cell1Ds)
    {
        const unsigned int e = cell1D.second[0];

        mesh.Cell1DSetState(e, true);
        mesh.Cell1DSetMarker(e, cell1D.second[1]);
        mesh.Cell1DInsertExtremes(e,
                                  cell1D.first.first,
                                  cell1D.first.second);
    }

    /// <li> Set Faces
    for (const Cell2D& cell2D : cell2Ds)
    {
        const unsigned int f = cell2D.id;

        mesh.Cell2DSetState(f, true);
        mesh.Cell2DSetMarker(f, cell2D.marker);

        const unsigned int numCellVertices = cell2D.vertices.size();
        const unsigned int numCellEdges = cell2D.edges.size();

        Gedim::Output::Assert(numCellVertices == numCellEdges);

        mesh.Cell2DInitializeVertices(f, numCellVertices);
        mesh.Cell2DInitializeEdges(f, numCellEdges);

        for (unsigned int v = 0; v < numCellVertices; v++)
            mesh.Cell2DInsertVertex(f, v, cell2D.vertices[v]);

        for (unsigned int e = 0; e < numCellEdges; e++)
            mesh.Cell2DInsertEdge(f, e, cell2D.edges[e]);
    }

    /// <li> Set Cells
    unsigned int cells2DIndex = 0;
    for (unsigned int c = 0; c < numberOfCells3D; c++)
    {
        const unsigned int numVertices = cell3DVertices[c].size();
        const unsigned int numEdges = cellEdges[c].size();
        const unsigned int numFaces = cellFaces[c].size();

        if(numVertices == 0)
            continue;

        mesh.Cell3DSetState(cells2DIndex, true);
        mesh.Cell3DSetMarker(cells2DIndex, 0);

        mesh.Cell3DInitializeVertices(cells2DIndex, numVertices);
        mesh.Cell3DInitializeEdges(cells2DIndex, numEdges);
        mesh.Cell3DInitializeFaces(cells2DIndex, numFaces);

        for (unsigned int v = 0; v < numVertices; v++)
            mesh.Cell3DInsertVertex(cells2DIndex, v, cell3DVertices[c][v]);

        for (unsigned int e = 0; e < numEdges; e++)
            mesh.Cell3DInsertEdge(cells2DIndex, e, cellEdges[c][e]);

        for (unsigned int f = 0; f < numFaces; f++)
            mesh.Cell3DInsertFace(cells2DIndex, f, cellFaces[c][f]);

        cells2DIndex++;

    }
}
// ************************************************************************* //
void VoroInterface::GenerateRandomPoints(const Eigen::MatrixXd& polyhedronVertices,
                                         const unsigned int& numPoints,
                                         voro::container& con)
{
    const double x_min = polyhedronVertices.row(0).minCoeff(), x_max = polyhedronVertices.row(0).maxCoeff();
    const double y_min = polyhedronVertices.row(1).minCoeff(), y_max = polyhedronVertices.row(1).maxCoeff();
    const double z_min = polyhedronVertices.row(2).minCoeff(), z_max = polyhedronVertices.row(2).maxCoeff();

    srand(4);
    for(unsigned int i=0; i < numPoints; i++) {
        double x = x_min + rnd() * ( x_max - x_min );
        double y = y_min + rnd() * ( y_max - y_min );
        double z = z_min + rnd() * ( z_max - z_min );
        con.put(i, x, y, z);
    }
}
// ************************************************************************* //
void VoroInterface::GenerateCartesianPoints(const Eigen::MatrixXd& polyhedronVertices,
                                            const unsigned int& numPoints,
                                            voro::container& con)
{
    Eigen::Vector3d parallelepipedLengthTangent = polyhedronVertices.col(1) - polyhedronVertices.col(0);
    Eigen::Vector3d parallelepipedHeightTangent = polyhedronVertices.col(3) - polyhedronVertices.col(0);
    Eigen::Vector3d parallelepipedWidthTangent = polyhedronVertices.col(4) - polyhedronVertices.col(0);

    vector<double> lengthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1,
                                                                                             0.0, 1.0, 1);
    vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1,
                                                                                             0.0, 1.0, 1);
    vector<double> widthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1,
                                                                                            0.0, 1.0, 1);

    const unsigned int numLengthPoints = lengthMeshCurvilinearCoordinates.size();
    const unsigned int numHeightPoints = heightMeshCurvilinearCoordinates.size();
    const unsigned int numWidthPoints = widthMeshCurvilinearCoordinates.size();

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int w = 0; w < numWidthPoints; w++)
    {
        for (unsigned int h = 0; h < numHeightPoints; h++)
        {
            for (unsigned int l = 0; l < numLengthPoints; l++)
            {
                const Eigen::Vector3d coordinate = polyhedronVertices.col(0) +
                                                   lengthMeshCurvilinearCoordinates[l] * parallelepipedLengthTangent +
                                                   heightMeshCurvilinearCoordinates[h] * parallelepipedHeightTangent +
                                                   widthMeshCurvilinearCoordinates[w] * parallelepipedWidthTangent;
                con.put(cell0DIndex, coordinate(0), coordinate(1), coordinate(2));
                cell0DIndex++;
            }
        }
    }
}
#endif
// ************************************************************************* //
}

