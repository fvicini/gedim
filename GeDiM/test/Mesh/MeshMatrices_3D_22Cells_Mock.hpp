#ifndef __MeshMatrices_3D_22Cells_Mock_H
#define __MeshMatrices_3D_22Cells_Mock_H

#include "MeshMatrices.hpp"

using namespace std;

namespace GedimUnitTesting
{
/// \brief Mesh generated with Tetgen with parameter 0.5
class MeshMatrices_3D_22Cells_Mock final
{
  public:
    Gedim::MeshMatrices Mesh;

    MeshMatrices_3D_22Cells_Mock()
    {
        Mesh.Dimension = 3;
        Mesh.NumberCell0D = 20;
        Mesh.Cell0DCoordinates = {
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            1.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            5.0000000000000000e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
            5.0000000000000000e-01, 0.0000000000000000e+00, 5.0000000000000000e-01, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e-01, 0.0000000000000000e+00,
            5.0000000000000000e-01, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            5.0000000000000000e-01, 1.0000000000000000e+00, 5.0000000000000000e-01, 1.0000000000000000e+00,
            1.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e-01, 1.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 5.0000000000000000e-01, 0.0000000000000000e+00,
            0.0000000000000000e+00, 5.0000000000000000e-01, 1.0000000000000000e+00, 1.0000000000000000e+00,
            5.0000000000000000e-01, 1.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e-01};
        Mesh.Cell0DMarkers = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 17, 19, 18};
        Mesh.ActiveCell0D = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        Mesh.UpdatedCell0Ds = {};
        Mesh.NumberCell0DNeighbourCell1D = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Mesh.Cell0DNeighbourCell1Ds = {};
        Mesh.NumberCell0DNeighbourCell2D = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Mesh.Cell0DNeighbourCell2Ds = {};
        Mesh.NumberCell0DNeighbourCell3D.resize(Mesh.NumberCell0D + 1, 0);
        Mesh.Cell0DNeighbourCell3Ds = {};
        Mesh.Cell0DDoublePropertyIds = {};
        Mesh.Cell0DDoublePropertyIndices = {};
        Mesh.Cell0DDoublePropertySizes = {};
        Mesh.Cell0DDoublePropertyValues = {};
        Mesh.NumberCell1D = 59;
        Mesh.Cell1DVertices = {18, 10, 18, 14, 8,  14, 8,  10, 8,  18, 14, 10, 19, 12, 19, 8,  13, 8,  13, 12,
                               13, 19, 8,  12, 16, 10, 16, 8,  14, 16, 18, 9,  13, 9,  13, 18, 8,  9,  19, 9,
                               1,  8,  1,  9,  1,  19, 12, 14, 12, 15, 8,  15, 15, 14, 17, 11, 17, 0,  8,  0,
                               8,  11, 8,  17, 0,  11, 15, 11, 15, 17, 13, 14, 19, 5,  12, 5,  5,  13, 15, 16,
                               18, 6,  13, 6,  6,  14, 17, 12, 4,  12, 4,  15, 4,  17, 16, 11, 16, 7,  14, 7,
                               7,  15, 3,  10, 3,  11, 3,  16, 10, 11, 2,  9,  2,  10, 2,  18, 9,  10};
        Mesh.Cell1DMarkers = {26, 26, 0,  21, 0,  26, 25, 25, 0,  22, 24, 25, 26, 0,  26, 24, 24, 24, 21, 24,
                              9,  10, 18, 22, 22, 0,  22, 23, 17, 9,  21, 25, 12, 23, 23, 22, 18, 13, 14, 23,
                              19, 14, 15, 25, 13, 16, 17, 23, 20, 15, 16, 11, 12, 20, 21, 10, 11, 19, 21};
        Mesh.ActiveCell1D = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        Mesh.Cell1DOriginalCell1Ds.resize(Mesh.NumberCell1D, std::numeric_limits<unsigned int>::max());
        Mesh.UpdatedCell1Ds = {};
        Mesh.NumberCell1DNeighbourCell2D = {0,   4,   8,   14,  19,  23,  26,  29,  33,  38,  42,  46,  51,  55,  59,
                                            63,  67,  70,  74,  79,  82,  84,  86,  88,  91,  95,  100, 104, 107, 109,
                                            111, 116, 120, 122, 125, 129, 133, 135, 137, 139, 143, 145, 147, 149, 152,
                                            154, 156, 158, 162, 164, 166, 168, 170, 172, 174, 177, 179, 181, 183, 186};
        Mesh.Cell1DNeighbourCell2Ds = {
            0,  1,  57, 58, 1,  2,  31, 39, 2,  3,  10, 20, 23, 32, 0,  3,  9,  56, 61, 0,  2,  12, 13, 1,  3,  8,  4,
            5,  34, 5,  6,  16, 17, 6,  7,  13, 14, 32, 4,  7,  35, 38, 4,  6,  19, 33, 5,  7,  20, 22, 48, 8,  9,  53,
            54, 9,  10, 37, 46, 8,  10, 36, 50, 11, 12, 58, 59, 11, 14, 19, 11, 13, 31, 40, 12, 14, 16, 18, 61, 15, 16,
            19, 17, 18, 15, 18, 15, 17, 20, 21, 38, 21, 22, 43, 45, 22, 23, 29, 30, 37, 21, 23, 36, 51, 24, 25, 28, 25,
            26, 26, 27, 24, 27, 30, 46, 56, 24, 26, 29, 48, 25, 27, 28, 30, 47, 28, 29, 42, 43, 31, 32, 38, 41, 33, 34,
            34, 35, 33, 35, 36, 37, 47, 49, 39, 40, 40, 41, 39, 41, 43, 44, 48, 44, 45, 42, 45, 42, 44, 46, 47, 52, 53,
            49, 50, 50, 51, 49, 51, 54, 55, 52, 55, 52, 54, 53, 55, 56, 59, 60, 57, 60, 57, 59, 58, 60, 61};
        std::replace(Mesh.Cell1DNeighbourCell2Ds.begin(),
                     Mesh.Cell1DNeighbourCell2Ds.end(),
                     static_cast<unsigned int>(62),
                     std::numeric_limits<unsigned int>::max());
        Mesh.NumberCell1DNeighbourCell3D.resize(Mesh.NumberCell1D + 1, 0);
        Mesh.Cell1DNeighbourCell3Ds = {};
        Mesh.Cell1DDoublePropertyIds = {};
        Mesh.Cell1DDoublePropertyIndices = {};
        Mesh.Cell1DDoublePropertySizes = {};
        Mesh.Cell1DDoublePropertyValues = {};
        Mesh.NumberCell2D = 62;
        Mesh.NumberCell2DVertices = {0,   3,   6,   9,   12,  15,  18,  21,  24,  27,  30,  33,  36,  39,  42,  45,
                                     48,  51,  54,  57,  60,  63,  66,  69,  72,  75,  78,  81,  84,  87,  90,  93,
                                     96,  99,  102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141,
                                     144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186};
        Mesh.Cell2DVertices = {
            18, 10, 8,  18, 14, 10, 8,  14, 18, 8,  10, 14, 19, 12, 13, 19, 8,  12, 13, 8,  19, 13, 12, 8,  16, 10, 14,
            16, 8,  10, 14, 8,  16, 18, 9,  13, 18, 8,  9,  13, 8,  18, 13, 9,  8,  19, 9,  1,  19, 8,  9,  1,  8,  19,
            1,  9,  8,  19, 13, 9,  12, 14, 8,  12, 15, 14, 8,  15, 12, 8,  14, 15, 17, 11, 8,  17, 0,  11, 8,  0,  17,
            8,  11, 0,  17, 11, 15, 15, 8,  17, 15, 11, 8,  18, 13, 14, 8,  14, 13, 19, 5,  13, 12, 5,  19, 12, 13, 5,
            16, 14, 15, 15, 8,  16, 14, 13, 12, 18, 6,  14, 13, 6,  18, 13, 14, 6,  17, 15, 4,  17, 12, 15, 4,  12, 17,
            4,  15, 12, 16, 11, 8,  16, 15, 11, 17, 12, 8,  16, 7,  15, 14, 7,  16, 14, 15, 7,  16, 11, 3,  16, 10, 11,
            3,  10, 16, 3,  11, 10, 10, 11, 8,  18, 10, 2,  18, 9,  10, 2,  9,  18, 2,  10, 9,  9,  10, 8};
        Mesh.NumberCell2DEdges = {0,   3,   6,   9,   12,  15,  18,  21,  24,  27,  30,  33,  36,  39,  42,  45,
                                  48,  51,  54,  57,  60,  63,  66,  69,  72,  75,  78,  81,  84,  87,  90,  93,
                                  96,  99,  102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141,
                                  144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186};
        Mesh.Cell2DEdges = {
            0,  3,  4,  1,  5,  0,  2,  1,  4,  3,  5,  2,  6,  9,  10, 7,  11, 6,  8,  7,  10, 9,  11, 8,  12, 5,  14,
            13, 3,  12, 2,  13, 14, 15, 16, 17, 4,  18, 15, 8,  4,  17, 16, 18, 8,  19, 21, 22, 7,  18, 19, 20, 7,  22,
            21, 18, 20, 10, 16, 19, 23, 2,  11, 24, 26, 23, 25, 24, 11, 2,  26, 25, 27, 30, 31, 28, 32, 27, 29, 28, 31,
            30, 32, 29, 27, 33, 34, 25, 31, 34, 33, 30, 25, 17, 35, 1,  2,  35, 8,  36, 38, 10, 37, 36, 6,  9,  38, 37,
            14, 26, 39, 25, 13, 39, 35, 9,  23, 40, 42, 1,  41, 40, 17, 35, 42, 41, 34, 45, 46, 43, 24, 34, 44, 43, 46,
            45, 24, 44, 47, 30, 13, 39, 33, 47, 43, 11, 31, 48, 50, 39, 49, 48, 14, 26, 50, 49, 47, 52, 53, 12, 54, 47,
            51, 12, 53, 52, 54, 51, 54, 30, 3,  0,  56, 57, 15, 58, 0,  55, 15, 57, 56, 58, 55, 58, 3,  18};
        Mesh.NumberCell2DSubdivision = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Mesh.Cell2DSubdivision = {};
        Mesh.Cell2DMarkers = {0,  26, 0,  0,  0,  25, 0,  0,  26, 0,  0,  24, 0,  0,  0,  24, 0, 25, 21, 24, 0,
                              22, 0,  0,  0,  23, 25, 21, 23, 0,  0,  0,  0,  24, 25, 22, 0,  0, 22, 26, 24, 22,
                              23, 0,  25, 22, 0,  23, 25, 23, 26, 22, 23, 0,  26, 21, 21, 26, 0, 24, 21, 21};
        Mesh.ActiveCell2D = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        Mesh.Cell2DOriginalCell2Ds.resize(Mesh.NumberCell2D, std::numeric_limits<unsigned int>::max());
        Mesh.UpdatedCell2Ds = {};
        Mesh.Cell2DDoublePropertyIds = {};
        Mesh.Cell2DDoublePropertyIndices = {};
        Mesh.Cell2DDoublePropertySizes = {};
        Mesh.Cell2DDoublePropertyValues = {};
        Mesh.NumberCell2DNeighbourCell3D = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Mesh.Cell2DNeighbourCell3Ds = {};
        Mesh.NumberCell3D = 22;
        Mesh.NumberCell3DVertices = {0,  4,  8,  12, 16, 20, 24, 28, 32, 36, 40, 44,
                                     48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88};
        Mesh.Cell3DVertices = {14, 8,  10, 18, 8,  13, 12, 19, 8,  14, 10, 16, 8,  13, 9,  18, 8,  1,  9,  19, 13, 8,
                               9,  19, 15, 8,  14, 12, 0,  8,  11, 17, 8,  15, 11, 17, 13, 8,  14, 18, 5,  12, 13, 19,
                               8,  15, 14, 16, 13, 8,  12, 14, 6,  13, 14, 18, 12, 4,  15, 17, 15, 8,  11, 16, 15, 8,
                               12, 17, 7,  14, 15, 16, 10, 3,  11, 16, 8,  10, 11, 16, 9,  2,  10, 18, 8,  9,  10, 18};
        Mesh.NumberCell3DEdges = {0,  6,  12, 18, 24, 30,  36,  42,  48,  54,  60, 66,
                                  72, 78, 84, 90, 96, 102, 108, 114, 120, 126, 132};
        Mesh.Cell3DEdges = {3,  5,  2,  0,  4,  1,  9,  11, 8,  6,  10, 7,  3,  5,  2,  12, 14, 13, 16, 18, 8,  15,
                            17, 4,  21, 18, 20, 19, 22, 7,  16, 18, 8,  7,  19, 10, 2,  26, 25, 23, 11, 24, 30, 32,
                            29, 27, 31, 28, 33, 30, 25, 27, 34, 31, 2,  35, 8,  1,  4,  17, 9,  38, 37, 6,  10, 36,
                            2,  26, 25, 14, 39, 13, 9,  11, 8,  23, 2,  35, 35, 42, 41, 17, 1,  40, 45, 24, 44, 34,
                            46, 43, 33, 30, 25, 47, 13, 39, 25, 24, 11, 43, 31, 34, 26, 50, 49, 14, 39, 48, 52, 54,
                            51, 47, 53, 12, 54, 30, 3,  12, 47, 13, 56, 58, 55, 0,  57, 15, 58, 3,  18, 15, 0,  4};
        Mesh.NumberCell3DFaces = {0,  4,  8,  12, 16, 20, 24, 28, 32, 36, 40, 44,
                                  48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88};
        Mesh.Cell3DFaces = {3,  0,  1,  2,  7,  4,  5,  6,  3,  8,  9,  10, 14, 11, 12, 13, 18, 15, 16, 17, 14, 16,
                            19, 6,  23, 20, 21, 22, 27, 24, 25, 26, 30, 28, 24, 29, 32, 2,  31, 13, 35, 4,  33, 34,
                            23, 36, 10, 37, 7,  20, 38, 32, 41, 31, 39, 40, 45, 42, 43, 44, 30, 46, 47, 37, 22, 48,
                            43, 29, 51, 36, 49, 50, 55, 52, 53, 54, 56, 53, 46, 9,  60, 57, 58, 59, 61, 58, 0,  12};
        Mesh.Cell3DMarkers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Mesh.ActiveCell3D = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        Mesh.Cell3DOriginalCell3Ds.resize(Mesh.NumberCell3D, std::numeric_limits<unsigned int>::max());
        Mesh.UpdatedCell3Ds = {};
        Mesh.Cell3DDoublePropertyIds = {};
        Mesh.Cell3DDoublePropertyIndices = {};
        Mesh.Cell3DDoublePropertySizes = {};
        Mesh.Cell3DDoublePropertyValues = {};
    }
};
} // namespace GedimUnitTesting

#endif // __MeshMatrices_3D_22Cells_Mock_H
