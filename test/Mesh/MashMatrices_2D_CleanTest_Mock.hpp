#ifndef __MashMatrices_2D_CleanTest_Mock_H
#define __MashMatrices_2D_CleanTest_Mock_H

#include "MeshMatrices.hpp"

using namespace std;

namespace GedimUnitTesting
{
  class MashMatrices_2D_CleanTest_Mock final
  {
    public:
      Gedim::MeshMatrices Mesh;

      MashMatrices_2D_CleanTest_Mock() {
        Mesh.Dimension = 2;
        Mesh.NumberCell0D = 53;
        Mesh.Cell0DCoordinates = {0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,3.5355338217215095e+00,1.0130785099704553e-15,0.0000000000000000e+00,3.5355337874108588e+00,2.3570225890964744e+00,0.0000000000000000e+00,7.5804043805915988e-09,2.3570226012651645e+00,0.0000000000000000e+00,1.8543476199006648e-09,5.7658391606561521e-01,0.0000000000000000e+00,2.0800404195068691e+00,5.5511151231257827e-16,0.0000000000000000e+00,2.1151591227841232e-01,5.1795204001768336e-01,0.0000000000000000e+00,3.3521159490867891e-01,4.8366379230571849e-01,0.0000000000000000e+00,1.4099316784565923e+00,1.8575308384662059e-01,0.0000000000000000e+00,1.8871664326627764e+00,5.3464364296950673e-02,0.0000000000000000e+00,3.5355338032537196e+00,1.2686730725730853e+00,0.0000000000000000e+00,2.0943046956479745e+00,2.0617074849731490e-01,0.0000000000000000e+00,2.1376853971985486e+00,2.3815185127568803e-01,0.0000000000000000e+00,2.1387026992538010e+00,2.3890182627431167e-01,0.0000000000000000e+00,2.2345242120255202e+00,3.0954332211916136e-01,0.0000000000000000e+00,2.3413215029355299e+00,3.8827637656418884e-01,0.0000000000000000e+00,2.4600756404097117e+00,4.7582425246682269e-01,0.0000000000000000e+00,2.5750410530181331e+00,5.6057900619761736e-01,0.0000000000000000e+00,2.6508407685105349e+00,6.1646004096852358e-01,0.0000000000000000e+00,2.7762345047921038e+00,7.0890275914609546e-01,0.0000000000000000e+00,2.8432311695860704e+00,7.5829401253756357e-01,0.0000000000000000e+00,3.4517584544032860e+00,1.2069122446907541e+00,0.0000000000000000e+00,3.5355337903192301e+00,2.1572275557136447e+00,0.0000000000000000e+00,3.1038529036714140e+00,1.9314860693452975e+00,0.0000000000000000e+00,2.6315482717016723e+00,1.6845009631307413e+00,0.0000000000000000e+00,2.3041916345247953e+00,1.5133143849555970e+00,0.0000000000000000e+00,2.2272498520444808e+00,1.4730787551694060e+00,0.0000000000000000e+00,2.1915091550528762e+00,1.4543886572005613e+00,0.0000000000000000e+00,1.9543090996102959e+00,1.3303482033874010e+00,0.0000000000000000e+00,1.8598288684546722e+00,1.2809410869891811e+00,0.0000000000000000e+00,1.8206193946448963e+00,1.2604370406639593e+00,0.0000000000000000e+00,1.7776480583129837e+00,1.2379657809097275e+00,0.0000000000000000e+00,1.7685144508252053e+00,1.2331894885280863e+00,0.0000000000000000e+00,1.7348480802286366e+00,1.2155841305561050e+00,0.0000000000000000e+00,1.6609879375780991e+00,1.1769600012489549e+00,0.0000000000000000e+00,5.5847096325525802e-01,6.0041415864133185e-01,0.0000000000000000e+00,3.5327816586370329e+00,9.9920072216264089e-16,0.0000000000000000e+00,3.1059988779366656e+00,2.9831311079983275e-01,0.0000000000000000e+00,2.5511754023576092e+00,6.8612425133636679e-01,0.0000000000000000e+00,2.5364177919433022e+00,6.9643954255650309e-01,0.0000000000000000e+00,2.4068652671461686e+00,7.8699431287967236e-01,0.0000000000000000e+00,2.2720510543939296e+00,8.8122690380134316e-01,0.0000000000000000e+00,2.1885577753059100e+00,9.3958713045986020e-01,0.0000000000000000e+00,2.0796979744760629e+00,1.0156780771782592e+00,0.0000000000000000e+00,1.9187554741459347e+00,1.1281738481055568e+00,0.0000000000000000e+00,1.8043588218012616e+00,1.2081349495586757e+00,0.0000000000000000e+00,1.6069835660445353e-01,2.3570226007120687e+00,0.0000000000000000e+00,1.7109228464639523e+00,1.2734449332879820e+00,0.0000000000000000e+00,1.6334541357544190e+00,1.3275941003267477e+00,0.0000000000000000e+00,1.5860898971780075e+00,1.3607008093818973e+00,0.0000000000000000e+00,1.5144373172539451e+00,1.4107846107454471e+00,0.0000000000000000e+00,1.1734756542422693e+00,1.6491103796285811e+00,0.0000000000000000e+00,9.9770224279681452e-01,1.7719726775829261e+00,0.0000000000000000e+00};
        Mesh.Cell0DMarkers = {2,2,2,2,2,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0};
        Mesh.ActiveCell0D = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        Mesh.UpdatedCell0Ds = {};
        Mesh.NumberCell0DNeighbourCell1D = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        Mesh.Cell0DNeighbourCell1Ds = {};
        Mesh.NumberCell0DNeighbourCell2D = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        Mesh.Cell0DNeighbourCell2Ds = {};
        Mesh.Cell0DDoublePropertyIds = {};
        Mesh.Cell0DDoublePropertyIndices = {};
        Mesh.Cell0DDoublePropertySizes = {};
        Mesh.Cell0DDoublePropertyValues = {};
        Mesh.NumberCell1D = 59;
        Mesh.Cell1DVertices = {3,4,4,0,0,5,4,6,6,7,7,8,8,9,9,5,1,10,9,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,10,10,22,22,2,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,7,5,36,36,1,36,37,37,18,18,38,38,39,39,40,40,41,41,42,42,43,43,44,44,45,45,32,46,3,2,46,32,47,47,48,48,49,49,50,50,51,51,52,52,46};
        Mesh.Cell1DAdjacency.resize(53, 53);
        Mesh.Cell1DAdjacency.reserve(59);
        Mesh.Cell1DAdjacency.insert(4, 0) = 2;
        Mesh.Cell1DAdjacency.insert(36, 1) = 39;
        Mesh.Cell1DAdjacency.insert(22, 2) = 23;
        Mesh.Cell1DAdjacency.insert(46, 3) = 51;
        Mesh.Cell1DAdjacency.insert(3, 4) = 1;
        Mesh.Cell1DAdjacency.insert(0, 5) = 3;
        Mesh.Cell1DAdjacency.insert(9, 5) = 8;
        Mesh.Cell1DAdjacency.insert(4, 6) = 4;
        Mesh.Cell1DAdjacency.insert(6, 7) = 5;
        Mesh.Cell1DAdjacency.insert(35, 7) = 37;
        Mesh.Cell1DAdjacency.insert(7, 8) = 6;
        Mesh.Cell1DAdjacency.insert(8, 9) = 7;
        Mesh.Cell1DAdjacency.insert(1, 10) = 9;
        Mesh.Cell1DAdjacency.insert(21, 10) = 21;
        Mesh.Cell1DAdjacency.insert(9, 11) = 10;
        Mesh.Cell1DAdjacency.insert(11, 12) = 11;
        Mesh.Cell1DAdjacency.insert(12, 13) = 12;
        Mesh.Cell1DAdjacency.insert(13, 14) = 13;
        Mesh.Cell1DAdjacency.insert(14, 15) = 14;
        Mesh.Cell1DAdjacency.insert(15, 16) = 15;
        Mesh.Cell1DAdjacency.insert(16, 17) = 16;
        Mesh.Cell1DAdjacency.insert(17, 18) = 17;
        Mesh.Cell1DAdjacency.insert(37, 18) = 41;
        Mesh.Cell1DAdjacency.insert(18, 19) = 18;
        Mesh.Cell1DAdjacency.insert(19, 20) = 19;
        Mesh.Cell1DAdjacency.insert(20, 21) = 20;
        Mesh.Cell1DAdjacency.insert(10, 22) = 22;
        Mesh.Cell1DAdjacency.insert(22, 23) = 24;
        Mesh.Cell1DAdjacency.insert(23, 24) = 25;
        Mesh.Cell1DAdjacency.insert(24, 25) = 26;
        Mesh.Cell1DAdjacency.insert(25, 26) = 27;
        Mesh.Cell1DAdjacency.insert(26, 27) = 28;
        Mesh.Cell1DAdjacency.insert(27, 28) = 29;
        Mesh.Cell1DAdjacency.insert(28, 29) = 30;
        Mesh.Cell1DAdjacency.insert(29, 30) = 31;
        Mesh.Cell1DAdjacency.insert(30, 31) = 32;
        Mesh.Cell1DAdjacency.insert(31, 32) = 33;
        Mesh.Cell1DAdjacency.insert(45, 32) = 50;
        Mesh.Cell1DAdjacency.insert(32, 33) = 34;
        Mesh.Cell1DAdjacency.insert(33, 34) = 35;
        Mesh.Cell1DAdjacency.insert(34, 35) = 36;
        Mesh.Cell1DAdjacency.insert(5, 36) = 38;
        Mesh.Cell1DAdjacency.insert(36, 37) = 40;
        Mesh.Cell1DAdjacency.insert(18, 38) = 42;
        Mesh.Cell1DAdjacency.insert(38, 39) = 43;
        Mesh.Cell1DAdjacency.insert(39, 40) = 44;
        Mesh.Cell1DAdjacency.insert(40, 41) = 45;
        Mesh.Cell1DAdjacency.insert(41, 42) = 46;
        Mesh.Cell1DAdjacency.insert(42, 43) = 47;
        Mesh.Cell1DAdjacency.insert(43, 44) = 48;
        Mesh.Cell1DAdjacency.insert(44, 45) = 49;
        Mesh.Cell1DAdjacency.insert(2, 46) = 52;
        Mesh.Cell1DAdjacency.insert(52, 46) = 59;
        Mesh.Cell1DAdjacency.insert(32, 47) = 53;
        Mesh.Cell1DAdjacency.insert(47, 48) = 54;
        Mesh.Cell1DAdjacency.insert(48, 49) = 55;
        Mesh.Cell1DAdjacency.insert(49, 50) = 56;
        Mesh.Cell1DAdjacency.insert(50, 51) = 57;
        Mesh.Cell1DAdjacency.insert(51, 52) = 58;
        Mesh.Cell1DAdjacency.makeCompressed();
        Mesh.Cell1DMarkers = {2,2,2,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0};
        Mesh.ActiveCell1D = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        Mesh.UpdatedCell1Ds = {};
        Mesh.NumberCell1DNeighbourCell2D = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118};
        Mesh.Cell1DNeighbourCell2Ds = {7,6,7,0,7,0,0,6,0,6,0,3,0,3,0,2,7,1,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,1,4,1,4,1,4,1,4,7,4,7,5,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,6,3,6,3,6,3,6,3,7,2,7,1,1,2,1,2,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,7,6,7,5,5,6,5,6,5,6,5,6,5,6,5,6,5,6};
        Mesh.Cell1DDoublePropertyIds = {};
        Mesh.Cell1DDoublePropertyIndices = {};
        Mesh.Cell1DDoublePropertySizes = {};
        Mesh.Cell1DDoublePropertyValues = {};
        Mesh.NumberCell2D = 7;
        Mesh.NumberCell2DVertices = {0,7,15,27,50,74,93,108};
        Mesh.Cell2DVertices = {0,5,9,8,7,6,4,1,10,21,20,19,18,37,36,17,16,15,14,13,12,11,9,5,36,37,18,8,9,11,12,13,14,15,16,17,18,38,39,40,41,42,43,44,45,32,33,34,35,7,19,20,21,10,22,23,24,25,26,27,28,29,30,31,32,45,44,43,42,41,40,39,38,18,2,46,52,51,50,49,48,47,32,31,30,29,28,27,26,25,24,23,22,3,4,6,7,35,34,33,32,47,48,49,50,51,52,46};
        Mesh.NumberCell2DEdges = {0,7,15,27,50,74,93,108};
        Mesh.Cell2DEdges = {2,7,6,5,4,3,1,8,20,19,18,17,40,39,38,15,14,13,12,11,10,9,7,37,39,40,16,6,9,10,11,12,13,14,15,16,41,42,43,44,45,46,47,48,49,33,34,35,36,5,18,19,20,21,23,24,25,26,27,28,29,30,31,32,49,48,47,46,45,44,43,42,41,17,51,58,57,56,55,54,53,52,32,31,30,29,28,27,26,25,24,23,22,0,3,4,36,35,34,33,52,53,54,55,56,57,58,50};
        Mesh.NumberCell2DSubdivision = {0,0,0,0,0,0,0,0};
        Mesh.Cell2DSubdivision = {};
        Mesh.Cell2DMarkers = {0,0,0,0,0,0,0};
        Mesh.ActiveCell2D = {1,1,1,1,1,1,1};
        Mesh.UpdatedCell2Ds = {};
        Mesh.Cell2DDoublePropertyIds = {};
        Mesh.Cell2DDoublePropertyIndices = {};
        Mesh.Cell2DDoublePropertySizes = {};
        Mesh.Cell2DDoublePropertyValues = {};
        Mesh.NumberCell3D = 0;
        Mesh.NumberCell3DVertices = {};
        Mesh.Cell3DVertices = {};
        Mesh.NumberCell3DEdges = {};
        Mesh.Cell3DEdges = {};
        Mesh.NumberCell3DFaces = {};
        Mesh.Cell3DFaces = {};
        Mesh.Cell3DMarkers = {};
        Mesh.ActiveCell3D = {};
        Mesh.UpdatedCell3Ds = {};
        Mesh.Cell3DDoublePropertyIds = {};
        Mesh.Cell3DDoublePropertyIndices = {};
        Mesh.Cell3DDoublePropertySizes = {};
        Mesh.Cell3DDoublePropertyValues = {};
      }
  };
}

#endif // __MashMatrices_2D_CleanTest_Mock_H
