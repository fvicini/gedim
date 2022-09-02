#ifndef __Quadrature_Gauss1D_H
#define __Quadrature_Gauss1D_H

#include "Eigen/Eigen"

namespace Gedim
{
  /// Gauss quadrature rule for segments
  class Quadrature_Gauss1D
  {
    public:
      static void FillPointsAndWeights(const unsigned int& order,
                                       Eigen::MatrixXd& points,
                                       Eigen::VectorXd& weights)
      {
        switch(order)
        {
          case 0:
          case 1:
          {
            points.setZero(3,1);
            points(0,0) = 5.0000000000000000e-01;
            weights.resize(1);
            weights[0] = 1.0000000000000000e+00;
          }
            break;
          case 2:
          case 3:
          {
            points.setZero(3,2);
            points(0,0) = 2.1132486540518713e-01;
            points(0,1) = 7.8867513459481287e-01;
            weights.resize(2);
            weights[0] = 5.0000000000000000e-01;
            weights[1] = 5.0000000000000000e-01;
          }
            break;
          case 4:
          case 5:
          {
            points.setZero(3,3);
            points(0,0) = 1.1270166537925830e-01;
            points(0,1) = 5.0000000000000000e-01;
            points(0,2) = 8.8729833462074170e-01;
            weights.resize(3);
            weights[0] = 2.7777777777777779e-01;
            weights[1] = 4.4444444444444442e-01;
            weights[2] = 2.7777777777777779e-01;
          }
            break;
          case 6:
          case 7:
          {
            points.setZero(3,4);
            points(0,0) = 6.9431844202973714e-02;
            points(0,1) = 3.3000947820757187e-01;
            points(0,2) = 6.6999052179242813e-01;
            points(0,3) = 9.3056815579702623e-01;
            weights.resize(4);
            weights[0] = 1.7392742256872692e-01;
            weights[1] = 3.2607257743127305e-01;
            weights[2] = 3.2607257743127305e-01;
            weights[3] = 1.7392742256872692e-01;
          }
            break;
          case 8:
          case 9:
          {
            points.setZero(3,5);
            points(0,0) = 4.6910077030668018e-02;
            points(0,1) = 2.3076534494715845e-01;
            points(0,2) = 5.0000000000000000e-01;
            points(0,3) = 7.6923465505284150e-01;
            points(0,4) = 9.5308992296933193e-01;
            weights.resize(5);
            weights[0] = 1.1846344252809454e-01;
            weights[1] = 2.3931433524968324e-01;
            weights[2] = 2.8444444444444444e-01;
            weights[3] = 2.3931433524968324e-01;
            weights[4] = 1.1846344252809454e-01;
          }
            break;
          case 10:
          case 11:
          {
            points.setZero(3,6);
            points(0,0) = 3.3765242898423975e-02;
            points(0,1) = 1.6939530676686776e-01;
            points(0,2) = 3.8069040695840156e-01;
            points(0,3) = 6.1930959304159849e-01;
            points(0,4) = 8.3060469323313224e-01;
            points(0,5) = 9.6623475710157603e-01;
            weights.resize(6);
            weights[0] = 8.5662246189585178e-02;
            weights[1] = 1.8038078652406930e-01;
            weights[2] = 2.3395696728634552e-01;
            weights[3] = 2.3395696728634552e-01;
            weights[4] = 1.8038078652406930e-01;
            weights[5] = 8.5662246189585178e-02;
          }
            break;
          case 12:
          case 13:
          {
            points.setZero(3,7);
            points(0,0) = 2.5446043828620757e-02;
            points(0,1) = 1.2923440720030277e-01;
            points(0,2) = 2.9707742431130141e-01;
            points(0,3) = 5.0000000000000000e-01;
            points(0,4) = 7.0292257568869854e-01;
            points(0,5) = 8.7076559279969723e-01;
            points(0,6) = 9.7455395617137919e-01;
            weights.resize(7);
            weights[0] = 6.4742483084434851e-02;
            weights[1] = 1.3985269574463832e-01;
            weights[2] = 1.9091502525255946e-01;
            weights[3] = 2.0897959183673470e-01;
            weights[4] = 1.9091502525255946e-01;
            weights[5] = 1.3985269574463832e-01;
            weights[6] = 6.4742483084434851e-02;
          }
            break;
          case 14:
          case 15:
          {
            points.setZero(3,8);
            points(0,0) = 1.9855071751231912e-02;
            points(0,1) = 1.0166676129318664e-01;
            points(0,2) = 2.3723379504183550e-01;
            points(0,3) = 4.0828267875217511e-01;
            points(0,4) = 5.9171732124782495e-01;
            points(0,5) = 7.6276620495816450e-01;
            points(0,6) = 8.9833323870681336e-01;
            points(0,7) = 9.8014492824876809e-01;
            weights.resize(8);
            weights[0] = 5.0614268145188129e-02;
            weights[1] = 1.1119051722668724e-01;
            weights[2] = 1.5685332293894363e-01;
            weights[3] = 1.8134189168918100e-01;
            weights[4] = 1.8134189168918100e-01;
            weights[5] = 1.5685332293894363e-01;
            weights[6] = 1.1119051722668724e-01;
            weights[7] = 5.0614268145188129e-02;
          }
            break;
          case 16:
          case 17:
          {
            points.setZero(3,9);
            points(0,0) = 1.5919880246186957e-02;
            points(0,1) = 8.1984446336682115e-02;
            points(0,2) = 1.9331428364970482e-01;
            points(0,3) = 3.3787328829809554e-01;
            points(0,4) = 5.0000000000000000e-01;
            points(0,5) = 6.6212671170190451e-01;
            points(0,6) = 8.0668571635029518e-01;
            points(0,7) = 9.1801555366331788e-01;
            points(0,8) = 9.8408011975381304e-01;
            weights.resize(9);
            weights[0] = 4.0637194180787206e-02;
            weights[1] = 9.0324080347428698e-02;
            weights[2] = 1.3030534820146772e-01;
            weights[3] = 1.5617353852000143e-01;
            weights[4] = 1.6511967750062989e-01;
            weights[5] = 1.5617353852000143e-01;
            weights[6] = 1.3030534820146772e-01;
            weights[7] = 9.0324080347428698e-02;
            weights[8] = 4.0637194180787206e-02;
          }
            break;
          case 18:
          case 19:
          {
            points.setZero(3,10);
            points(0,0) = 1.3046735741414128e-02;
            points(0,1) = 6.7468316655507732e-02;
            points(0,2) = 1.6029521585048778e-01;
            points(0,3) = 2.8330230293537639e-01;
            points(0,4) = 4.2556283050918442e-01;
            points(0,5) = 5.7443716949081558e-01;
            points(0,6) = 7.1669769706462361e-01;
            points(0,7) = 8.3970478414951222e-01;
            points(0,8) = 9.3253168334449232e-01;
            points(0,9) = 9.8695326425858587e-01;
            weights.resize(10);
            weights[0] = 3.3335672154344069e-02;
            weights[1] = 7.4725674575290293e-02;
            weights[2] = 1.0954318125799102e-01;
            weights[3] = 1.3463335965499817e-01;
            weights[4] = 1.4776211235737644e-01;
            weights[5] = 1.4776211235737644e-01;
            weights[6] = 1.3463335965499817e-01;
            weights[7] = 1.0954318125799102e-01;
            weights[8] = 7.4725674575290293e-02;
            weights[9] = 3.3335672154344069e-02;
          }
            break;
          case 20:
          case 21:
          {
            points.setZero(3,11);
            points(0,0) = 1.0885670926971514e-02;
            points(0,1) = 5.6468700115952342e-02;
            points(0,2) = 1.3492399721297532e-01;
            points(0,3) = 2.4045193539659410e-01;
            points(0,4) = 3.6522842202382755e-01;
            points(0,5) = 5.0000000000000000e-01;
            points(0,6) = 6.3477157797617245e-01;
            points(0,7) = 7.5954806460340585e-01;
            points(0,8) = 8.6507600278702468e-01;
            points(0,9) = 9.4353129988404771e-01;
            points(0,10) = 9.8911432907302843e-01;
            weights.resize(11);
            weights[0] = 2.7834283558086832e-02;
            weights[1] = 6.2790184732452306e-02;
            weights[2] = 9.3145105463867131e-02;
            weights[3] = 1.1659688229599524e-01;
            weights[4] = 1.3140227225512333e-01;
            weights[5] = 1.3646254338895031e-01;
            weights[6] = 1.3140227225512333e-01;
            weights[7] = 1.1659688229599524e-01;
            weights[8] = 9.3145105463867131e-02;
            weights[9] = 6.2790184732452306e-02;
            weights[10] = 2.7834283558086832e-02;
          }
            break;
          case 22:
          case 23:
          {
            points.setZero(3,12);
            points(0,0) = 9.2196828766403782e-03;
            points(0,1) = 4.7941371814762546e-02;
            points(0,2) = 1.1504866290284765e-01;
            points(0,3) = 2.0634102285669126e-01;
            points(0,4) = 3.1608425050090994e-01;
            points(0,5) = 4.3738329574426554e-01;
            points(0,6) = 5.6261670425573440e-01;
            points(0,7) = 6.8391574949909006e-01;
            points(0,8) = 7.9365897714330869e-01;
            points(0,9) = 8.8495133709715235e-01;
            points(0,10) = 9.5205862818523745e-01;
            points(0,11) = 9.9078031712335957e-01;
            weights.resize(12);
            weights[0] = 2.3587668193255914e-02;
            weights[1] = 5.3469662997659213e-02;
            weights[2] = 8.0039164271673111e-02;
            weights[3] = 1.0158371336153296e-01;
            weights[4] = 1.1674626826917740e-01;
            weights[5] = 1.2457352290670139e-01;
            weights[6] = 1.2457352290670139e-01;
            weights[7] = 1.1674626826917740e-01;
            weights[8] = 1.0158371336153296e-01;
            weights[9] = 8.0039164271673111e-02;
            weights[10] = 5.3469662997659213e-02;
            weights[11] = 2.3587668193255914e-02;
          }
            break;
          case 24:
          case 25:
          {
            points.setZero(3,13);
            points(0,0) = 7.9084726407059325e-03;
            points(0,1) = 4.1200800388511039e-02;
            points(0,2) = 9.9210954633345061e-02;
            points(0,3) = 1.7882533027982989e-01;
            points(0,4) = 2.7575362448177654e-01;
            points(0,5) = 3.8477084202243261e-01;
            points(0,6) = 5.0000000000000000e-01;
            points(0,7) = 6.1522915797756739e-01;
            points(0,8) = 7.2424637551822346e-01;
            points(0,9) = 8.2117466972017006e-01;
            points(0,10) = 9.0078904536665494e-01;
            points(0,11) = 9.5879919961148896e-01;
            points(0,12) = 9.9209152735929407e-01;
            weights.resize(13);
            weights[0] = 2.0242002382657939e-02;
            weights[1] = 4.6060749918864226e-02;
            weights[2] = 6.9436755109893625e-02;
            weights[3] = 8.9072990380972869e-02;
            weights[4] = 1.0390802376844425e-01;
            weights[5] = 1.1314159013144862e-01;
            weights[6] = 1.1627577661543695e-01;
            weights[7] = 1.1314159013144862e-01;
            weights[8] = 1.0390802376844425e-01;
            weights[9] = 8.9072990380972869e-02;
            weights[10] = 6.9436755109893625e-02;
            weights[11] = 4.6060749918864226e-02;
            weights[12] = 2.0242002382657939e-02;
          }
            break;
          case 26:
          case 27:
          {
            points.setZero(3,14);
            points(0,0) = 6.8580956515938429e-03;
            points(0,1) = 3.5782558168213241e-02;
            points(0,2) = 8.6399342465117490e-02;
            points(0,3) = 1.5635354759415726e-01;
            points(0,4) = 2.4237568182092295e-01;
            points(0,5) = 3.4044381553605513e-01;
            points(0,6) = 4.4597252564632817e-01;
            points(0,7) = 5.5402747435367183e-01;
            points(0,8) = 6.5955618446394482e-01;
            points(0,9) = 7.5762431817907705e-01;
            points(0,10) = 8.4364645240584268e-01;
            points(0,11) = 9.1360065753488251e-01;
            points(0,12) = 9.6421744183178681e-01;
            points(0,13) = 9.9314190434840621e-01;
            weights.resize(14);
            weights[0] = 1.7559730165875930e-02;
            weights[1] = 4.0079043579880104e-02;
            weights[2] = 6.0759285343951593e-02;
            weights[3] = 7.8601583579096773e-02;
            weights[4] = 9.2769198738968911e-02;
            weights[5] = 1.0259923186064780e-01;
            weights[6] = 1.0763192673157890e-01;
            weights[7] = 1.0763192673157890e-01;
            weights[8] = 1.0259923186064780e-01;
            weights[9] = 9.2769198738968911e-02;
            weights[10] = 7.8601583579096773e-02;
            weights[11] = 6.0759285343951593e-02;
            weights[12] = 4.0079043579880104e-02;
            weights[13] = 1.7559730165875930e-02;
          }
            break;
          case 28:
          case 29:
          {
            points.setZero(3,15);
            points(0,0) = 6.0037409897573113e-03;
            points(0,1) = 3.1363303799647024e-02;
            points(0,2) = 7.5896708294786397e-02;
            points(0,3) = 1.3779113431991497e-01;
            points(0,4) = 2.1451391369573058e-01;
            points(0,5) = 3.0292432646121831e-01;
            points(0,6) = 3.9940295300128276e-01;
            points(0,7) = 5.0000000000000000e-01;
            points(0,8) = 6.0059704699871730e-01;
            points(0,9) = 6.9707567353878175e-01;
            points(0,10) = 7.8548608630426942e-01;
            points(0,11) = 8.6220886568008503e-01;
            points(0,12) = 9.2410329170521366e-01;
            points(0,13) = 9.6863669620035298e-01;
            points(0,14) = 9.9399625901024269e-01;
            weights.resize(15);
            weights[0] = 1.5376620998058635e-02;
            weights[1] = 3.5183023744054062e-02;
            weights[2] = 5.3579610233585970e-02;
            weights[3] = 6.9785338963077162e-02;
            weights[4] = 8.3134602908496960e-02;
            weights[5] = 9.3080500007781106e-02;
            weights[6] = 9.9215742663555789e-02;
            weights[7] = 1.0128912096278064e-01;
            weights[8] = 9.9215742663555789e-02;
            weights[9] = 9.3080500007781106e-02;
            weights[10] = 8.3134602908496960e-02;
            weights[11] = 6.9785338963077162e-02;
            weights[12] = 5.3579610233585970e-02;
            weights[13] = 3.5183023744054062e-02;
            weights[14] = 1.5376620998058635e-02;
          }
            break;
          case 30:
          case 31:
          {
            points.setZero(3,16);
            points(0,0) = 5.2995325041750307e-03;
            points(0,1) = 2.7712488463383700e-02;
            points(0,2) = 6.7184398806084122e-02;
            points(0,3) = 1.2229779582249850e-01;
            points(0,4) = 1.9106187779867811e-01;
            points(0,5) = 2.7099161117138632e-01;
            points(0,6) = 3.5919822461037054e-01;
            points(0,7) = 4.5249374508118129e-01;
            points(0,8) = 5.4750625491881877e-01;
            points(0,9) = 6.4080177538962946e-01;
            points(0,10) = 7.2900838882861363e-01;
            points(0,11) = 8.0893812220132189e-01;
            points(0,12) = 8.7770220417750155e-01;
            points(0,13) = 9.3281560119391593e-01;
            points(0,14) = 9.7228751153661630e-01;
            points(0,15) = 9.9470046749582497e-01;
            weights.resize(16);
            weights[0] = 1.3576229705877048e-02;
            weights[1] = 3.1126761969323947e-02;
            weights[2] = 4.7579255841246393e-02;
            weights[3] = 6.2314485627766938e-02;
            weights[4] = 7.4797994408288368e-02;
            weights[5] = 8.4578259697501268e-02;
            weights[6] = 9.1301707522461792e-02;
            weights[7] = 9.4725305227534251e-02;
            weights[8] = 9.4725305227534251e-02;
            weights[9] = 9.1301707522461792e-02;
            weights[10] = 8.4578259697501268e-02;
            weights[11] = 7.4797994408288368e-02;
            weights[12] = 6.2314485627766938e-02;
            weights[13] = 4.7579255841246393e-02;
            weights[14] = 3.1126761969323947e-02;
            weights[15] = 1.3576229705877048e-02;
          }
            break;
          default:
            throw std::runtime_error("Wrong initialization of 1D Gauss quadrature. Order not supported.");
        }
      }
  };
}

#endif // __Quadrature_Gauss1D_H
