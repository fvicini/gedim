#ifndef __Quadrature_GaussLobatto1D_H
#define __Quadrature_GaussLobatto1D_H

#include "Eigen/Eigen"

namespace Gedim
{
  /// Gauss quadrature rule for segments
  class Quadrature_GaussLobatto1D
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
            points.setZero(3,2);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 1.0000000000000000e+00;
            weights.resize(2);
            weights[0] = 5.0000000000000000e-01;
            weights[1] = 5.0000000000000000e-01;
          }
            break;
          case 2:
          case 3:
          {
            points.setZero(3,3);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 5.0000000000000000e-01;
            points(0,2) = 1.0000000000000000e+00;
            weights.resize(3);
            weights[0] = 1.6666666666666649e-01;
            weights[1] = 6.6666666666666652e-01;
            weights[2] = 1.6666666666666649e-01;
          }
            break;
          case 4:
          case 5:
          {
            points.setZero(3,4);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 2.7639320225002106e-01;
            points(0,2) = 7.2360679774997894e-01;
            points(0,3) = 1.0000000000000000e+00;
            weights.resize(4);
            weights[0] = 8.3333333333333329e-02;
            weights[1] = 4.1666666666666663e-01;
            weights[2] = 4.1666666666666663e-01;
            weights[3] = 8.3333333333333329e-02;
          }
            break;
          case 6:
          case 7:
          {
            points.setZero(3,5);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 1.7267316464601151e-01;
            points(0,2) = 5.0000000000000000e-01;
            points(0,3) = 8.2732683535398843e-01;
            points(0,4) = 1.0000000000000000e+00;
            weights.resize(5);
            weights[0] = 5.0000000000000003e-02;
            weights[1] = 2.7222222222222220e-01;
            weights[2] = 3.5555555555555557e-01;
            weights[3] = 2.7222222222222220e-01;
            weights[4] = 5.0000000000000003e-02;
          }
            break;
          case 8:
          case 9:
          {
            points.setZero(3,6);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 1.1747233803526763e-01;
            points(0,2) = 3.5738424175967742e-01;
            points(0,3) = 6.4261575824032258e-01;
            points(0,4) = 8.8252766196473242e-01;
            points(0,5) = 1.0000000000000000e+00;
            weights.resize(6);
            weights[0] = 3.3333333333333333e-02;
            weights[1] = 1.8923747814892350e-01;
            weights[2] = 2.7742918851774317e-01;
            weights[3] = 2.7742918851774317e-01;
            weights[4] = 1.8923747814892350e-01;
            weights[5] = 3.3333333333333333e-02;
          }
            break;
          case 10:
          case 11:
          {
            points.setZero(3,7);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 8.4888051860716518e-02;
            points(0,2) = 2.6557560326464291e-01;
            points(0,3) = 5.0000000000000000e-01;
            points(0,4) = 7.3442439673535709e-01;
            points(0,5) = 9.1511194813928354e-01;
            points(0,6) = 1.0000000000000000e+00;
            weights.resize(7);
            weights[0] = 2.3809523809523808e-02;
            weights[1] = 1.3841302368078298e-01;
            weights[2] = 2.1587269060493131e-01;
            weights[3] = 2.4380952380952381e-01;
            weights[4] = 2.1587269060493131e-01;
            weights[5] = 1.3841302368078298e-01;
            weights[6] = 2.3809523809523808e-02;
          }
            break;
          case 12:
          case 13:
          {
            points.setZero(3,8);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 6.4129925745196714e-02;
            points(0,2) = 2.0414990928342885e-01;
            points(0,3) = 3.9535039104876057e-01;
            points(0,4) = 6.0464960895123943e-01;
            points(0,5) = 7.9585009071657109e-01;
            points(0,6) = 9.3587007425480329e-01;
            points(0,7) = 1.0000000000000000e+00;
            weights.resize(8);
            weights[0] = 1.7857142857142856e-02;
            weights[1] = 1.0535211357175302e-01;
            weights[2] = 1.7056134624175218e-01;
            weights[3] = 2.0622939732935194e-01;
            weights[4] = 2.0622939732935194e-01;
            weights[5] = 1.7056134624175218e-01;
            weights[6] = 1.0535211357175302e-01;
            weights[7] = 1.7857142857142856e-02;
          }
            break;
          case 14:
          case 15:
          {
            points.setZero(3,9);
            points(0,0) = 0.0000000000000000e+00;
            points(0,1) = 5.0121002294269912e-02;
            points(0,2) = 1.6140686024463113e-01;
            points(0,3) = 3.1844126808691092e-01;
            points(0,4) = 5.0000000000000000e-01;
            points(0,5) = 6.8155873191308913e-01;
            points(0,6) = 8.3859313975536887e-01;
            points(0,7) = 9.4987899770573003e-01;
            points(0,8) = 1.0000000000000000e+00;
            weights.resize(9);
            weights[0] = 1.3888888888888888e-02;
            weights[1] = 8.2747680780402760e-02;
            weights[2] = 1.3726935625008085e-01;
            weights[3] = 1.7321425548652317e-01;
            weights[4] = 1.8575963718820862e-01;
            weights[5] = 1.7321425548652317e-01;
            weights[6] = 1.3726935625008085e-01;
            weights[7] = 8.2747680780402760e-02;
            weights[8] = 1.3888888888888888e-02;
          }
            break;
          default:
            throw std::runtime_error("Wrong initialization of 1D Gauss Lobatto quadrature. Order not supported.");
        }
      }
  };
}

#endif // __Quadrature_GaussLobatto1D_H
