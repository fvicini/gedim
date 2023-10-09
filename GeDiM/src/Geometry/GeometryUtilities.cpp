#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"
#include "CommonUtilities.hpp"
#include <unordered_set>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // **************************************************************
  GeometryUtilities::GeometryUtilities(const GeometryUtilitiesConfig& configuration) :
    _configuration(configuration)
  {
  }
  GeometryUtilities::~GeometryUtilities()
  {
  }
  // ***************************************************************************
  vector<double> GeometryUtilities::EquispaceCoordinates(const double& step,
                                                         const bool& insertExtremes) const
  {
    Output::Assert(IsValue1DPositive(step) &&
                   Compare1DValues(step, 1.0) != CompareTypes::SecondBeforeFirst);

    return insertExtremes ?
          EquispaceCoordinates(static_cast<unsigned int>(1.0 / step + 0.5) + 1, 0.0, 1.0, true) :
          EquispaceCoordinates(static_cast<unsigned int>(1.0 / step + 0.5) - 1, 0.0, 1.0, false);
  }
  // ***************************************************************************
  std::vector<double> GeometryUtilities::EquispaceCoordinates(const unsigned int& size,
                                                              const double& origin,
                                                              const double& end,
                                                              const bool& insertExtremes) const
  {
    if (size == 0)
      return { };
    else if (size == 1 && insertExtremes)
      throw invalid_argument("size is not valid with false insertExtremes");

    const VectorXd generated = insertExtremes ? VectorXd::LinSpaced(size,
                                                                    origin,
                                                                    end) :
                                                VectorXd::LinSpaced(size + 2,
                                                                    origin,
                                                                    end);

    vector<double> coordinates;
    if (insertExtremes)
      coordinates.resize(generated.size());
    else
      coordinates.resize(generated.size() - 2);

    for (unsigned int c = 0; c < coordinates.size(); c++)
      coordinates[c] = insertExtremes ? generated[c] : generated[c + 1];

    return coordinates;
  }
  // ***************************************************************************
  std::vector<double> GeometryUtilities::RandomCoordinates(const unsigned int size,
                                                           const bool insertExtremes,
                                                           const unsigned int seed) const
  {
    if (size == 0)
      return { };
    else if (size == 1 && insertExtremes)
      throw invalid_argument("size is not valid with false insertExtremes");

    const unsigned int max_value = 1.0 / (3.0 * _configuration.Tolerance);

    Gedim::Output::Assert(size <= max_value);

    const unsigned int num_random_points = insertExtremes ? (size - 2) :
                                                            size;
    std::vector<unsigned int> random_points = num_random_points == 0 ? std::vector<unsigned int>() :
                                                                       Gedim::Utilities::RandomArrayNoRepetition(num_random_points,
                                                                                                                 max_value,
                                                                                                                 seed);
    std::sort(random_points.begin(),
              random_points.end());

    list<double> coordinates;

    if (insertExtremes)
      coordinates.push_back(0.0);

    for (const unsigned int random_point : random_points)
    {
      if (insertExtremes &&
          (random_point == 0 ||
           random_point == max_value))
        continue;

      coordinates.push_back(static_cast<double>(random_point) /
                            static_cast<double>(max_value));
    }

    if (insertExtremes)
      coordinates.push_back(1.0);

    if (coordinates.size() < size)
    {
      const double middlePoint = 0.5 * (*coordinates.begin() +
                                        *std::next(coordinates.begin()));
      coordinates.insert(std::next(coordinates.begin()),
                         middlePoint);
    }

    if (coordinates.size() < size)
    {
      const double middlePoint = 0.5 * (*coordinates.rbegin() +
                                        *std::next(coordinates.rbegin()));
      coordinates.insert(std::prev(coordinates.end()),
                         middlePoint);
    }

    return std::vector<double>(coordinates.begin(),
                               coordinates.end());
  }
  // ***************************************************************************
  GeometryUtilities::CompareTypes GeometryUtilities::CompareValues(const double& first,
                                                                   const double& second,
                                                                   const double& tolerance) const
  {
    const double max_tolerance = std::max(abs(tolerance),
                                          std::numeric_limits<double>::epsilon());

    double relativeValue = (abs(first) <= max_tolerance ||
                            abs(second) <= max_tolerance) ? 1.0 :
                                                            abs(first);
    double difference = second - first;

    if (abs(difference) <= max_tolerance * relativeValue)
      return CompareTypes::Coincident;
    else if (difference < -max_tolerance * relativeValue)
      return CompareTypes::SecondBeforeFirst;
    else
      return CompareTypes::FirstBeforeSecond;
  }
  // ***************************************************************************
  Matrix3d GeometryUtilities::PlaneRotationMatrix(const Eigen::Vector3d& planeNormal) const
  {
    Matrix3d Q;
    Q.setIdentity();

    Output::Assert(Compare1DValues(planeNormal.norm(), 1.0) == CompareTypes::Coincident);

    // if planeNormal is already oriented as z-axis the return the identity
    const double n_xy = sqrt(planeNormal.x() * planeNormal.x() + planeNormal.y() * planeNormal.y());
    if (IsValue1DZero(n_xy))
      return Q;

    Q.col(0)<< -planeNormal.y() / n_xy, planeNormal.x() / n_xy, 0.0;
    Q.col(1)<< -planeNormal.x() * planeNormal.z() / n_xy, -planeNormal.y() * planeNormal.z() / n_xy, n_xy;
    Q.col(2)<< planeNormal;

    return Q;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::ConvexHull(const Eigen::MatrixXd& points,
                                                     const bool& includeCollinear) const
  {
    Output::Assert(points.rows() == 3 && PointsAre2D(points));

    list<unsigned int> convexHull;
    const unsigned int numPoints = points.cols();

    if (numPoints == 1)
      return vector<unsigned int> { 0 };

    struct pt final
    {
        Eigen::Vector3d Point;
        unsigned int id;
    };

    vector<pt> structPoints(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
      structPoints[p] = { points.col(p), p };

    pt p0 = { Eigen::Vector3d(std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max(), 0.0),
              0
            };
    for (const auto& p : structPoints)
    {
      if (IsValue1DGreater(p0.Point.y(), p.Point.y()))
        p0 = p;

      if (Are1DValuesEqual(p0.Point.y(), p.Point.y()) &&
          IsValue1DGreater(p0.Point.x(), p.Point.x()))
        p0 = p;
    }

    sort(structPoints.begin(),
         structPoints.end(),
         [&p0, this](const pt& a, const pt& b)
    {
      const double norm_b_a = (a.Point - b.Point).norm();
      const double squared_norm_b_p0 = (p0.Point - b.Point).squaredNorm();

      const double orientation_polar = PolarAngle(a.Point,
                                                  b.Point,
                                                  p0.Point,
                                                  norm_b_a,
                                                  sqrt(squared_norm_b_p0));

      if (IsValue1DZero(orientation_polar))
      {
        const double squared_norm_a_p0 = (p0.Point - a.Point).squaredNorm();
        return IsValue1DGreater(squared_norm_b_p0,
                                squared_norm_a_p0);
      }

      return IsValue1DNegative(orientation_polar);
    });

    if (includeCollinear)
    {
      int i = (int)structPoints.size() - 1;

      double norm_b_a = (structPoints[i].Point - structPoints.back().Point).norm();
      double norm_b_p0 = (p0.Point - structPoints.back().Point).norm();
      double orientation_polar = PolarAngle(structPoints[i].Point,
                                            structPoints.back().Point,
                                            p0.Point,
                                            norm_b_a,
                                            norm_b_p0);

      while (i >= 0 &&
             IsValue1DZero(orientation_polar))
      {
        i--;

        norm_b_a = (structPoints[i].Point - structPoints.back().Point).norm();
        norm_b_p0 = (p0.Point - structPoints.back().Point).norm();
        orientation_polar = PolarAngle(structPoints[i].Point,
                                       structPoints.back().Point,
                                       p0.Point,
                                       norm_b_a,
                                       norm_b_p0);
      }

      reverse(structPoints.begin() + i + 1,
              structPoints.end());
    }

    vector<pt> st;
    for (unsigned int i = 0; i < structPoints.size(); i++)
    {
      while (st.size() > 1)
      {
        const double norm_b_a = (st.back().Point - structPoints[i].Point).norm();
        const double norm_b_c = (st[st.size() - 2].Point - structPoints[i].Point).norm();
        const double orientation_polar = PolarAngle(st.back().Point,
                                                    structPoints.at(i).Point,
                                                    st.at(st.size() - 2).Point,
                                                    norm_b_a,
                                                    norm_b_c);
        if (IsValue1DNegative(orientation_polar) ||
            (includeCollinear && IsValue1DZero(orientation_polar)))
          break;

        st.pop_back();
      }

      st.push_back(structPoints[i]);
    }

    structPoints = st;

    vector<unsigned int> result(structPoints.size());
    for (unsigned int p = 0; p < structPoints.size(); p++)
      result[result.size() - 1 - p] = structPoints[p].id;

    return result;
  }
  //  // ***************************************************************************
  MatrixXd GeometryUtilities::ExtractPoints(const Eigen::MatrixXd& points,
                                            const vector<unsigned int>& filter) const
  {
    Eigen::MatrixXd extraction(3, filter.size());
    for (unsigned int c = 0; c < filter.size(); c++)
      extraction.col(c) = points.col(filter[c]);

    return extraction;
  }
  // ***************************************************************************
  vector<Matrix3d> GeometryUtilities::ExtractTriangulationPoints(const Eigen::MatrixXd& points,
                                                                 const vector<unsigned int>& pointsTriangulation) const
  {
    const unsigned int numTriangles = pointsTriangulation.size() / 3;
    vector<Matrix3d> triangulations(numTriangles);

    for (unsigned int t = 0; t < numTriangles; t++)
    {
      Eigen::Matrix3d& triangleVertices = triangulations[t];
      triangleVertices.col(0)<< points.col(pointsTriangulation[3 * t]);
      triangleVertices.col(1)<< points.col(pointsTriangulation[3 * t + 1]);
      triangleVertices.col(2)<< points.col(pointsTriangulation[3 * t + 2]);
    }

    return triangulations;
  }
  // ***************************************************************************
  vector<Matrix3d> GeometryUtilities::ExtractTriangulationPointsByInternalPoint(const Eigen::MatrixXd& points,
                                                                                const Eigen::Vector3d& internalPoint,
                                                                                const vector<unsigned int>& pointsTriangulation) const
  {
    const unsigned int numTriangles = pointsTriangulation.size() / 3;
    vector<Matrix3d> triangulations(numTriangles);

    for (unsigned int t = 0; t < numTriangles; t++)
    {
      Eigen::Matrix3d& triangleVertices = triangulations[t];
      triangleVertices.col(0)<< internalPoint;
      triangleVertices.col(1)<< points.col(pointsTriangulation[3 * t + 1]);
      triangleVertices.col(2)<< points.col(pointsTriangulation[3 * t + 2]);
    }

    return triangulations;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::CreateEllipse(const double& axisMajorLength,
                                            const double& axisMinorLength,
                                            const unsigned int& resolution) const
  {
    const vector<double> ellipseXPoints = EquispaceCoordinates(resolution,
                                                               0.0,
                                                               axisMajorLength,
                                                               false);
    const unsigned int numInternalPoints = ellipseXPoints.size();
    vector<double> ellipseYPoints(numInternalPoints, 0.0);

    for (unsigned int p = 0; p < numInternalPoints; p++)
    {
      ellipseYPoints[p] =  axisMinorLength * sqrt(1.0 -
                                                  ellipseXPoints.at(p) *
                                                  ellipseXPoints.at(p) /
                                                  (axisMajorLength *
                                                   axisMajorLength));
    }

    Eigen::MatrixXd vertices(3, 4 + 4 * numInternalPoints);

    vertices.col(0)<< Vector3d(axisMajorLength, 0.0, 0.0);
    vertices.col(numInternalPoints + 1)<< Vector3d(0.0, axisMinorLength, 0.0);
    vertices.col(2 * (numInternalPoints + 1))<< Vector3d(-axisMajorLength, 0.0, 0.0);
    vertices.col(3 * (numInternalPoints + 1))<< Vector3d(0.0, -axisMinorLength, 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(v + 1)<< Vector3d(ellipseXPoints[numInternalPoints - 1 - v], ellipseYPoints[numInternalPoints - 1 - v], 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(numInternalPoints + 1 + v + 1)<< Vector3d(-ellipseXPoints[v], ellipseYPoints[v], 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(2 * (numInternalPoints + 1) + v + 1)<< Vector3d(-ellipseXPoints[numInternalPoints - 1 - v], -ellipseYPoints[numInternalPoints - 1 - v], 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(3 * (numInternalPoints + 1) + v + 1)<< Vector3d(ellipseXPoints[v], -ellipseYPoints[v], 0.0);

    return vertices;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::CreateTriangle(const Eigen::Vector3d& p1,
                                             const Eigen::Vector3d& p2,
                                             const Eigen::Vector3d& p3) const
  {
    MatrixXd vertices(3, 3);
    vertices.col(0)<< p1;
    vertices.col(1)<< p2;
    vertices.col(2)<< p3;

    return vertices;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::CreateParallelogram(const Eigen::Vector3d& origin,
                                                  const Eigen::Vector3d& lengthVector,
                                                  const Eigen::Vector3d& widthVector) const
  {
    MatrixXd vertices(3, 4);
    vertices.col(0)<< origin;
    vertices.col(1)<< origin + lengthVector;
    vertices.col(2)<< origin + lengthVector + widthVector;
    vertices.col(3)<< origin + widthVector;
    return vertices;
  }
  // ***************************************************************************
}
