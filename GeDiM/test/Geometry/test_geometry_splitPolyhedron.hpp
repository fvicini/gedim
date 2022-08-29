#ifndef __TEST_GEOMETRY_SPLITPOLYHEDRON_H
#define __TEST_GEOMETRY_SPLITPOLYHEDRON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // no action
      {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                             Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                             Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                                                             Eigen::Vector3d(0.0, 1.0, 0.0));
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtility.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                              polyhedron.Edges);
        const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtility.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                      polyhedron.Faces);
        const vector<vector<bool>> polyhedronFaceEdgeDirections = geometryUtility.PolyhedronFaceEdgeDirections(polyhedron.Vertices,
                                                                                                               polyhedron.Edges,
                                                                                                               polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents = geometryUtility.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                                                                               polyhedron.Edges,
                                                                                                               polyhedron.Faces,
                                                                                                               polyhedronFaceEdgeDirections,
                                                                                                               polyhedronEdgeTangents);

        const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.5);
        const Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
        const Eigen::MatrixXd planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtility.PlaneTranslation(planeNormal,
                                                                                  planeOrigin);

        {
          using namespace Gedim;
          cerr<< "polyhedron.Vertices:\n"<< polyhedron.Vertices<< endl;
          cerr<< "polyhedron.Edges:\n"<< polyhedron.Edges<< endl;
          cerr<< "polyhedron.Faces:\n"<< polyhedron.Faces<< endl;
        }

        geometryUtility.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                 polyhedron.Edges,
                                                 polyhedron.Faces,
                                                 polyhedronFaceVertices,
                                                 polyhedronFaceEdgesTangents,
                                                 planeNormal,
                                                 planeOrigin,
                                                 planeRotationMatrix,
                                                 planeTranslation);

        //ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::NoAction);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_SPLIT_H
