/** @file tutorial1.cpp @brief A tilted level-set pyramid falling onto a plane. */
#include "solver/solvers.h"

int main()
{
    using namespace fundem;
    using math::Vec3;

    // 1. Pure LSDEM supports CPU or GPU; select the mode before initialization.
    LSDEM simulation;
    simulation.setExecutionMode(executionMode::CPU);
    const int materialIndex = simulation.addMaterial(material{1.0e7, 0.25, 0.6, 2000.0});

    // 2. Build surface nodes before the SDF grid. Grid construction computes
    // mass properties and shifts the geometry into its centroid frame.
    levelset::Pyramid pyramid{0.04, 0.03, 0.05};
    pyramid.buildSurfaceNode(0.01);
    pyramid.buildLSGrid(24);
    const int pyramidIndex = simulation.addGeometry(pyramid);

    levelset::PlaneWall plane{Vec3::unitZ(), 0.4};
    plane.buildLSGrid(0.02, 2);
    const int planeIndex = simulation.addGeometry(plane);

    // 3. Assign material/geometry before setting the particle's centroid pose.
    LSParticle body;
    body.setMaterial(simulation.materials(), materialIndex);
    body.setGeometry(simulation.geometries(), pyramidIndex);
    body.setPosition({0.0, 0.0, 0.06});
    body.setOrientation(math::Quaternion::fromUnitAxisAngle(Vec3::unitY(), math::degreesToRadians(20.0)));
    simulation.addLSParticle(body);

    LSParticle wall;
    wall.setMaterial(simulation.materials(), materialIndex);
    wall.setGeometry(simulation.geometries(), planeIndex);
    simulation.addLSParticle(wall); // The zero-volume plane is fixed.

    // 4. LS contact requires normal/shear stiffness per nodal area, in N/m^3.
    materialPair contact{materialIndex, materialIndex};
    contact.setLevelSetStiffness(1.0e8, 5.0e7);
    contact.setFrictionCoefficients(0.3);
    simulation.setMaterialPair(contact);

    // 5. Set the search region and run settings directly; rebuild after editing.
    simulation.setBoundary({-0.25, -0.25, -0.25}, {0.25, 0.25, 0.25});
    simulation.setGravity({0.0, 0.0, -9.81});
    simulation.setTimeStep(1.0e-4);
    simulation.setOutputDirectory("tutorial1_files");
    simulation.setOutputStepInterval(100);
    simulation.setVTUOutputEnabled(true);
    simulation.setVTUOutputFormat(vtuFormat::binary);
    simulation.solve(3000); // 0.3 s, 31 frames including the initial state.
    return 0;
}
