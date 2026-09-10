/** Verifies one bond in axial, shear, bending, and torsional deformation modes. */
#include "interaction/bond.h"
#include "material/material.h"
#include "particle/particle.h"
#include "validationCommon.h"

#include <algorithm>
#include <exception>

namespace
{

using namespace fundem;
using Real = math::Real;
using Vec3 = math::Vec3;
using Quaternion = math::Quaternion;

constexpr Real equivalentLength = 1.0;
constexpr Real normalStiffness = 1.0e5;
constexpr Real shearStiffness = 2.0e4;
constexpr Real bendingStiffness = 3.0e3;
constexpr Real torsionalStiffness = 4.0e3;

/** Resultant loads and modal elastic energies from one bond evaluation. */
struct response {
    Vec3 force_{Vec3::zero()};        ///< Force acting on the master endpoint.
    Vec3 masterTorque_{Vec3::zero()}; ///< Torque acting on the master body.
    Vec3 slaveTorque_{Vec3::zero()};  ///< Torque acting on the slave body.
    Real normalEnergy_{0.0};          ///< Axial elastic energy.
    Real shearEnergy_{0.0};           ///< Shear elastic energy.
    Real bendingEnergy_{0.0};         ///< Bending elastic energy.
    Real torsionalEnergy_{0.0};       ///< Torsional elastic energy.
};

/** Builds a fresh reference bond, applies endpoint states, and evaluates it. */
response evaluate(const Vec3& masterDisplacement, const Quaternion& masterOrientation, const Vec3& slaveDisplacement, const Quaternion& slaveOrientation)
{
    materialContainer materials;
    materials.host().push_back(material{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1000.0});

    particleContainer particles;
    particle master;
    master.setPosition(Vec3::zero());
    master.setRadius(0.1);
    master.setMaterial(materials, 0);
    particles.host().push_back(master);

    particle slave;
    slave.setPosition({equivalentLength, 0.0, 0.0});
    slave.setRadius(0.1);
    slave.setMaterial(materials, 0);
    particles.host().push_back(slave);

    bond connection{equivalentLength};
    validation::require(connection.setConnection(particles, 0, 1, -Vec3::unitX()), "Cannot set the bond reference connection.");
    validation::require(connection.setStiffness(normalStiffness, shearStiffness, bendingStiffness, torsionalStiffness), "Cannot set the bond stiffness.");

    bondContainer bonds;
    bonds.host().push_back(connection);
    particles.host()[0].setPosition(masterDisplacement);
    particles.host()[0].setOrientation(masterOrientation);
    particles.host()[1].setPosition(Vec3{equivalentLength, 0.0, 0.0} + slaveDisplacement);
    particles.host()[1].setOrientation(slaveOrientation);

    bond& storedBond = bonds.host().front();
    validation::require(storedBond.updateMasterSlaveParticle(particles), "Cannot refresh bond endpoint particles.");
    validation::require(storedBond.calculateForce(), "Cannot calculate the bond response.");
    return {storedBond.force(),
            storedBond.masterTorque(),
            storedBond.slaveTorque(),
            storedBond.normalElasticEnergy(),
            storedBond.shearElasticEnergy(),
            storedBond.bendingElasticEnergy(),
            storedBond.torsionalElasticEnergy()};
}

/** Compares signed numerical and reference values by magnitude. */
Real normalizedMagnitudeError(Real value, Real reference) noexcept { return std::abs(std::abs(value) - std::abs(reference)) / std::abs(reference); }

} // namespace

int main(int argc, char** argv)
{
    try
    {
        constexpr Real displacement = 1.0e-4;
        constexpr Real angle = 1.0e-4;
        const Quaternion identity = Quaternion::identity();

        const response axial = evaluate(Vec3::zero(), identity, displacement * Vec3::unitX(), identity);
        const response shear = evaluate(Vec3::zero(), identity, displacement * Vec3::unitY(), identity);
        const response bending = evaluate(Vec3::zero(), Quaternion::fromUnitAxisAngle(Vec3::unitZ(), -0.5 * angle), Vec3::zero(), Quaternion::fromUnitAxisAngle(Vec3::unitZ(), 0.5 * angle));
        const response torsion = evaluate(Vec3::zero(), Quaternion::fromUnitAxisAngle(Vec3::unitX(), -0.5 * angle), Vec3::zero(), Quaternion::fromUnitAxisAngle(Vec3::unitX(), 0.5 * angle));

        const Real axialError = normalizedMagnitudeError(axial.force_.x, normalStiffness * displacement);
        const Real axialEnergyError = normalizedMagnitudeError(axial.normalEnergy_, 0.5 * normalStiffness * displacement * displacement);
        const Real shearError = normalizedMagnitudeError(shear.force_.y, shearStiffness * displacement);
        const Real bendingError = normalizedMagnitudeError(bending.masterTorque_.z, bendingStiffness * angle);
        const Real torsionalError = normalizedMagnitudeError(torsion.masterTorque_.x, torsionalStiffness * angle);
        const Real maximumError = std::max({axialError, axialEnergyError, shearError, bendingError, torsionalError});
        constexpr Real tolerance = 2.0e-4;

        auto output = validation::openDataFile(validation::outputDirectory(argc, argv, "bondModes"), "bondModes.dat");
        output << "# mode numerical expected relativeError elasticEnergy\n";
        output << "axialForce " << axial.force_.x << ' ' << normalStiffness * displacement << ' ' << axialError << ' ' << axial.normalEnergy_ << '\n';
        output << "shearForce " << shear.force_.y << ' ' << shearStiffness * displacement << ' ' << shearError << ' ' << shear.shearEnergy_ << '\n';
        output << "bendingTorque " << bending.masterTorque_.z << ' ' << bendingStiffness * angle << ' ' << bendingError << ' ' << bending.bendingEnergy_ << '\n';
        output << "torsionalTorque " << torsion.masterTorque_.x << ' ' << torsionalStiffness * angle << ' ' << torsionalError << ' ' << torsion.torsionalEnergy_ << '\n';

        validation::require(maximumError <= tolerance, "A single-bond deformation mode exceeded tolerance.");
        validation::printPass("bondModes", maximumError, tolerance);
    }
    catch (const std::exception& error)
    {
        std::cerr << "bondModes: FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
