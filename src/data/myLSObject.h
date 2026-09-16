/**
 * @file myLSObject.h
 * @brief Declares level-set input objects and mesh-to-grid construction utilities.
 */
#pragma once

#include "CudaTypes.h"
#include "math/Vector3.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace fundem::levelset
{

/**
 * Abstract input geometry that generates a signed-distance grid and triangulated
 * surface representation consumable by `LSGeometryContainer`.
 */
class LSInfo
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    LSInfo() = default;
    LSInfo(const LSInfo&) = default;
    LSInfo& operator=(const LSInfo&) = default;
    LSInfo(LSInfo&&) noexcept = default;
    LSInfo& operator=(LSInfo&&) noexcept = default;
    virtual ~LSInfo() = default;

    /** Reports whether the concrete shape has enough valid parameters to build a grid. */
    bool valid() const noexcept { return isValid(); }
    /** Builds a signed-distance grid from a target number of cells across the diameter. */
    void buildLSGrid(int resolutionPerDiameter = 50);
    /** Builds a signed-distance grid with explicit spacing and padding cells. */
    void buildLSGrid(Real spacing, int paddingSize = 2);
    /** Reverses the sign convention of every stored signed-distance value. */
    void reverseSDFSign() noexcept;
    /** Writes the generated Cartesian signed-distance grid as a VTI file. */
    void outputGridVTI(const std::string& fileName) const;

    /** Interpolates the signed distance at a point in the geometry-local frame. */
    Real signedDistance(const Vec3& point) const;
    /** Returns a representative bounding radius used to size the grid. */
    virtual Real radius() const noexcept;
    /** Returns the average tributary area of one generated surface node. */
    virtual Real meanSurfaceNodeArea() const noexcept;

    const Vec3& gridNodeOrigin() const noexcept { return gridNodeOrigin_; }
    const int3& gridNodeSize3D() const noexcept { return gridNodeSize_; }
    Real gridNodeSpacing() const noexcept { return gridNodeSpacing_; }
    const std::vector<Real>& gridNodeSFD() const noexcept { return gridNodeSignedDistance_; }
    const std::vector<Vec3>& surfaceNodePosition() const noexcept { return surfaceNodePositions_; }
    const std::vector<int3>& surfaceNodeConnectivity() const noexcept { return surfaceTriangles_; }

protected:
    /** Clears generated grid data while retaining analytic shape parameters. */
    void clearGrid() noexcept;
    /** Builds and projects an icosphere sampling of an implicit surface. */
    void buildImplicitSurfaceNode(int subdivisionLevel);
    /** Subdivides each stored surface triangle once. */
    void subdivideSurface();
    /** Reverses inconsistent triangles so normals point out of the solid. */
    void orientSurfaceOutward();

    std::vector<Vec3> surfaceNodePositions_; ///< Geometry-local surface-node positions.
    std::vector<int3> surfaceTriangles_;     ///< Triangle connectivity into `surfaceNodePositions_`.
    bool configured_{false};                 ///< Whether shape parameters define a valid object.

private:
    /** Validates the concrete shape parameters. */
    virtual bool isValid() const noexcept = 0;
    /** Returns the lower local-space extent used to build the Cartesian grid. */
    virtual Vec3 boundingBoxMin() const noexcept = 0;
    /** Returns the upper local-space extent used to build the Cartesian grid. */
    virtual Vec3 boundingBoxMax() const noexcept = 0;
    /** Evaluates analytic or mesh signed distance in the local frame. */
    virtual Real evaluateSFD(const Vec3& point) const noexcept = 0;
    /** Projects a sampling guess onto the represented surface. */
    virtual Vec3 projectToSurface(const Vec3& point) const noexcept { return point; }
    /** Returns a diagnostic name used in validation errors. */
    virtual const char* objectName() const noexcept = 0;

    /** Creates the base icosahedron used by implicit-surface sampling. */
    void buildIcosahedron();
    /** Samples the current shape and extracts node-area data. */
    void buildLSGridKernel(int paddingSize);

    Vec3 gridNodeOrigin_{Vec3::zero()};        ///< Position of grid node (0,0,0) in the local frame.
    int3 gridNodeSize_{0, 0, 0};               ///< Number of nodes along each Cartesian direction.
    Real gridNodeSpacing_{0.0};                ///< Uniform grid spacing.
    std::vector<Real> gridNodeSignedDistance_; ///< Flattened signed-distance samples.
};

using LSObject = LSInfo;

/** Triangle-soup geometry loaded from memory or an OBJ file. */
class TriangleMesh final : public LSInfo
{
public:
    TriangleMesh();
    /** Creates a triangle mesh from indexed geometry and builds acceleration data. */
    TriangleMesh(const std::vector<Vec3>& vertexPositions, const std::vector<int3>& triangles);
    ~TriangleMesh() override;
    TriangleMesh(const TriangleMesh&) = delete;
    TriangleMesh& operator=(const TriangleMesh&) = delete;
    TriangleMesh(TriangleMesh&&) noexcept;
    TriangleMesh& operator=(TriangleMesh&&) noexcept;

    /** Replaces the mesh with triangles loaded from an OBJ file. */
    void loadOBJ(const std::string& fileName);
    /** Subdivides every surface triangle once and rebuilds acceleration data. */
    void fineMesh();

private:
    struct accelerationData;

    bool isValid() const noexcept override;
    Vec3 boundingBoxMin() const noexcept override;
    Vec3 boundingBoxMax() const noexcept override;
    Real evaluateSFD(const Vec3& point) const noexcept override;
    const char* objectName() const noexcept override { return "Triangle Mesh"; }
    /** Validates connectivity, orients triangles, and rebuilds acceleration data. */
    void initializeMesh();

    std::unique_ptr<accelerationData> acceleration_; ///< Owned mesh-query acceleration structure.
};

/**
 * Creates a closed, outward-facing mesh by smoothly deforming an icosphere.
 * Heights are radial offsets from the positive base radius; both must be finite,
 * minimumSurfaceHeight <= maximumSurfaceHeight, and radius + minimumSurfaceHeight > 0.
 * Sampled vertex radii span [radius + minimumSurfaceHeight, radius + maximumSurfaceHeight].
 * The result is star-shaped about the origin, not a generator of arbitrary topology.
 * A fixed seed reproduces the mesh; equal heights produce a spherical mesh.
 * subdivisionLevel is non-negative and controls surface resolution only. The
 * TriangleMesh provides mesh signed distances; callers build the LS grid separately.
 */
std::unique_ptr<TriangleMesh> makeRandomShape(math::Real radius,
                                            math::Real minimumSurfaceHeight,
                                            math::Real maximumSurfaceHeight,
                                            int subdivisionLevel = 3,
                                            std::uint64_t seed = 0);

/** Analytic sphere input geometry centered at the local origin. */
class Sphere final : public LSInfo
{
public:
    Sphere() = default;
    /** Creates an analytic sphere with positive @p radius when valid. */
    explicit Sphere(Real radius);

    void setParameter(Real radius);
    /** Creates an icosphere surface with the requested subdivision level. */
    void buildSurfaceNode(int subdivisionLevel = 4);
    Real radius() const noexcept override { return radius_; }

private:
    bool isValid() const noexcept override;
    Vec3 boundingBoxMin() const noexcept override;
    Vec3 boundingBoxMax() const noexcept override;
    Real evaluateSFD(const Vec3& point) const noexcept override;
    Vec3 projectToSurface(const Vec3& point) const noexcept override;
    const char* objectName() const noexcept override { return "Sphere"; }

    Real radius_{0.0}; ///< Sphere radius.
};

/** Analytic superellipsoid with independent semi-axes and shape exponents. */
class Superellipsoid final : public LSInfo
{
public:
    Superellipsoid() = default;
    /** Creates a superellipsoid from semi-axes and equatorial/polar exponents. */
    Superellipsoid(Real radiusX, Real radiusY, Real radiusZ, Real equatorialExponent, Real polarExponent);

    void setParameter(Real radiusX, Real radiusY, Real radiusZ, Real equatorialExponent, Real polarExponent);
    /** Creates a projected icosphere sampling of the surface. */
    void buildSurfaceNode(int subdivisionLevel = 4);
    Real radius() const noexcept override;

private:
    bool isValid() const noexcept override;
    Vec3 boundingBoxMin() const noexcept override;
    Vec3 boundingBoxMax() const noexcept override;
    Real evaluateSFD(const Vec3& point) const noexcept override;
    Vec3 projectToSurface(const Vec3& point) const noexcept override;
    /** Evaluates the local implicit-function gradient used for projection. */
    Vec3 evaluateImplicitGradient(const Vec3& point) const noexcept;
    const char* objectName() const noexcept override { return "Superellipsoid"; }

    Real radiusX_{0.0};            ///< Semi-axis length along x.
    Real radiusY_{0.0};            ///< Semi-axis length along y.
    Real radiusZ_{0.0};            ///< Semi-axis length along z.
    Real equatorialExponent_{0.0}; ///< Equatorial shape exponent.
    Real polarExponent_{0.0};      ///< Polar shape exponent.
};

/** Finite square plane representing an infinite-mass wall. */
class PlaneWall final : public LSInfo
{
public:
    /** Creates a finite square wall oriented by @p outwardNormal. */
    PlaneWall(const Vec3& outwardNormal, Real size);
    void setParameter(const Vec3& outwardNormal, Real size);

private:
    bool isValid() const noexcept override;
    Vec3 boundingBoxMin() const noexcept override;
    Vec3 boundingBoxMax() const noexcept override;
    Real evaluateSFD(const Vec3& point) const noexcept override;
    const char* objectName() const noexcept override { return "Plane Wall"; }
    /** Generates the square's surface nodes and triangle connectivity. */
    void buildSurfaceNode();

    Vec3 outwardNormal_{Vec3::unitZ()}; ///< Unit normal pointing out of the solid.
    Real size_{0.0};                    ///< Side length of the sampled square surface.
};

/** Axis-aligned rectangular wall represented by its six faces. */
class BoxWall : public LSInfo
{
public:
    /** Creates an axis-aligned box wall with full dimensions @p size. */
    explicit BoxWall(const Vec3& size);
    void setParameter(const Vec3& size);
    const Vec3& boxSize3D() const noexcept { return size_; }

protected:
    const char* objectName() const noexcept override { return "Box Wall"; }
    /** Generates the closed six-face box surface. */
    void buildSurfaceNode();

private:
    bool isValid() const noexcept override;
    Vec3 boundingBoxMin() const noexcept override;
    Vec3 boundingBoxMax() const noexcept override;
    Real evaluateSFD(const Vec3& point) const noexcept override;

    Vec3 size_{Vec3::zero()}; ///< Full box dimensions.
};

/** Closed box particle with surface nodes sampled at a requested spacing. */
class BoxParticle final : public BoxWall
{
public:
    using BoxWall::BoxWall;
    /** Resamples the closed box surface at approximately @p surfaceSpacing. */
    void buildSurfaceNode(Real surfaceSpacing);
    Real meanSurfaceNodeArea() const noexcept override;

private:
    const char* objectName() const noexcept override { return "Box Particle"; }
};

/** Open-ended cylindrical wall aligned between two local-frame endpoints. */
class CylinderWall final : public LSInfo
{
public:
    /** Creates an open cylinder between the two rim centers. */
    CylinderWall(const Vec3& bottomCenter, const Vec3& topCenter, Real radius);
    void setParameter(const Vec3& bottomCenter, const Vec3& topCenter, Real radius);
    void setCircumferentialSegments(int segmentCount);

private:
    bool isValid() const noexcept override;
    Vec3 boundingBoxMin() const noexcept override;
    Vec3 boundingBoxMax() const noexcept override;
    Real evaluateSFD(const Vec3& point) const noexcept override;
    const char* objectName() const noexcept override { return "Cylinder Wall"; }
    /** Generates the lateral wall with @p segmentCount circumferential segments. */
    void buildSurfaceNode(int segmentCount);

    Vec3 bottomCenter_{0.0, 0.0, -0.5}; ///< Center of the bottom rim.
    Vec3 topCenter_{0.0, 0.0, 0.5};     ///< Center of the top rim.
    Real radius_{0.0};                  ///< Cylinder radius.
    int segmentCount_{360};             ///< Circumferential surface resolution.
};

/** Open-ended conical-frustum wall aligned between two local-frame endpoints. */
class ConeWall final : public LSInfo
{
public:
    /** Creates an open conical frustum between the two rim centers. */
    ConeWall(const Vec3& bottomCenter, const Vec3& topCenter, Real bottomRadius, Real topRadius);
    void setParameter(const Vec3& bottomCenter, const Vec3& topCenter, Real bottomRadius, Real topRadius);
    void setCircumferentialSegments(int segmentCount);

private:
    bool isValid() const noexcept override;
    Vec3 boundingBoxMin() const noexcept override;
    Vec3 boundingBoxMax() const noexcept override;
    Real evaluateSFD(const Vec3& point) const noexcept override;
    const char* objectName() const noexcept override { return "Cone Wall"; }
    /** Generates the lateral frustum with @p segmentCount circumferential segments. */
    void buildSurfaceNode(int segmentCount);

    Vec3 bottomCenter_{0.0, 0.0, -0.5}; ///< Center of the bottom rim.
    Vec3 topCenter_{0.0, 0.0, 0.5};     ///< Center of the top rim.
    Real bottomRadius_{0.0};            ///< Radius of the bottom rim.
    Real topRadius_{0.0};               ///< Radius of the top rim.
    int segmentCount_{360};             ///< Circumferential surface resolution.
};

} // namespace fundem::levelset
