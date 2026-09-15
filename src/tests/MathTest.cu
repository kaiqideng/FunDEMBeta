#include "math/Quaternion.h"
#include "TestSupport.h"

#include <cuda_runtime.h>

#include <limits>

namespace
{

using fundem::math::Mat3;
using fundem::math::Quaternion;
using fundem::math::Real;
using fundem::math::Vec3;

struct Evaluation {
    Vec3 normalized;
    Mat3 inverse;
    Quaternion orientation;
    Quaternion orientationInverse;
    Vec3 rotated;
    Vec3 recovered;
    Mat3 rotationMatrix;
    Quaternion integrated;
    int successMask{0};
};

FUNDEM_MATH_HD Evaluation evaluate() noexcept
{
    using namespace fundem::math;

    Evaluation result;

    Vec3 large{1.0e300, -1.0e300, 0.0};
    if (tryNormalize(large))
    {
        result.normalized = large;
        result.successMask |= 1;
    }

    const Mat3 matrix{4.0, 1.0, 2.0, 0.0, 3.0, -1.0, 0.0, 0.0, 2.0};
    if (tryInverse(matrix, result.inverse))
    {
        result.successMask |= 2;
    }

    if (Quaternion::tryFromAxisAngle(Vec3::unitZ(), halfPi, result.orientation))
    {
        result.successMask |= 4;
    }

    if ((result.successMask & 4) != 0 && tryInverse(result.orientation, result.orientationInverse))
    {
        result.successMask |= 8;
    }

    if ((result.successMask & 4) != 0 && tryRotate(result.orientation, Vec3::unitX(), result.rotated))
    {
        result.successMask |= 16;
    }

    if ((result.successMask & 4) != 0 && tryInverseRotate(result.orientation, result.rotated, result.recovered))
    {
        result.successMask |= 32;
    }

    if ((result.successMask & 4) != 0 && tryRotationMatrix(result.orientation, result.rotationMatrix))
    {
        result.successMask |= 64;
    }

    if ((result.successMask & 4) != 0 && tryIntegrateWorldAngularVelocity(result.orientation, {0.1, 0.2, 0.3}, 0.01, result.integrated))
    {
        result.successMask |= 128;
    }

    return result;
}

__global__ void evaluateKernel(Evaluation* result) { *result = evaluate(); }

void verify(const Evaluation& result)
{
    using namespace fundem::math;

    FUNDEM_TEST_REQUIRE(result.successMask == 255);
    FUNDEM_TEST_REQUIRE(nearlyEqual(norm(result.normalized), 1.0));
    FUNDEM_TEST_REQUIRE(nearlyEqual(result.normalized.x, 0.70710678118654752440));
    FUNDEM_TEST_REQUIRE(nearlyEqual(result.normalized.y, -0.70710678118654752440));

    const Mat3 matrix{4.0, 1.0, 2.0, 0.0, 3.0, -1.0, 0.0, 0.0, 2.0};
    FUNDEM_TEST_REQUIRE(nearlyEqual(matrix * result.inverse, Mat3::identity()));

    FUNDEM_TEST_REQUIRE(nearlyEqual(result.rotated, Vec3::unitY()));
    FUNDEM_TEST_REQUIRE(nearlyEqual(result.recovered, Vec3::unitX()));
    FUNDEM_TEST_REQUIRE(nearlyEqual(result.rotationMatrix * Vec3::unitX(), result.rotated));
    FUNDEM_TEST_REQUIRE(representsSameRotation(result.orientation * result.orientationInverse, Quaternion::identity()));
    FUNDEM_TEST_REQUIRE(nearlyEqual(norm(result.integrated), 1.0));
}

} // namespace

int main()
{
    using namespace fundem::math;

    constexpr Vec3 constexprCross = cross(Vec3::unitX(), Vec3::unitY());
    static_assert(constexprCross == Vec3::unitZ());
    static_assert(Mat3::identity() * Vec3{1.0, 2.0, 3.0} == Vec3{1.0, 2.0, 3.0});
    static_assert(degreesToRadians(180.0) == pi);

    const Evaluation hostResult = evaluate();
    verify(hostResult);

    Vec3 zeroVector = Vec3::zero();
    FUNDEM_TEST_REQUIRE(!tryNormalize(zeroVector));
    FUNDEM_TEST_REQUIRE(zeroVector == Vec3::zero());

    Mat3 unchangedMatrix = Mat3::identity();
    FUNDEM_TEST_REQUIRE(!tryInverse(Mat3::zero(), unchangedMatrix));
    FUNDEM_TEST_REQUIRE(unchangedMatrix == Mat3::identity());

    const Mat3 anisotropic = Mat3::diagonal({1.0, 1.0e-8, 1.0e-8});
    Mat3 anisotropicInverse;
    FUNDEM_TEST_REQUIRE(tryInverse(anisotropic, anisotropicInverse));
    FUNDEM_TEST_REQUIRE(nearlyEqual(anisotropic * anisotropicInverse, Mat3::identity()));

    const Mat3 illConditioned = Mat3::diagonal({1.0, 1.0, 1.0e-14});
    FUNDEM_TEST_REQUIRE(!tryInverse(illConditioned, unchangedMatrix));

    Quaternion unchangedQuaternion{2.0, 3.0, 4.0, 5.0};
    FUNDEM_TEST_REQUIRE(!tryInverse(Quaternion::zero(), unchangedQuaternion));
    FUNDEM_TEST_REQUIRE(unchangedQuaternion == Quaternion(2.0, 3.0, 4.0, 5.0));

    Quaternion invalidAxisResult{2.0, 3.0, 4.0, 5.0};
    FUNDEM_TEST_REQUIRE(!Quaternion::tryFromAxisAngle(Vec3::zero(), halfPi, invalidAxisResult));
    FUNDEM_TEST_REQUIRE(invalidAxisResult == Quaternion(2.0, 3.0, 4.0, 5.0));

    Quaternion rotationVectorResult;
    FUNDEM_TEST_REQUIRE(Quaternion::tryFromRotationVector({0.0, 0.0, halfPi}, rotationVectorResult));
    FUNDEM_TEST_REQUIRE(representsSameRotation(rotationVectorResult, hostResult.orientation));

    Quaternion bodyIntegrated;
    FUNDEM_TEST_REQUIRE(tryIntegrateBodyAngularVelocity(Quaternion::identity(), {0.1, 0.2, 0.3}, 0.01, bodyIntegrated));
    FUNDEM_TEST_REQUIRE(nearlyEqual(norm(bodyIntegrated), 1.0));

    Mat3 invalidRotationResult = Mat3::zero();
    FUNDEM_TEST_REQUIRE(!tryRotationMatrix(Quaternion::zero(), invalidRotationResult));
    FUNDEM_TEST_REQUIRE(invalidRotationResult == Mat3::zero());

    const Real infinity = std::numeric_limits<Real>::infinity();
    const Real notANumber = std::numeric_limits<Real>::quiet_NaN();
    FUNDEM_TEST_REQUIRE(nearlyEqual(infinity, infinity));
    FUNDEM_TEST_REQUIRE(!nearlyEqual(notANumber, notANumber));
    FUNDEM_TEST_REQUIRE(!nearlyZero(notANumber));

    const Vec3 first{1.0, 2.0, 3.0};
    const Vec3 second{-2.0, 4.0, 1.0};
    FUNDEM_TEST_REQUIRE(crossProductMatrix(first) * second == cross(first, second));
    FUNDEM_TEST_REQUIRE(outerProduct(first, second) * Vec3::unitY() == first * second.y);

    Evaluation* deviceResult = nullptr;
    FUNDEM_TEST_REQUIRE_CUDA(cudaMalloc(reinterpret_cast<void**>(&deviceResult), sizeof(Evaluation)));

    evaluateKernel<<<1, 1>>>(deviceResult);
    FUNDEM_TEST_REQUIRE_CUDA(cudaGetLastError());

    Evaluation copiedResult;
    FUNDEM_TEST_REQUIRE_CUDA(cudaMemcpy(&copiedResult, deviceResult, sizeof(Evaluation), cudaMemcpyDeviceToHost));
    FUNDEM_TEST_REQUIRE_CUDA(cudaFree(deviceResult));

    verify(copiedResult);
    FUNDEM_TEST_REQUIRE(nearlyEqual(copiedResult.normalized, hostResult.normalized));
    FUNDEM_TEST_REQUIRE(nearlyEqual(copiedResult.inverse, hostResult.inverse));
    FUNDEM_TEST_REQUIRE(representsSameRotation(copiedResult.orientation, hostResult.orientation));
}
