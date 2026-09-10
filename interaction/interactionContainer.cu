#include "interactionContainer.h"

#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/system/cuda/execution_policy.h>

namespace fundem
{

void interactionContainer::buildPrefixSum(neighborRangeContainer& ranges, cudaStream_t stream)
{
    const int valueCount = static_cast<int>(ranges.deviceSize());
    if (valueCount <= 0)
    {
        return;
    }
    thrust::inclusive_scan(thrust::cuda::par.on(stream),
                           thrust::device_pointer_cast(ranges.device<neighborRange::countField>()),
                           thrust::device_pointer_cast(ranges.device<neighborRange::countField>() + valueCount),
                           thrust::device_pointer_cast(ranges.device<neighborRange::prefixSumField>()));
    host_device_detail::checkCuda(cudaGetLastError(), "interaction prefix sum");
}

} // namespace fundem
