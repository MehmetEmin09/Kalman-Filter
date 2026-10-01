#pragma once
#include "core/Types.h"
#include <cstdint>

namespace OUKB
{

    void kalmanFilterInit(const float *groundPressureSamplesHpa, const float *groundTempSamplesCelsius,
                           const float *groundAccelXG, const float *groundAccelYG, const float *groundAccelZG,
                           int sampleCount);

    void kalmanFilterTick(const RawSensorData &raw, FilteredState &out, uint32_t nowUs);

    void kalmanFilterReset();

    void kalmanFilterCaptureSutGroundPressure(float pressurePa);

} // namespace OUKB