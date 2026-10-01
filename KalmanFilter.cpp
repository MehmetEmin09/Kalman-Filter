#include "KalmanFilter.h"
#include "Kalman.h"
#include <cmath>

namespace OUKB
{

    namespace
    {
        bool gInitialized = false;
        uint32_t gLastTickUs = 0;
        float gMaxAltitude = -1e9f;
        float gVelX = 0.0f;
        float gVelY = 0.0f;
        bool gSutGroundCaptured = false;
        float gBootGroundPressure = 0.0f; // when we stop the SUT test, this ground pressure will be set as P_GROUND
    }

    void kalmanFilterInit(const float *groundPressureSamplesHpa, const float *groundTempSamplesCelsius,
                           const float *groundAccelXG, const float *groundAccelYG, const float *groundAccelZG,
                           int sampleCount)
    {
        calibrateGround(groundPressureSamplesHpa, groundTempSamplesCelsius, sampleCount);

        gBootGroundPressure = P_ground;

        if (sampleCount > 0)
        {
            float ax = groundAccelXG[0] * GRAVITY;
            float ay = groundAccelYG[0] * GRAVITY;
            float az = groundAccelZG[0] * GRAVITY;
            calibrateGravityFromFirstSample(ax, ay, az);
        }

        kalmanInit();
        gInitialized = true;
        gLastTickUs = 0;
        gMaxAltitude = -1e9f;
        gVelX = 0.0f;
        gVelY = 0.0f;
    }

    void kalmanFilterTick(const RawSensorData &raw, FilteredState &out, uint32_t nowUs)
    {
        if (!gInitialized || !raw.imuValid)
        {
            out.valid = false;
            return;
        }

        float dt = (gLastTickUs == 0) ? 0.01f : (float)(nowUs - gLastTickUs) / 1000000.0f;
        gLastTickUs = nowUs;

        decomposition();
        generateSigmaPoints();
        prediction(dt);

        float measAltitude = raw.baroValid ? calculateISAAltitude(raw.pressurePa) : x_state[3];



        float measAx = raw.accelBodyXG * GRAVITY;
        float measAy = raw.accelBodyYG * GRAVITY;
        float measAz = raw.accelBodyZG * GRAVITY;

        update(measAx, measAy, measAz, measAltitude,
               raw.orientationQx, raw.orientationQy, raw.orientationQz, raw.orientationQw);

        out.altitudeM = x_state[3];
        out.accelXMs2 = x_state[0];
        out.accelYMs2 = x_state[1];
        out.accelZMs2 = x_state[2];
        out.quaternion_x = x_state[5];
        out.quaternion_y = x_state[6];
        out.quaternion_z = x_state[7];
        out.quaternion_w = x_state[8];

        gVelX += out.accelXMs2 * dt;
        gVelY += out.accelYMs2 * dt;

        out.vertical_vel = x_state[4];
        out.horizontalVelocityMs = sqrtf(gVelX * gVelX + gVelY * gVelY);

        if (out.altitudeM > gMaxAltitude)
        {
            gMaxAltitude = out.altitudeM;
        }
        out.maxAltitudeM = gMaxAltitude;
        out.timestampUs = nowUs;

        bool finiteOk = std::isfinite(out.altitudeM) && std::isfinite(out.accelXMs2) &&
                         std::isfinite(out.accelYMs2) && std::isfinite(out.accelZMs2) &&
                         std::isfinite(out.quaternion_x) && std::isfinite(out.quaternion_y) &&
                         std::isfinite(out.quaternion_z) && std::isfinite(out.quaternion_w) &&
                         std::isfinite(out.vertical_vel) && std::isfinite(out.horizontalVelocityMs);

        out.valid = finiteOk;

        if (!finiteOk)
{
    kalmanReset();
    gLastTickUs = 0;
}
        
    }

    void kalmanFilterReset()
{
    kalmanReset();          // Kalman.cpp'deki: x_state, P, flying-detector, gravity_calibrated
    gLastTickUs = 0;
    gMaxAltitude = -1e9f;
    gVelX = 0.0f;
    gVelY = 0.0f;
    gSutGroundCaptured = false;

    if (gInitialized) {
        P_ground = gBootGroundPressure;
    }
}

void kalmanFilterCaptureSutGroundPressure(float pressurePa)
{
    if (!gSutGroundCaptured)
    {
        setGroundPressureOverride(pressurePa);
        gSutGroundCaptured = true;
    }
}

} // namespace OUKB