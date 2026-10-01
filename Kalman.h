#pragma once
#include <cstdint>
#include <cmath>

namespace OUKB
{
    constexpr uint8_t STATE_DIM = 9;

    constexpr uint8_t MEAS_DIM = 8;
    constexpr uint8_t SIGMA_COUNT = 2 * STATE_DIM + 1;

    extern float GRAVITY;
    extern bool gravity_calibrated;

    constexpr float LAPSE_RATE = 0.0065f;
    constexpr float ISA_EXPONENT = 0.190263f;

    extern float alpha;
    extern float beta;
    extern float kappa;
    extern float lambda;

    extern float Wm[SIGMA_COUNT];
    extern float Wc[SIGMA_COUNT];

    extern float x_state[STATE_DIM];
    extern float P[STATE_DIM][STATE_DIM];
    extern float R[MEAS_DIM][MEAS_DIM];
    extern float Q[STATE_DIM][STATE_DIM];

    extern float A[STATE_DIM][STATE_DIM];
    extern float L[STATE_DIM][STATE_DIM];
    extern float sigma_points[SIGMA_COUNT][STATE_DIM];
    extern float sigma_points_pred[SIGMA_COUNT][STATE_DIM];

    extern float P_ground;
    extern float T_ground;

    void calibrateGround(const float *pressure_samples_hpa, const float *temp_samples_celsius, int n);
    void calibrateGravityFromFirstSample(float ax, float ay, float az);
    void normalizeQuaternion(float &qx, float &qy, float &qz, float &qw);
    bool isFlying(float sigma[SIGMA_COUNT][STATE_DIM]);
    void quaternionToRotationMatrix(float qx, float qy, float qz, float qw, float Rm[3][3]);
    void kalmanInit();
    void decomposition();
    void generateSigmaPoints();
    void prediction(float dt);

    void update(float meas_ax, float meas_ay, float meas_az, float meas_altitude,
                float meas_qx, float meas_qy, float meas_qz, float meas_qw);

    void kalmanReset();

    float calculateISAAltitude(float raw_pressure);

    void setGroundPressureOverride(float pressurePa);
}