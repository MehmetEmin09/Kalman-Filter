#include "Kalman.h"
#include <algorithm>
#include <cmath>

namespace OUKB
{
    float GRAVITY = 9.80665f;
    bool gravity_calibrated = false;

    float alpha = 0.5f;
    float beta = 2.0f;
    float kappa = 0.0f;
    float lambda;

    float Wm[SIGMA_COUNT];
    float Wc[SIGMA_COUNT];

    float x_state[STATE_DIM] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f};

   float P[STATE_DIM][STATE_DIM] = {
        {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.05f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.05f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.05f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.05f}};

    float R[MEAS_DIM][MEAS_DIM] = {
        {0.032f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.032f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.032f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.01f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.01f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.01f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.01f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.01f}};

    
    float Q[STATE_DIM][STATE_DIM] = {
        {0.001f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},   //ax
        {0.0f, 0.001f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},   //ay
        {0.0f, 0.0f, 0.001f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},   //az
        {0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},     //h -> 0.001f -> 2.0f
        {0.0f, 0.0f, 0.0f, 0.0f, 1.5f, 0.0f, 0.0f, 0.0f, 0.0f},     //v -> 0.001f -> 1.5f
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0001f, 0.0f, 0.0f, 0.0f},  //qx
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0001f, 0.0f, 0.0f},  //qy
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0001f, 0.0f},  //qz
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0001f}}; //qw

    float A[STATE_DIM][STATE_DIM];
    float L[STATE_DIM][STATE_DIM];

    float sigma_points[SIGMA_COUNT][STATE_DIM];
    float sigma_points_pred[SIGMA_COUNT][STATE_DIM];

    float P_ground = 101320.0f;
    float T_ground = 288.15f;

    void calibrateGround(const float *pressure_samples_hpa, const float *temp_samples_celsius, int n)
    {
        float p_sum = 0.0f, t_sum = 0.0f;
        for (int i = 0; i < n; i++) { p_sum += pressure_samples_hpa[i]; t_sum += temp_samples_celsius[i]; }
        P_ground = (p_sum / (float)n) * 100.0f;
        T_ground = (t_sum / (float)n) + 273.15f;
    }

    void calibrateGravityFromFirstSample(float ax, float ay, float az)
    {
        if (gravity_calibrated) return;
        float measuredG = sqrt(ax * ax + ay * ay + az * az);
        if (measuredG > 9.0f && measuredG < 10.6f) GRAVITY = measuredG;
        gravity_calibrated = true;
    }

    static const float FLIGHT_ACCEL_DEVIATION_THRESHOLD = 3.0f;
    static const int FLIGHT_CONFIRM_COUNT = 5;
    static const int LANDING_CONFIRM_COUNT = 50;
    static bool flying_state = false;
    static int above_threshold_streak = 0;
    static int below_threshold_streak = 0;

    bool isFlying(float sigma[SIGMA_COUNT][STATE_DIM])
    {
        float ax = sigma[0][0], ay = sigma[0][1], az = sigma[0][2];
        float accel_mag = sqrt(ax * ax + ay * ay + az * az);
        float deviation = fabs(accel_mag - GRAVITY);
        if (deviation > FLIGHT_ACCEL_DEVIATION_THRESHOLD) { above_threshold_streak++; below_threshold_streak = 0; }
        else { below_threshold_streak++; above_threshold_streak = 0; }
        if (!flying_state && above_threshold_streak >= FLIGHT_CONFIRM_COUNT) flying_state = true;
        else if (flying_state && below_threshold_streak >= LANDING_CONFIRM_COUNT) flying_state = false;
        return flying_state;
    }

    void normalizeQuaternion(float &qx, float &qy, float &qz, float &qw)
    {
        float mag = sqrt(qx * qx + qy * qy + qz * qz + qw * qw);
        if (mag > 1e-6f) { qx /= mag; qy /= mag; qz /= mag; qw /= mag; }
        else { qx = 0.0f; qy = 0.0f; qz = 0.0f; qw = 1.0f; }
    }

    void quaternionToRotationMatrix(float qx, float qy, float qz, float qw, float Rm[3][3])
    {
        Rm[0][0] = 1 - 2 * (qy * qy + qz * qz);
        Rm[0][1] = 2 * (qx * qy - qz * qw);
        Rm[0][2] = 2 * (qx * qz + qy * qw);
        Rm[1][0] = 2 * (qx * qy + qz * qw);
        Rm[1][1] = 1 - 2 * (qx * qx + qz * qz);
        Rm[1][2] = 2 * (qy * qz - qx * qw);
        Rm[2][0] = 2 * (qx * qz - qy * qw);
        Rm[2][1] = 2 * (qy * qz + qx * qw);
        Rm[2][2] = 1 - 2 * (qx * qx + qy * qy);
    }

    void kalmanInit()
    {
        lambda = (alpha * alpha * (STATE_DIM + kappa)) - STATE_DIM;
        Wm[0] = lambda / (STATE_DIM + lambda);
        Wc[0] = Wm[0] + (1.0f - (alpha * alpha) + beta);
        float weight_others = 1.0f / (2.0f * (STATE_DIM + lambda));
        for (int row = 1; row < SIGMA_COUNT; row++) { Wm[row] = weight_others; Wc[row] = weight_others; }
    }

    void decomposition()
    {
        float scale = STATE_DIM + lambda;
        for (int row = 0; row < STATE_DIM; row++)
            for (int col = 0; col < STATE_DIM; col++) { A[row][col] = scale * P[row][col]; L[row][col] = 0.0f; }

        for (int row = 0; row < STATE_DIM; row++)
        {
            for (int col = 0; col <= row; col++)
            {
                float sum = 0.0f;
                for (int inner = 0; inner < col; inner++) sum += L[row][inner] * L[col][inner];
                if (row == col) { float val = A[row][row] - sum; L[row][col] = (val > 0.0f) ? sqrt(val) : 0.0f; }
                else { L[row][col] = (L[col][col] > 1e-6f) ? (A[row][col] - sum) / L[col][col] : 0.0f; }
            }
        }
    }

    void generateSigmaPoints()
    {
        for (int col = 0; col < STATE_DIM; col++) sigma_points[0][col] = x_state[col];
        for (int row = 0; row < STATE_DIM; row++)
            for (int col = 0; col < STATE_DIM; col++) sigma_points[row + 1][col] = x_state[col] + L[col][row];
        for (int row = 0; row < STATE_DIM; row++)
            for (int col = 0; col < STATE_DIM; col++) sigma_points[row + 1 + STATE_DIM][col] = x_state[col] - L[col][row];

        for (int row = 0; row < SIGMA_COUNT; row++)
        {
            float qx = sigma_points[row][5], qy = sigma_points[row][6], qz = sigma_points[row][7], qw = sigma_points[row][8];
            normalizeQuaternion(qx, qy, qz, qw);
            if (row > 0)
            {
                float dot = qx * sigma_points[0][5] + qy * sigma_points[0][6] + qz * sigma_points[0][7] + qw * sigma_points[0][8];
                if (dot < 0.0f) { qx = -qx; qy = -qy; qz = -qz; qw = -qw; }
            }
            sigma_points[row][5] = qx; sigma_points[row][6] = qy; sigma_points[row][7] = qz; sigma_points[row][8] = qw;
        }
    }

    void prediction(float dt)
    {
        for (int row = 0; row < SIGMA_COUNT; row++)
        {
            float ax_body = sigma_points[row][0], ay_body = sigma_points[row][1], az_body = sigma_points[row][2];
            float altitude = sigma_points[row][3], vz = sigma_points[row][4];
            float qx = sigma_points[row][5], qy = sigma_points[row][6], qz = sigma_points[row][7], qw = sigma_points[row][8];

            float Rm[3][3];
            quaternionToRotationMatrix(qx, qy, qz, qw, Rm);
            float aWorldZ = Rm[2][0] * ax_body + Rm[2][1] * ay_body + Rm[2][2] * az_body;
            float azLinWorld = aWorldZ - GRAVITY;

            sigma_points_pred[row][0] = ax_body;
            sigma_points_pred[row][1] = ay_body;
            sigma_points_pred[row][2] = az_body;
            sigma_points_pred[row][3] = altitude + vz * dt;
            sigma_points_pred[row][4] = vz + azLinWorld * dt;
            sigma_points_pred[row][5] = qx;
            sigma_points_pred[row][6] = qy;
            sigma_points_pred[row][7] = qz;
            sigma_points_pred[row][8] = qw;
        }

        for (int col = 0; col < STATE_DIM; col++)
        {
            x_state[col] = 0.0f;
            for (int row = 0; row < SIGMA_COUNT; row++) x_state[col] += Wm[row] * sigma_points_pred[row][col];
        }
        normalizeQuaternion(x_state[5], x_state[6], x_state[7], x_state[8]);

        for (int row = 0; row < STATE_DIM; row++)
            for (int col = 0; col < STATE_DIM; col++) P[row][col] = Q[row][col];

        for (int inner = 0; inner < SIGMA_COUNT; inner++)
        {
            for (int row = 0; row < STATE_DIM; row++)
            {
                float wc_diff_row = Wc[inner] * (sigma_points_pred[inner][row] - x_state[row]);
                for (int col = 0; col < STATE_DIM; col++)
                    P[row][col] += wc_diff_row * (sigma_points_pred[inner][col] - x_state[col]);
            }
        }
    }

    static bool invert4x4(const float M[4][4], float Minv[4][4])
    {
        float aug[4][8];
        for (int r = 0; r < 4; r++)
        {
            for (int c = 0; c < 4; c++) aug[r][c] = M[r][c];
            for (int c = 0; c < 4; c++) aug[r][4 + c] = (r == c) ? 1.0f : 0.0f;
        }
        for (int col = 0; col < 4; col++)
        {
            int pivotRow = col;
            float maxVal = fabs(aug[col][col]);
            for (int r = col + 1; r < 4; r++) if (fabs(aug[r][col]) > maxVal) { maxVal = fabs(aug[r][col]); pivotRow = r; }
            if (maxVal < 1e-9f) return false;
            if (pivotRow != col) for (int c = 0; c < 8; c++) std::swap(aug[col][c], aug[pivotRow][c]);
            float pivot = aug[col][col];
            for (int c = 0; c < 8; c++) aug[col][c] /= pivot;
            for (int r = 0; r < 4; r++)
            {
                if (r == col) continue;
                float factor = aug[r][col];
                for (int c = 0; c < 8; c++) aug[r][c] -= factor * aug[col][c];
            }
        }
        for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) Minv[r][c] = aug[r][4 + c];
        return true;
    }

    void update(float meas_ax, float meas_ay, float meas_az, float meas_altitude,
                float meas_qx, float meas_qy, float meas_qz, float meas_qw)
    {
        float dot = meas_qx * x_state[5] + meas_qy * x_state[6] + meas_qz * x_state[7] + meas_qw * x_state[8];
        if (dot < 0.0f) { meas_qx = -meas_qx; meas_qy = -meas_qy; meas_qz = -meas_qz; meas_qw = -meas_qw; }

        const float CHI2_GATE_SUSPECT = 3.84f;
        const float CHI2_GATE_REJECT = 10.83f;

        float z_scalar[4] = {meas_ax, meas_ay, meas_az, meas_altitude};
        int scalar_state_map[4] = {0, 1, 2, 3};

        for (int i = 0; i < 4; i++)
        {
            int s = scalar_state_map[i];
            float y = z_scalar[i] - x_state[s];
            float S_i = P[s][s] + R[i][i];
            if (S_i < 1e-9f) continue;

            float gate_scale = 1.0f;

            if (i == 3) {
    float normalized_error = (y * y) / S_i;
    if (isFlying(sigma_points_pred)) {
        if (normalized_error > CHI2_GATE_REJECT) {
            gate_scale = 0.0f;
        } else if (normalized_error > CHI2_GATE_SUSPECT) {
            gate_scale = CHI2_GATE_SUSPECT / normalized_error;
        }
    }
}

y *= gate_scale;

float P_s_row[STATE_DIM];
for (int k = 0; k < STATE_DIM; k++) P_s_row[k] = P[s][k];

float K[STATE_DIM];
for (int j = 0; j < STATE_DIM; j++) {
    K[j] = gate_scale * (P_s_row[j] / S_i);
    x_state[j] += K[j] * y;
}

for (int j = 0; j < STATE_DIM; j++) {
    float k_j = K[j];
    for (int k = 0; k < STATE_DIM; k++) P[j][k] -= k_j * P_s_row[k];
}
        }

        int q_idx[4] = {5, 6, 7, 8};
        float z_q[4] = {meas_qx, meas_qy, meas_qz, meas_qw};
        float y_q[4];
        for (int i = 0; i < 4; i++) y_q[i] = z_q[i] - x_state[q_idx[i]];

        float S_q[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) S_q[i][j] = P[q_idx[i]][q_idx[j]] + R[4 + i][4 + j];

        float S_q_inv[4][4];
        if (invert4x4(S_q, S_q_inv))
        {
            float P_cols_q[STATE_DIM][4];
            for (int j = 0; j < STATE_DIM; j++)
                for (int i = 0; i < 4; i++) P_cols_q[j][i] = P[j][q_idx[i]];

            float K[STATE_DIM][4];
            for (int j = 0; j < STATE_DIM; j++)
                for (int i = 0; i < 4; i++)
                {
                    float sum = 0.0f;
                    for (int k = 0; k < 4; k++) sum += P_cols_q[j][k] * S_q_inv[k][i];
                    K[j][i] = sum;
                }

            for (int j = 0; j < STATE_DIM; j++)
            {
                float sum = 0.0f;
                for (int i = 0; i < 4; i++) sum += K[j][i] * y_q[i];
                x_state[j] += sum;
            }

            for (int j = 0; j < STATE_DIM; j++)
                for (int k = 0; k < STATE_DIM; k++)
                {
                    float sum = 0.0f;
                    for (int i = 0; i < 4; i++) sum += K[j][i] * P_cols_q[k][i];
                    P[j][k] -= sum;
                }
        }

        normalizeQuaternion(x_state[5], x_state[6], x_state[7], x_state[8]);

        for (int row = 0; row < STATE_DIM; row++)
        {
            for (int col = row + 1; col < STATE_DIM; col++)
            {
                float avg = (P[row][col] + P[col][row]) * 0.5f;
                P[row][col] = avg; P[col][row] = avg;
            }
            if (P[row][row] < 0.0f) P[row][row] = 1e-4f;
        }
    }

    
    
    void kalmanReset()
{
    x_state[0] = 0.0f; x_state[1] = 0.0f; x_state[2] = 0.0f;
    x_state[3] = 0.0f; x_state[4] = 0.0f;
    x_state[5] = 0.0f; x_state[6] = 0.0f; x_state[7] = 0.0f; x_state[8] = 1.0f;

    for (int row = 0; row < STATE_DIM; row++)
        for (int col = 0; col < STATE_DIM; col++)
            P[row][col] = 0.0f;

    P[0][0] = 1.0f; P[1][1] = 1.0f; P[2][2] = 1.0f;
    P[3][3] = 1.0f; P[4][4] = 1.0f;
    P[5][5] = 0.05f; P[6][6] = 0.05f; P[7][7] = 0.05f; P[8][8] = 0.05f;

    gravity_calibrated = false;      // yoksa yeni senaryonun ilk örneği kalibre edilmez, eski GRAVITY değeri kalır
    flying_state = false;
    above_threshold_streak = 0;
    below_threshold_streak = 0;      // yoksa isFlying() önceki uçuşun streak'iyle başlar → apogee'deki chi2 gating yanlış davranır
}
    
    
    
    float calculateISAAltitude(float raw_pressure)
{
    if (!std::isfinite(raw_pressure) || raw_pressure <= 0.0f || P_ground <= 0.0f) return x_state[3];
    float pressure_ratio = raw_pressure / P_ground;
    return (T_ground / LAPSE_RATE) * (1.0f - pow(pressure_ratio, ISA_EXPONENT));
}

void setGroundPressureOverride(float pressurePa)
{
    P_ground = pressurePa;
}

void setGroundTemperatureOverride(float tempKelvin)
{
    T_ground = tempKelvin;
}
}