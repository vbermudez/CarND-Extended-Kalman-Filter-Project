#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  // Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
			        0, 1, 0, 0;
  ekf_.x_ = VectorXd(4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  ekf_.Q_ = MatrixXd(4, 4);

  noise_ax_ = 9;
  noise_ay_ = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Initialize(const MeasurementPackage &measurement_pack) {
  if (is_initialized_) {
    return;
  }

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    float r = measurement_pack.raw_measurements_[0]; // Range
    float b = measurement_pack.raw_measurements_[1]; // Bearing
    float rd = measurement_pack.raw_measurements_[2]; // Radial velocity
    float x = r * cos(b);
    float y = r * sin(b);
    float vx = rd * cos(b);
    float vy = rd * sin(b);

    ekf_.x_ << x, y, vx, vy;
  } else {
    ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
  }

  previous_timestamp_ = measurement_pack.timestamp_;
  is_initialized_ = true;
}

void FusionEKF::Predict(const MeasurementPackage &measurement_pack) {
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // ms

  previous_timestamp_ = measurement_pack.timestamp_;

  MatrixXd Q_noise = MatrixXd(2, 2);
  MatrixXd G = MatrixXd(4, 2);
  double dt_2 = dt * dt;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  Q_noise << noise_ax_, 0,
        0, noise_ay_;
  G << dt_2 / 2, 0,
       0, dt_2 / 2,
       dt, 0,
       0, dt;
  ekf_.Q_ = G * Q_noise * G.transpose();
  ekf_.Predict();
}

void FusionEKF::Update(const MeasurementPackage &measurement_pack) {
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  Initialize(measurement_pack);

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   Predict(measurement_pack);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  Update(measurement_pack);

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
