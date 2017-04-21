#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  // state dimensions
  n_x_ = 5;
  n_aug_ = 7;

  // Number of sigma points
  n_sigma_points_ = 2 * n_aug_ + 1;

  lambda_ = 3 - n_aug_;

  // Augmented sigma points initialization
  Xsig_aug_ = MatrixXd(n_aug_, n_sigma_points_);
  Xsig_aug_.fill(0.0);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // Initialization
  if (!is_initialized_) {

    P_.fill(0.0);
    for(int i=0; i<5; i++) P_(i,i) = 1;

    x_.fill(0.0);

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = fabs(meas_package.raw_measurements_[2]);
      x_ << rho *cos(phi), rho * sin(phi), rho_dot, 0, 0;
    }

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Prediction
  float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; // expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(delta_t);
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  CreateAugmentedSigmaPointsMatrix();
  PredictAugmentedSigmaPointsMatrix(delta_t);
  PredictMeanAndCovariance();
}

void UKF::CreateAugmentedSigmaPointsMatrix() {
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug_.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

}

void UKF::PredictAugmentedSigmaPointsMatrix(double delta_t) {

  double p_x=0, p_y=0, v=0, yaw=0, yawd=0, nu_a=0, nu_yawdd=0;

  for (int i = 0; i < n_sigma_points_; i++) {
    //extract values for better readability
    p_x = Xsig_aug_(0, i);
    p_y = Xsig_aug_(1, i);
    v = Xsig_aug_(2, i);
    yaw = Xsig_aug_(3, i);
    yawd = Xsig_aug_(4, i);
    nu_a = Xsig_aug_(5, i);
    nu_yawdd = Xsig_aug_(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    //add noise and write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);  // px_p
    Xsig_pred_(1, i) = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);  // py_p
    Xsig_pred_(2, i) = v + nu_a * delta_t; // v_p
    Xsig_pred_(3, i) = yaw + yawd*delta_t + 0.5 * nu_yawdd * delta_t * delta_t;  // yaw_p
    Xsig_pred_(4, i) = yawd + nu_yawdd*delta_t;  // yawd_p

  }
}

void UKF::PredictMeanAndCovariance() {

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
