#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools tools;

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
  std_a_ = 0.9;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;

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

  is_initialized_ = false;

  // state dimensions
  n_x_ = 5;
  n_aug_ = 7;

  // measurement dimensions
  // radar: rho, phi, r_dot / laser: px, py
  n_z_radar_ = 3;
  n_z_laser_ = 2;

  // Number of sigma points
  n_sigma_points_ = 2 * n_aug_ + 1;

  lambda_ = 3 - n_aug_;

  // Predicted sigma points matrix initialization
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points_);
  Xsig_pred_.fill(0.0);

  // set weights for compensating lambda
  weights_ = VectorXd(n_sigma_points_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sigma_points_; i++) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  // set laser measurement matrix
  H_ = MatrixXd(n_z_laser_, n_x_);
  H_.fill(0.0);
  H_(0, 0) = 1;
  H_(1, 1) = 1;

  // measurement noise covariance matrix
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ *std_laspy_;

  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ *std_radphi_, 0,
              0, 0, std_radrd_ *std_radrd_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // Initialization
  if (!is_initialized_) {
    // initialize state vector x
    x_.fill(0.0);

    // initialize state covariance matrix P
    P_.fill(0.0);
    for(int i = 0; i < 5; i++) {
      P_(i, i) = 1;
    }
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
  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug = CreateAugmentedSigmaPointsMatrix();
  PredictAugmentedSigmaPointsMatrix(delta_t, Xsig_aug);
  PredictMeanAndCovariance();
}

MatrixXd UKF::CreateAugmentedSigmaPointsMatrix() {
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
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_points_);
  Xsig_aug.fill(0.0);

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  return Xsig_aug;
}

void UKF::PredictAugmentedSigmaPointsMatrix(double delta_t, MatrixXd Xsig_aug) {
  double p_x, p_y, v, yaw, yawd, nu_a, nu_yawdd;
  double px_p, py_p;

  for (int i = 0; i < n_sigma_points_; i++) {
    // extract values for better readability
    p_x = Xsig_aug(0, i);
    p_y = Xsig_aug(1, i);
    v = Xsig_aug(2, i);
    yaw = Xsig_aug(3, i);
    yawd = Xsig_aug(4, i);
    nu_a = Xsig_aug(5, i);
    nu_yawdd = Xsig_aug(6, i);

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    // add noise and write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);  // px_p
    Xsig_pred_(1, i) = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);  // py_p
    Xsig_pred_(2, i) = v + nu_a * delta_t; // v_p
    Xsig_pred_(3, i) = yaw + yawd * delta_t + 0.5 * nu_yawdd * delta_t * delta_t;  // yaw_p
    Xsig_pred_(4, i) = yawd + nu_yawdd * delta_t;  // yawd_p
  }
}

void UKF::PredictMeanAndCovariance() {
  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) { // iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) { // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    tools.AngleNormalization(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_laser_;
  VectorXd y = z - H_ * x_;
  UpdateStateFromLidar(S, y);
  NIS_laser_ = y.transpose() * S.inverse() * y;
}

/**
 * new estimate
 */
void UKF::UpdateStateFromLidar(MatrixXd S, VectorXd y) {
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ += K * y;
  P_ = (MatrixXd::Identity(n_x_, n_x_) - K * H_) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd Zsig_pred = CreatePredictedSigmaPointsMatrix();
  VectorXd z_pred = CreateMeanPredictedMeasurements(Zsig_pred);
  MatrixXd S = PredictRadarMeasurement(Zsig_pred, z_pred);
  UpdateStateFromRadar(z, Zsig_pred, z_pred, S);
  VectorXd y = z - z_pred;
  NIS_radar_ = y.transpose() * S.inverse() * y;
}

MatrixXd UKF::CreatePredictedSigmaPointsMatrix() {
  // Predicted sigma points initialization
  MatrixXd Zsig_pred(n_z_radar_, n_sigma_points_);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sigma_points_; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_pred(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_pred(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_pred(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  return Zsig_pred;
}

VectorXd UKF::CreateMeanPredictedMeasurements(MatrixXd Zsig_pred) {
  VectorXd z_pred(n_z_radar_);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {
    z_pred += weights_(i) * Zsig_pred.col(i);
  }
  return z_pred;
}

MatrixXd UKF::PredictRadarMeasurement(MatrixXd Zsig_pred, VectorXd z_pred) {
  MatrixXd S(n_z_radar_, n_z_radar_);  // measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;  // residual
    tools.AngleNormalization(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  // add measurement noise covariance matrix and return
  return S + R_radar_;
}

/**
 * new estimate
 */
void UKF::UpdateStateFromRadar(VectorXd raw_measurements, MatrixXd Zsig_pred, VectorXd z_pred, MatrixXd S) {
  MatrixXd Tc(n_x_, n_z_radar_);  // cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) { // 2n+1 sigma points
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;  // residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_;  // state difference
    tools.AngleNormalization(z_diff(1));
    tools.AngleNormalization(x_diff(3));
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K = Tc * S.inverse();  // Kalman gain
  VectorXd z_diff = raw_measurements - z_pred;  // residual
  tools.AngleNormalization(z_diff(1));
  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += K * S * K.transpose();
}
