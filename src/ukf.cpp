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

  // radar measurement dimension (r, phi, r_dot)
  n_z_radar_ = 3;

  n_z_laser_ = 2;

  // Number of sigma points
  n_sigma_points_ = 2 * n_aug_ + 1;

  lambda_ = 3 - n_aug_;

  // Predicted sigma points matrix initialization
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points_);
  Xsig_pred_.fill(0.0);

  // Augmented sigma points initialization
  Xsig_aug_ = MatrixXd(n_aug_, n_sigma_points_);
  Xsig_aug_.fill(0.0);

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

  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ *std_laspy_;
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_ * std_radr_, 0, 0, 0, std_radphi_ *std_radphi_, 0, 0, 0,
      std_radrd_ *std_radrd_;

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

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
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
  VectorXd z = meas_package.raw_measurements_;

  MatrixXd S = H_ * P_ * H_.transpose() + R_laser_;
  VectorXd y = z - H_ * x_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // new estimate
  x_ += K * y;
  P_ = (MatrixXd::Identity(n_x_, n_x_) - K * H_) * P_;

  NIS_laser_ = y.transpose() * S.inverse() * y;

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

  // Predicted sigma points initialization
  MatrixXd Zsig_pred(n_z_radar_, n_sigma_points_);
  Zsig_pred.fill(0.0);
  // mean predicted measurement
  VectorXd z_pred(n_z_radar_);
  z_pred.fill(0.0);
  // measurement covariance matrix S
  MatrixXd S(n_z_radar_, n_z_radar_);
  S.fill(0.0);
  PredictRadarMeasurement(&Zsig_pred, &z_pred, &S);
  VectorXd z = meas_package.raw_measurements_;
  UpdateStateFromRadar(z, Zsig_pred, z_pred, S);
  VectorXd y = z - z_pred;
  NIS_radar_ = y.transpose() * S.inverse() * y;
}

void UKF::PredictRadarMeasurement(MatrixXd *Zsig_pred, VectorXd *z_pred, MatrixXd *S) {

  // transform sigma points into measurement space
  double px=0, py=0, v=0, yaw=0;
  for (int i = 0; i < n_sigma_points_; i++) {
    px = Xsig_pred_(0, i);
    py = Xsig_pred_(1, i);
    v = Xsig_pred_(2, i);
    yaw = Xsig_pred_(3, i);
    double rho = sqrt(px*px + py*py);
    // avoid zero division
    if (fabs(rho) < 0.0001) {
      rho = 0.0001;
    }
    // measurement model
    (*Zsig_pred)(0, i) = rho;
    (*Zsig_pred)(1, i) = atan2(py, px);  // phi
    (*Zsig_pred)(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / rho;  // r_dot
  }

  // mean predicted measurement
  for (int i = 0; i < n_sigma_points_; i++) {
    *z_pred += weights_(i) * Zsig_pred->col(i);
  }
  // calculate measurement covariance matrix S
  MatrixXd z_diff;
  z_diff.fill(0.0);

  // residual
  z_diff = Zsig_pred->colwise() - *z_pred;
  for (int i = 0; i < n_sigma_points_; i++) {
    // angle normalization
    while (z_diff.col(i)(1) > M_PI) z_diff.col(i)(1) -= 2. * M_PI;
    while (z_diff.col(i)(1) < -M_PI) z_diff.col(i)(1) += 2. * M_PI;

    *S += weights_(i) * z_diff.col(i) * z_diff.col(i).transpose();
  }

  *S += R_radar_;

}

void UKF::UpdateStateFromRadar(VectorXd raw_measurements, MatrixXd Zsig_pred,
                        VectorXd z_pred, MatrixXd S) {
  // calculate cross correlation matrix
  MatrixXd Tc(n_x_, n_z_radar_);
  Tc.fill(0.0);

  for (int i = 0; i < n_sigma_points_; i++) { // 2n+1 sigma points

    // residual
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z = raw_measurements;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += K * S * K.transpose();
}
