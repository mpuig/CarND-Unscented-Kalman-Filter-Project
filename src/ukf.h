#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

private:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  long previous_timestamp_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Number of sigma points
  int n_sigma_points_;

  ///* Radar measurement dimension
  int n_z_radar_;

  ///* Lidar measurement dimension
  int n_z_laser_;

  ///* Predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* Measurement matrix for laser
  MatrixXd H_;

  ///* Measurement noise laser
  MatrixXd R_laser_;

  ///* Measurement noise radar
  MatrixXd R_radar_;

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
   * Create the augmented sigma points matrix.
   */
  MatrixXd CreateAugmentedSigmaPointsMatrix();

  /**
   * Create the predicted sigma points matrix.
   */
  MatrixXd CreatePredictedSigmaPointsMatrix();

  /**
   * Create the mean predicted measurement.
   */
  VectorXd CreateMeanPredictedMeasurements(MatrixXd Zsig_pred);

  /**
   * Predicts the sigma points.
   * @param delta_t Time between k and k+1 in s
   */
  void PredictAugmentedSigmaPointsMatrix(double delta_t, MatrixXd Xsig_out);

  /**
   * Predicts the state, and the state covariance matrix.
   */
  void PredictMeanAndCovariance();

  /**
   * Predicts the Radar measurements
   */
  MatrixXd PredictRadarMeasurement(MatrixXd Zsig_pred, VectorXd z_pred);

  /**
   * Updates the State from Lidar measurements
   */
  void UpdateStateFromLidar(MatrixXd S, VectorXd y);

  /**
   * Updates the State from Radar measurements
   */
  void UpdateStateFromRadar(VectorXd raw_measurements, MatrixXd Zsig_pred, VectorXd z_pred, MatrixXd S);
};

#endif /* UKF_H */
