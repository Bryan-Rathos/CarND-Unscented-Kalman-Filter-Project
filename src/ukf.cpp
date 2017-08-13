#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() 
{
  // Set to false initailly. Set to true in the first call to ProcessMeasurement() 
  is_initialized_ = false;

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

  // Inititialize timestamp, in microseconds
  time_us_ = 0.0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Predicted sigma point matrix
  Xsig_pred_ =  MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*********************
  * Initialization 
  *********************/
  if((meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ ) || 
     (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ ) )
  {
    if(!is_initialized_)
    {
      x_ << 1, 1, 1, 1, 1;

      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true)
      {
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
      }

      if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true)
      {
        float rho      = meas_package.raw_measurements_(0);
        float phi      = meas_package.raw_measurements_(1);
        float rho _dot = meas_package.raw_measurements_(3);

        x_(0) = rho * cos(phi);
        x_(1) = rho * sin(phi);
      }

      // Record current timestamp
      time_us_ = meas_package.timestamp_;

      // Done initializing, use sensor measurements from next time stamp
      is_initialized_ = true;        
      return;
    }
    
    /**********************
    * Prediction
    **********************/
    //compute the time elapsed between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
    time_us_ = measurement_pack.timestamp_;

    Prediction(dt);

    /**********************
    * Update
    **********************/
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      UpdateLidar(meas_package);
    }
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      UpdateRadar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
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
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

