# Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../data/output.txt`

## Running modes:

It's possible to run the project in three different modes:

  1. LIDAR and RADAR combined
  2. LIDAR alone
  3. RADAR alone

To change the mode, `use_laser_` and `use_radar_` variables can be defined as `true` or `false`.

## Tuning Process Noise
The values for the process noise `std_a_` and `std_yawdd_` were both initially set to 30. Once all the methods where implemented, a few experiments were executed to adjust its values in order to get the Kalman filter working
The final values are:

```
  std_a = 0.9
  std_yawdd = 0.4
```

## RMSE

Once good values were found, and using the file "obj_pose-laser-radar-synthetic-input.txt" as input, the project was executed with `use_laser_ = false' to generate the file 'data/output-no-laser.txt' and with `use_radar_ = false' to generate the file 'data/output-no-radar.txt'.

The RMSE obtained running with laser and radar was:

```
 0.0664948
 0.0827918
 0.333513
 0.217135
```

The RMSE obtained running only with laser was:

```
 0.0958565
 0.0921837
 0.605361
 0.234694
```

The RMSE obtained running only with radar was:

```
 0.148082
 0.195942
 0.277767
 0.311838
```

## Results visualization

A Jupyter notebook has been used to visualize the obtained results: [UKF Visualization](UKF%20Visualization.ipynb)
