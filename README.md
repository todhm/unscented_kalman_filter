# Unscented Kalman Filter Project
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

[image1]: ./result_pic/screen_ds1.png
[image2]: ./result_pic/screen_ds2.png

This is part of udacity self-driving car nanodegree project part2. In this project unscented kalman filter model was built based on c++  to estimate the state of a moving object in the car simulator provided by Udacity.  By lidar and radar measurements. This code consists of following processes.
* Read the data from project simulator which provide radar and laser measurement of moving object near by car in the simulator.
* Estimate the position and velocity of the object by kalman filter model.
* Represent the estimation and RMSE value on the simulator.

The source code of this project consists with followings.
* main.cpp: Reads the data and send the processed data to simulator.
* ukf.h/ukf.cpp: code to predict and update state of moving object based on laser and radar sensor and the fuse both sensor. 
* tools.cpp: tools to calculate RMSE. 

The visualization of this project was made with the Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases).

---

## Basic Build Instructions

I used xcode as main IDE. You can execute project with ide_profiles/xcode/UnscentedKF.xcodeproj file. 

Also you can execute project with following steps.
0. Install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems and [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) for windows.
1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
4. Run it: `./UnscentedKF`  
5. Run the simulator and watch the results.
---

## Results
* The result of model was measured by the RMSE value and NIS value. 
#### RMSE

|    |Fused sensor | only laser | only radar |
|:--:|:----------:|:----------:|:----------:|
| px |   0.0838   |   0.1473   |   0.2302   |
| py |   0.1067   |   0.1153   |   0.3464   |
| vx |   0.3632   |   0.6383   |   0.5835   |
| vy |   0.2574   |   0.5346   |   0.8040   |

#### ratio of NIS that exeed 95% distribution.
|    |radar       | laser     |
|:--:|:----------:|:----------:|
|percent |  0.044   |   0  | 
* Here are the visualization of result in the simulator.

#### Dataset1 result in the simulator.
[![alt text][image1]](https://youtu.be/NueBEqH7qqE)


#### Dataset2 result in the simulator.
[![alt text][image2]](https://youtu.be/YLHAFDCLDt4)

---
## Discussion  
#### Here are the crucial tuning point that was important to improve the reduce RMSE of estimations.
* It was crucial to normalize angle between -pi and pi value.
* It was hard to adjust noise parameter to get adequate RMSE and NIS value. 