[//]: # (Image References) 
[FC]: ./data/no-log/flow_chart.png
[P20]: ./data/no-log/20.png
[P2000]: ./data/no-log/2000p.png
[RS]: ./data/no-log/results.png
[VT]: ./data/no-log/vtune.png
# Overview
This repository contains all the code needed to complete the final project for the Localization course in Udacity's Self-Driving Car Nanodegree.

## Project Introduction
Your robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.

## Project Implementation
The C++ program for localization was implemented using following major steps:

![Flow Chart][FC]

#### 1. Initializatoin :
A noisy measurement from GPS sensor was received and used to initialize the position of vehicle. This measurement included the x coordinate, y coordinate (both in m) and the theta (orientation) of vehicle in radian. Particle filter algorithm uses particles to represent the location of vehicle. Many (20 -2000 in my results) particles were created and initialized to locations taken from normal distribution with mean equal to the location received from GPS and standard deviation equal to the GPS measurement uncertainty. 

#### 2. Prediction : 
Global map of environment is initialized. This map is represented by a list x and y coordinates of landmarks in the environment. the prediction approach is used "bicycle model" with following equations.

	xf=x0+θ˙v[sin(θ0+θ˙(dt))−sin(θ0)]
	yf=y0+θ˙v[cos(θ0)−cos(θ0+θ˙(dt))]
	θf=θ0+θ˙(dt)

#### 3. Update:  
 Location update is done with the 4 steps listed in below
    
    Step 1: Transform observations from vehicle co-ordinates to map co-ordinates.
	using Homogenous Transformation
	xm=xp+(cosθ×xc)−(sinθ×yc)
	ym=yp+(sinθ×xc)+(cosθ×yc)       
    Step 2:  keep only map landmarks which are in the sensor_range of current
     particle, and save them into predictions vector.     
    Step 3: Associate observations to predicted landmarks by nearest neighbor algorithm.
    Step 4: get the weight of each particle by Multivariate Gaussian distribution.
    Step 5: Normalize the weights of all particles.
#### 4. Resample :  
After prediction step, the vehicle implements Update step. In this step, particles are assigned with weights corresponding to their prediction.

## Project Results
couple tests were done for different number of particles. Some hotspot functions are also identifed via a perf profiler.

here is a result for a test with 20 particles 

![20 particles][P20]

here is a result for a test with 2000 particle

![2000 paricles][P2000]

many tests with different number of particles are conducted. the system time keeps within 48.8sec when there are few than 500 partciles,and the error rate doesn't improve too much with > 500 particles. nearest neighbor algorithms may not work well when there are too many particles, so I used 500 particles in the end.

![result][RS]

Due to sensitivity for running latency, I also used a profiler to identify some hotspot functions.
the main 2 top hotspot functions are:

###### ParticleFilter::updateWeights (used  8% of total CPU time)
###### ParticleFilter::dataAssociation (used  5% of total CPU time)

moreover, the program uses single-thread for computation as below diagram for a profiler.

![profiling_result][VT]

Therefore, may try to accelerate the performance by mulit-threading those two hotspot functions as the next step.

## Running the Code
This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and intall uWebSocketIO for either Linux or Mac systems. For windows you can use either Docker, VMware, or even Windows 10 Bash on Ubuntu to install uWebSocketIO.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./particle_filter

Alternatively some scripts have been included to streamline this process, these can be leveraged by executing the following in the top directory of the project:

1. ./clean.sh
2. ./build.sh
3. ./run.sh

Tips for setting up your environment can be found [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)

Note that the programs that need to be written to accomplish the project are src/particle_filter.cpp, and particle_filter.h

The program main.cpp has already been filled out, but feel free to modify it.

Here is the main protcol that main.cpp uses for uWebSocketIO in communicating with the simulator.

INPUT: values provided by the simulator to the c++ program

// sense noisy position data from the simulator

["sense_x"] 

["sense_y"] 

["sense_theta"] 

// get the previous velocity and yaw rate to predict the particle's transitioned state

["previous_velocity"]

["previous_yawrate"]

// receive noisy observation data from the simulator, in a respective list of x/y values

["sense_observations_x"] 

["sense_observations_y"] 


OUTPUT: values provided by the c++ program to the simulator

// best particle values used for calculating the error evaluation

["best_particle_x"]

["best_particle_y"]

["best_particle_theta"] 

//Optional message data used for debugging particle's sensing and associations

// for respective (x,y) sensed positions ID label 

["best_particle_associations"]

// for respective (x,y) sensed positions

["best_particle_sense_x"] <= list of sensed x positions

["best_particle_sense_y"] <= list of sensed y positions


Your job is to build out the methods in `particle_filter.cpp` until the simulator output says:

```
Success! Your particle filter passed!
```

# Implementing the Particle Filter
The directory structure of this repository is as follows:

```
root
|   build.sh
|   clean.sh
|   CMakeLists.txt
|   README.md
|   run.sh
|
|___data
|   |   
|   |   map_data.txt
|   
|   
|___src
    |   helper_functions.h
    |   main.cpp
    |   map.h
    |   particle_filter.cpp
    |   particle_filter.h
```

The only file you should modify is `particle_filter.cpp` in the `src` directory. The file contains the scaffolding of a `ParticleFilter` class and some associated methods. Read through the code, the comments, and the header file `particle_filter.h` to get a sense for what this code is expected to do.

If you are interested, take a look at `src/main.cpp` as well. This file contains the code that will actually be running your particle filter and calling the associated methods.

## Inputs to the Particle Filter
You can find the inputs to the particle filter in the `data` directory. 

#### The Map*
`map_data.txt` includes the position of landmarks (in meters) on an arbitrary Cartesian coordinate system. Each row has three columns
1. x position
2. y position
3. landmark id

### All other data the simulator provides, such as observations and controls.

> * Map data provided by 3D Mapping Solutions GmbH.

## Success Criteria
If your particle filter passes the current grading code in the simulator (you can make sure you have the current version at any time by doing a `git pull`), then you should pass! 

The things the grading code is looking for are:


1. **Accuracy**: your particle filter should localize vehicle position and yaw to within the values specified in the parameters `max_translation_error` and `max_yaw_error` in `src/main.cpp`.

2. **Performance**: your particle filter should complete execution within the time of 100 seconds.

## How to write a README
A well written README file can enhance your project and portfolio.  Develop your abilities to create professional README files by completing [this free course](https://www.udacity.com/course/writing-readmes--ud777).



