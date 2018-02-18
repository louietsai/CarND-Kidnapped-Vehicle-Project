/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	/* Init every particle. 
	   assign x, y, theta with a value from normal distribution
	   assign 1 as initial weight 
	*/
	num_particles = 20;
	default_random_engine gen;

	normal_distribution<double> distributionX(x, std[0]);
	normal_distribution<double> distributionY(y, std[1]);
	normal_distribution<double> distributionTheta(theta, std[2]);

	int i;
	for (i = 0; i < num_particles; i++) {
	  Particle particle_;
	  particle_.id = i;
	  particle_.x = distributionX(gen);
	  particle_.y = distributionY(gen);
	  particle_.theta = distributionTheta(gen);
	  particle_.weight = 1.0;

	  particles.push_back(particle_);
	  weights.push_back(particle_.weight);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	/*
		predict x,y, theta for each particles by following equations
		xf=x0+θ˙v[sin(θ0+θ˙(dt))−sin(θ0)]
		yf=y0+θ˙v[cos(θ0)−cos(θ0+θ˙(dt))]
		θf=θ0+θ˙(dt)
	*/
	default_random_engine gen;

	int i;
	for (i = 0; i < num_particles; i++) {
	  double pX = particles[i].x;
	  double pY = particles[i].y;
	  double pTheta = particles[i].theta;

	  double predX;
	  double predY;
	  double predTheta;
	  // check if yaw is zero or not
	  if (fabs(yaw_rate) < 0.0001) {
	    predX = pX + velocity * cos(pTheta) * delta_t;
	    predY = pY + velocity * sin(pTheta) * delta_t;
	    predTheta = pTheta;
	  } else {
	    predX = pX + (velocity/yaw_rate) * (sin(pTheta + (yaw_rate * delta_t)) - sin(pTheta));
	    predY = pY + (velocity/yaw_rate) * (cos(pTheta) - cos(pTheta + (yaw_rate * delta_t)));
	    predTheta = pTheta + (yaw_rate * delta_t);
	  }
	  // normalize the predicted x,y,theta , and save to particles
	  normal_distribution<double> distributionX(predX, std_pos[0]);
	  normal_distribution<double> distributionY(predY, std_pos[1]);
	  normal_distribution<double> distributionTheta(predTheta, std_pos[2]);

	  particles[i].x = distributionX(gen);
	  particles[i].y = distributionY(gen);
	  particles[i].theta = distributionTheta(gen);  
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, double sensor_range) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	/*Associate observations in map co-ordinates to predicted landmarks using nearest neighbor algorithm.*/
	int i, j;
	// loop over the observation landmark from current position first
	for (i = 0; i < observations.size(); i++) {
		double lowestDist = sensor_range * sqrt(2);
		int closestLandmarkId = -1;
		double obsX = observations[i].x;
		double obsY = observations[i].y;
		// loop over predicted landmark for current particle
		for (j = 0; j < predicted.size(); j++) {
		  double predX = predicted[j].x;
		  double predY = predicted[j].y;
		  int predId = predicted[j].id;
		  double currentDist = dist(obsX, obsY, predX, predY);

		  if (currentDist < lowestDist) {
		    lowestDist = currentDist;
		    closestLandmarkId = predId;
		  }
		}
		//assocaite the closet predicted landmark with related observed landmark
		observations[i].id = closestLandmarkId;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	//   TODO complete

  int i, j;
  /* used for normalizing weights of all particles within the range
    of [0, 1]*/
  double wNormalizer = 0.0;

  for (i = 0; i < num_particles; i++) {
    double pX = particles[i].x;
    double pY = particles[i].y;
    double pTheta = particles[i].theta;


    /*Step 1: Transform observations from vehicle co-ordinates to map co-ordinates.
	using Homogenous Transformation
	xm=xp+(cosθ×xc)−(sinθ×yc)
	ym=yp+(sinθ×xc)+(cosθ×yc)
    */

    vector<LandmarkObs> transformedObserv;

    //Transform each observations.
    for (j = 0; j < observations.size(); j++) {
      LandmarkObs tObs;
      tObs.id = j;
      tObs.x = pX + (cos(pTheta) * observations[j].x) - (sin(pTheta) * observations[j].y);
      tObs.y = pY + (sin(pTheta) * observations[j].x) + (cos(pTheta) * observations[j].y);
      transformedObserv.push_back(tObs);
    }

    /*Step 2:  keep only map landmarks which are in the sensor_range of current
     particle, and save them into predictions vector.*/
    vector<LandmarkObs> predictedLandmarks;
    for (j = 0; j < map_landmarks.landmark_list.size(); j++) {
      Map::single_landmark_s landmark_ = map_landmarks.landmark_list[j];
      if ((fabs((pX - landmark_.x_f)) <= sensor_range) && (fabs((pY - landmark_.y_f)) <= sensor_range)) {
        predictedLandmarks.push_back(LandmarkObs {landmark_.id_i, landmark_.x_f, landmark_.y_f});
      }
    }

    /*Step 3: Associate observations to predicted landmarks by nearest neighbor algorithm.*/
    dataAssociation(predictedLandmarks, transformedObserv, sensor_range);

    /*Step 4: get the weight of each particle by Multivariate Gaussian distribution.*/
    //Reset the weight of particle to 1.0
    particles[i].weight = 1.0;

    double sigmaX = std_landmark[0];
    double sigmaY = std_landmark[1];
    double sigmaX2 = pow(sigmaX, 2);
    double sigmaY2 = pow(sigmaY, 2);
    double normalizer = (1.0/(2.0 * M_PI * sigmaX * sigmaY));
    int k, l;

    /*get the weight of transformed particle based on the multivariate Gaussian probability function*/
    for (k = 0; k < transformedObserv.size(); k++) {
      double transObsX = transformedObserv[k].x;
      double transObsY = transformedObserv[k].y;
      double transObsId = transformedObserv[k].id;
      double multi_prob = 1.0;

      for (l = 0; l < predictedLandmarks.size(); l++) {
        double predLandmarkX = predictedLandmarks[l].x;
        double predLandmarkY = predictedLandmarks[l].y;
        double predLandmarkId = predictedLandmarks[l].id;

        if (transObsId == predLandmarkId) {
          multi_prob = normalizer * exp(-1.0 * ((pow((transObsX - predLandmarkX), 2)/(2.0 * sigmaX2)) + (pow((transObsY - predLandmarkY), 2)/(2.0 * sigmaY2))));
          particles[i].weight *= multi_prob;
        }
      }
    }
    wNormalizer += particles[i].weight;
  }

  /*Step 5: Normalize the weights of all particles.*/
  for (int i = 0; i < particles.size(); i++) {
    particles[i].weight /= wNormalizer;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	/* 0215 TODO: use a random method, may adopt other method later*/
	vector<Particle> resampledParticles;

	// Create a generator to be used for generating random particle index and beta value
	default_random_engine gen;

	//Generate random particle index
	uniform_int_distribution<int> particle_index(0, num_particles - 1);

	int Index = particle_index(gen);

	double randomBigWeight = 0.0;

	double maxWeight2 = 2.0 * *max_element(weights.begin(), weights.end());

	for (int i = 0; i < particles.size(); i++) {
	  uniform_real_distribution<double> random_weight(0.0, maxWeight2);
	  randomBigWeight += random_weight(gen);
	  // save the particle if its weight is more than random big weight
	  while (randomBigWeight > weights[Index]) {
	    randomBigWeight -= weights[Index];
	    Index = (Index + 1) % num_particles;
	  }
	  resampledParticles.push_back(particles[Index]);
	}
	particles = resampledParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

void ParticleFilter::dumpParticles() {

	for (int i = 0; i < particles.size(); i++) {
                double particle_x = particles[i].x;
                double particle_y = particles[i].y;
                double particle_theta = particles[i].theta;
                double particle_weight = particles[i].weight;
                int particle_id = particles[i].id;
                cout << " p"<< i << " ,id:" << particle_id << " ,weight:"<< particle_weight << " ,x:"<< particle_x << " ,y:"<< particle_y<<endl;

	}
}
