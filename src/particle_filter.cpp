/*
 * particle_filter.cpp
 *
 *  Created on: June 21, 2017
 *  Author: Junsheng Fu
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>

#include "particle_filter.h"

using namespace std;

const int NUMBER_OF_PARTICLES = 100;
const double INITIAL_WEIGHT = 1.0;
const double THRESH = 0.0001;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  if (is_initialized) {
    return;
  }
  
  // Generate normal distributions
  num_particles = NUMBER_OF_PARTICLES;
  particles.resize(NUMBER_OF_PARTICLES);
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  default_random_engine gen;

  for(auto& partic: particles){
    partic.x = dist_x(gen);
    partic.y = dist_y(gen);
    partic.theta = dist_theta(gen);
    partic.weight = INITIAL_WEIGHT;
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/
	
  const bool STRAIGHT = fabs(yaw_rate) < THRESH;
  const double vel_delta = velocity * delta_t;
  const double yaw_delta = yaw_rate*delta_t;
  const double vel_yaw = velocity / yaw_rate;
  default_random_engine gen;

  // generate random Gaussian noise
  normal_distribution<double> noise_x(0, std_pos[0]);
  normal_distribution<double> noise_y(0, std_pos[1]);
  normal_distribution<double> noise_theta(0, std_pos[2]);

  for(auto& partic: particles){

    // add measurements to each particle
    if(STRAIGHT){  
      partic.x += vel_delta * cos(partic.theta);
      partic.y += vel_delta * sin(partic.theta);

    } else{
      partic.x += vel_yaw * ( sin( partic.theta + yaw_delta ) - sin(partic.theta) );
      partic.y += vel_yaw * ( cos( partic.theta ) - cos( partic.theta + yaw_delta ) );
      partic.theta += yaw_rate*delta_t;
    }

    // Add noise
    partic.x += noise_x(gen);
    partic.y += noise_y(gen);
    partic.theta += noise_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	// observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	// implement this method and use it as a helper during the updateWeights phase.
	
	double  low_dist, curr_dist;
	for(auto& observ: observations){
		low_dist  = numeric_limits<float>::max();
		for(const auto& pred: predicted){
		  curr_dist = dist(observ.x, observ.y, pred.x, pred.y);
		  if( low_dist  > curr_dist){
			low_dist  = curr_dist;
			observ.id = pred.id;
		  }
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	// more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	// according to the MAP'S coordinate system. You will need to transform between the two systems.
	// Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	// The following is a good resource for the theory:
	// https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	// and the following is a good resource for the actual equation to implement (look at equation 
	// 3.33
	// http://planning.cs.uiuc.edu/node99.html
	
	//Constant used for the new weights
	const double std0 = std_landmark[0];
	const double std1 = std_landmark[1];
	const double pstd0 = 2 * pow(std0, 2);
	const double pstd1 = 2 * pow(std1, 2);
	const double mpistd = 2 * M_PI * std0 * std1;
	double curr_dist, cos_theta, sin_theta, x_term, y_term, w;

	// loop through all particles
	for(auto& partic: particles){
		partic.weight = 1.0;

		// step 1: Filter valid landmarks which are in the sensor_range
		vector<LandmarkObs> predictions;
		for(const auto& list_lmark: map_landmarks.landmark_list){
		  curr_dist = dist(partic.x, partic.y, list_lmark.x_f, list_lmark.y_f);
		  if( curr_dist < sensor_range){ 
			predictions.push_back(LandmarkObs{list_lmark.id_i, list_lmark.x_f, list_lmark.y_f});
		  }
		}

		// step 2: Convert coordinates from vehicle to map
		vector<LandmarkObs> observations_map;
		cos_theta = cos(partic.theta);
		sin_theta = sin(partic.theta);

		for(const auto& obs: observations){
		  LandmarkObs tmp;
		  tmp.x = obs.x * cos_theta - obs.y * sin_theta + partic.x;
		  tmp.y = obs.x * sin_theta + obs.y * cos_theta + partic.y;
		  observations_map.push_back(tmp);
		}

		// step 3: Associate landmark index for each observation
		dataAssociation(predictions, observations_map);

		// step 4: Accumulate weights for particle
		for(const auto& obs_m: observations_map){

		  Map::single_landmark_s landmark = map_landmarks.landmark_list.at(obs_m.id-1);
		  x_term = pow(obs_m.x - landmark.x_f, 2) / pstd0;
		  y_term = pow(obs_m.y - landmark.y_f, 2) / pstd1;
		  w = exp(-(x_term + y_term)) / mpistd;
		  partic.weight *=  w;
		}
		weights.push_back(partic.weight);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	// http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// generate distribution accordeviceing to weights

	// Create a generator for random particle
	default_random_engine gen;
	//random_device rdevice;
    //mt19937 gen(rdevice());
  
	discrete_distribution<> dist(weights.begin(), weights.end());

	//Generate random re-sample particle
	vector<Particle> resampled_particles;
	resampled_particles.resize(num_particles);

	// Re-sample base on weight
	for(int i=0; i<num_particles; i++){
		int current_index = dist(gen);
		resampled_particles[i] = particles[current_index];
	}

	// Set re sample particle
	particles = resampled_particles;

	// Clear particle weight
	weights.clear();
}

Particle ParticleFilter::SetAssociations(Particle particle, vector<int> associations, vector<double> sense_x, vector<double> sense_y)
{
	// particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
	
	// Clear the previous associations
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