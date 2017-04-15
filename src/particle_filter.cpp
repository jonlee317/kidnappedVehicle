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

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	num_particles = 1000;
	Particle myParticle;
	std::default_random_engine gen;

	for (int i=0; i<num_particles; i++) {
		myParticle.id = i;
		myParticle.x = dist_x(gen);
		myParticle.y = dist_y(gen);
		myParticle.theta = dist_theta(gen);
		myParticle.weight = 1;
		particles.push_back(myParticle);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Creating noise generators
	std::default_random_engine gen;
	std::normal_distribution<double> dist_xpos(0,std_pos[0]);
	std::normal_distribution<double> dist_ypos(0,std_pos[1]);
	std::normal_distribution<double> dist_thetapos(0,std_pos[2]);

	for (int i=0; i<particles.size(); i++) {
		if (yaw_rate == 0) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}
		else {
			particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
			particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;
		}
		// Adding noise
		particles[i].x += dist_xpos(gen);
		particles[i].y += dist_ypos(gen);
		particles[i].theta += dist_thetapos(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, std::vector<LandmarkObs>& the_chosen) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	double distance_min;
	double distance_cur;
	double xdiff;
	double ydiff;

	// for each observation we check the distance of each landmark in range (a.k.a. predicted)

	for (int i = 0; i<observations.size(); i++) {
		// We first initialize the distance to something large so that any initial distance will be smaller
		distance_min = 10000;
		for (int j=0; j<predicted.size(); j++) {
			// calculate the distance
			xdiff = observations[i].x - predicted[j].x;
			ydiff = observations[i].y - predicted[j].y;
			distance_cur = sqrt(xdiff*xdiff+ydiff*ydiff);

			if (distance_cur < distance_min) {
				// Updating the current distance minimum
				distance_min = distance_cur;
				// Adding this measurement to my chosen landmarks
				the_chosen[i] = predicted[j];
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	// Initializing some useful variables
	double multi_gauss;
	double my_pi = 3.14159;

	double sigx = std_landmark[0];
	double sigy = std_landmark[1];
	double sigx_sq = sigx*sigx;
	double sigy_sq = sigy*sigy;

	// initializing weights
	double my_weight=1;
	weights.clear();

	// initializing my vectors used
	std::vector<LandmarkObs> in_range;
	std::vector<LandmarkObs> observed_convert;
	std::vector<LandmarkObs> chosen_landmarks;

	// now that we are done initializing
	// Let's loop through each particle
	for (int i = 0; i < particles.size(); i++) {
		// Find the landmarks in range aka predicted landmarks
		for (int l = 0; l < map_landmarks.landmark_list.size(); l++) {
			double px_lx = particles[i].x-map_landmarks.landmark_list[l].x_f;
			double py_ly = particles[i].y-map_landmarks.landmark_list[l].y_f;
			double curr_dist = sqrt(px_lx*px_lx+py_ly*py_ly);
			// Checking if current distance is within the sensor range and add to vector in_range
			if (curr_dist <= sensor_range) {
				LandmarkObs pre_land;
				pre_land.id = map_landmarks.landmark_list[l].id_i;
				pre_land.x = map_landmarks.landmark_list[l].x_f;
				pre_land.y = map_landmarks.landmark_list[l].y_f;
				in_range.push_back(pre_land);
			}
		}

		// Convert each of the observed landmarks for each particle i and add to vector observed_convert
		for (int o = 0; o < observations.size(); o++) {
			LandmarkObs obs_conv;
			obs_conv.x = particles[i].x + (observations[o].x*cos(particles[i].theta)-observations[o].y*sin(particles[i].theta));
			obs_conv.y = particles[i].y + (observations[o].y*cos(particles[i].theta)+observations[o].x*sin(particles[i].theta));
			observed_convert.push_back(obs_conv);
		}

		// initializing chosen landmarks with number of observed converted
		// this vector will be later updated in the dataAssociation method
		chosen_landmarks = observed_convert;

		// setting chosen landmarks to be the min distance between the converted observed and the landmarks in range
		dataAssociation(in_range, observed_convert, chosen_landmarks);

		// Finding the deviation between the converted observations and the chosen landmarks
		for (int k=0; k<observed_convert.size(); k++) {
			double x_dev = chosen_landmarks[k].x-observed_convert[k].x;
			double x_dev_sq = x_dev*x_dev;
			double y_dev = chosen_landmarks[k].y-observed_convert[k].y;
			double y_dev_sq = y_dev*y_dev;

			// calculate the gaussian probablity and multiply by current weight of particle
			multi_gauss = (0.5/(my_pi*sigx*sigy))*exp(-((0.5*x_dev_sq/sigx_sq) + (0.5*y_dev_sq/sigy_sq)));
			my_weight *= multi_gauss;
		}

		// update particle weight and add weight to weights vector
		particles[i].weight = my_weight;
		weights.push_back(my_weight);

		// Done with current particle!
		// Reinitializing weight and clearing vectors for next particle
		my_weight=1;
		in_range.clear();
		observed_convert.clear();
		chosen_landmarks.clear();
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   [http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution]

	std::default_random_engine gen;
	// This discrete distribution will generate a random integer according to the weight
	// larger weight will have higher chance of being chosen!
	// This weight is provided by the vector weights
  std::discrete_distribution<> distribution(weights.begin(), weights.end());

	// Creating a bag of particles called resampled_particles
	std::vector<Particle> resampled_particles;
	resampled_particles.clear();
  for (int i=0; i<num_particles; i++) {
    int new_index = distribution(gen);
    resampled_particles.push_back(particles[new_index]);
  }

	// redefining our original set of particles with this new bag of particles
	particles = resampled_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
