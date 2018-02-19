/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if (is_initialized) {
		return;
	}

	std::default_random_engine generator;
	std::normal_distribution<double> x_dist(x, std[0]);
	std::normal_distribution<double> y_dist(y, std[1]);
	std::normal_distribution<double> theta_dist(theta, std[2]);

	// Initialize particles and weights
	particles = vector<Particle>(num_particles);
	weights = vector<double>(num_particles);
	for (int i = 0; i < num_particles; i++) {
		particles[i].x = x_dist(generator);
		particles[i].y = y_dist(generator);
		particles[i].theta = theta_dist(generator);
		particles[i].weight = 1.0;
		weights[i] = 1.0;
	}

	is_initialized = true;
	return;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	std::default_random_engine generator;
	std::normal_distribution<double> x_noise(0.0, std_pos[0]);
	std::normal_distribution<double> y_noise(0.0, std_pos[1]);
	std::normal_distribution<double> theta_noise(0.0, std_pos[2]);

	for (int i = 0; i < particles.size(); i++) {
		double curt_x = particles[i].x;		
		double curt_y = particles[i].y;
		double curt_theta = particles[i].theta;

		double next_x = 0;
		double next_y = 0;
		double next_theta = 0;

		if (fabs(yaw_rate) < 1e-5) {
			next_x = curt_x + velocity * delta_t * cos(curt_theta);
			next_y = curt_y + velocity * delta_t * sin(curt_theta);
			next_theta = curt_theta;
		} else {
			next_x = curt_x + (velocity / yaw_rate) * (sin(curt_theta + yaw_rate * delta_t) - sin(curt_theta));
			next_y = curt_y + (velocity / yaw_rate) * (cos(curt_theta) - cos(curt_theta + yaw_rate * delta_t));
			next_theta = curt_theta + yaw_rate * delta_t;
		}

		particles[i].x = next_x + x_noise(generator);
		particles[i].y = next_y + y_noise(generator);
		particles[i].theta = next_theta + theta_noise(generator);
	}
}

void ParticleFilter::dataAssociation(const Map &map_landmarks, std::vector<LandmarkObs>& projected) {

	std::vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;

	for (int i = 0; i < projected.size(); i++) {
		double x_proj = projected[i].x;
		double y_proj = projected[i].y;

		double min_dist = std::numeric_limits<double>::max();
		int min_idx = -1;
		for (int j = 0; j < landmarks.size(); j++) {
			double x_landmark = landmarks[j].x_f;
			double y_landmark = landmarks[j].y_f;

			double curt_dist = dist(x_proj, y_proj, x_landmark, y_landmark);
			if (curt_dist < min_dist) {
				min_dist = curt_dist;
				min_idx = j;
			}
		}

		projected[i].id = min_idx;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
	const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	// calculate normalization term
	gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));

	for (int i = 0; i < particles.size(); i++) {
		double x_part = particles[i].x;
		double y_part = particles[i].y;
		double theta_part = particles[i].theta;

		weights[i] = 1.0; //particles[i].weight;

		// Project observations to map coordinates
		vector<LandmarkObs> projected(observations.size());
		for (int j = 0; j < observations.size(); j++) {
			double x_obs = observations[j].x;
			double y_obs = observations[j].y;

			double x_obs_map = x_part + (cos(theta_part) * x_obs) - (sin(theta_part) * y_obs);
			double y_obs_map = y_part + (sin(theta_part) * x_obs) + (cos(theta_part) * y_obs);

			projected[j].x = x_obs_map;
			projected[j].y = y_obs_map;
		}

		dataAssociation(map_landmarks, projected);

		for (int j = 0; j < projected.size(); j++) {
			Map::single_landmark_s associated_landmark = map_landmarks.landmark_list[projected[j].id];

			double exponent = pow((projected[j].x - associated_landmark.x_f), 2)/(2 * pow(std_landmark[0], 2)) 
			                + pow((projected[j].y - associated_landmark.y_f), 2)/(2 * pow(std_landmark[1], 2));

			weights[i] *= gauss_norm * pow(M_E, -exponent);
		}

	}

	// Normalize weights
	double sum = 0.0;
	for (int i = 0; i < weights.size(); i++) {
		sum += weights[i];
	}

	for (int i = 0; i < weights.size(); i++) {
		particles[i].weight = weights[i] / sum;
		weights[i] = particles[i].weight;
	}
}

void ParticleFilter::resample() {
	 std::random_device rd;
	 std::mt19937 gen(rd());
	 std::discrete_distribution<> distribution(weights.begin(), weights.end());

	 std::vector<Particle> resampled_particles = std::vector<Particle>(num_particles);
	 for(int i = 0; i < particles.size(); i++) {
		 resampled_particles[i] = particles[distribution(gen)];
		 weights[i] = particles[distribution(gen)].weight;
	 }

	 particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
