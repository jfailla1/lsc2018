// Projectile
// 
// Need 'Eigen' library
//
#include <iostream>
#include <fstream>
#include <Eigen/Core>
using namespace std;

int main(int argc, char **argv) {
  
	double speed;
	double angle;
	cout << "Enter Speed: ";
	cin >> speed;
	cout << "Enter Angle ";
	cin >> angle;

	Eigen::Vector2d velocity;
	velocity << (speed*cos(angle*M_PI/180)),(speed*sin(angle*M_PI/180));
	  
	// Initialize
  	const double g = -9.8; // m/s^2
 	const double drag = 1;
  	const double weight = 1;
  	Eigen::Vector2d position;
  	position.setZero();
  	Eigen::Vector2d acceleration;
  	acceleration << 0.0, g;
	Eigen::Vector2d G;
	G << 0.0, g;
	
  	// Simulation
  	std::ofstream fout;
  	fout.open("hw06.csv");
  	const double dt = 0.01;
  	double simulationTime = 0.0;
  	while(position[1] >= 0.0){
    		fout << simulationTime << "," << position[0] << "," << position[1] << "," << velocity[0] << "," << velocity[1] << std::endl;
	 	position += velocity * dt;
		velocity += acceleration * dt;
		simulationTime += dt;
		acceleration = G - (drag*velocity.norm()*velocity/weight); 	
	}
  	fout.close();
  	return 0;
}
