// Projectile
// 
// Need 'Eigen' library
//
#include <iostream>
#include <fstream>
#include <Eigen/Core>

int main(int argc, char **argv) {
  //// Ask user to input initial values
  ///////////////////////
  //// YOUR CODE HERE
  ///////////////////////

  ///// Define velocity and Compute Initial Velocity
  //////////////////////
  //// YOUR CODE HERE
  ///////////////////////

  // Initialize
  const double g = -9.8; // m/s^2
  const double drag = 1;
  const double weight = 1;
  Eigen::Vector2d position;
  position.setZero();
  Eigen::Vector2d acceleration;
  acceleration << 0.0, g;

  // Simulation
  std::ofstream fout;
  fout.open("hw06.txt");
  const double dt = 0.01;
  double simulationTime = 0.0;
  while(position[1] >= 0.0){
    fout << simulationTime << " " << position[0] << " " << position[1] << " " << velocity[0] << " " << velocity[1] << std::endl;
    /// Update
    ////////////////////
    /// YOUR CODE HERE
    ////////////////////
  }
  fout.close();
  return 0;
}