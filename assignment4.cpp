/*Description: This code creates lepton objects and manipulates them.
This includes using overloaded operators to sum vectors, finding
the dot product of the momentum of two leptons, as well as copying and 
moving data from one lepton object to another via various means. 
Author: Hassan Hashmi
Date: 16/04/2024*/ 

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip> 
#include <vector>

const double speed_of_light = 1; //in natural units 

class lepton
{
private:
  std::string particle_type;
  double rest_mass;          
  int charge;                
  std::vector<double>* momentum_vector {nullptr};          
public:
  lepton() = default;
  lepton(const std::string &input_particle_type, double input_mass, int input_charge,
         double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z);
  ~lepton();
  void set_particle_type(const std::string &input_particle_type) {particle_type = input_particle_type;}
  void set_rest_mass(double input_mass) {rest_mass = input_mass;}
  void set_charge(int input_charge) 
  {
    if(input_charge == 1 || input_charge == -1)
      charge = input_charge;
    else
      std::cout<<"Invalid charge. Must be 1 (for particles) or -1 (for anti-particles)"<<std::endl; 
  }
  void set_momentum_magnitude(double input_energy) 
  {
    if(input_energy < 0) 
      std::cout<<"The energy value cannot be negative. "<<std::endl;
    else
      (*momentum_vector)[0] = input_energy / speed_of_light;
  }
  void set_momentum_x(double input_momentum_x) {(*momentum_vector)[1] = input_momentum_x;}
  void set_momentum_y(double input_momentum_y) {(*momentum_vector)[2] = input_momentum_y;}
  void set_momentum_z(double input_momentum_z) {(*momentum_vector)[3] = input_momentum_z;}
  std::string get_particle_type() const {return particle_type;}
  double get_rest_mass() const {return rest_mass;}
  int get_charge() const {return charge;}
  double get_momentum_magnitude() const {return (*momentum_vector)[0];}
  double get_momentum_x() const {return (*momentum_vector)[1];}
  double get_momentum_y() const {return (*momentum_vector)[2];}
  double get_momentum_z() const {return (*momentum_vector)[3];}
  friend void print_lepton_data(const lepton &input_lepton);  //declare friend function to access data members outside class
  std::vector<double> operator+(const lepton &particle) const; 
  double dotProduct(const lepton &particle) const; 
  lepton(const lepton &particle);     //copy constructor 
  lepton &operator=(const lepton&);  //copy assingment operator 
  lepton(lepton&&);                 //move constructor 
  lepton &operator =(lepton&&);    //move assingment operator 
};

lepton::lepton(const std::string &input_particle_type, double input_mass, int input_charge,
               double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z):
  particle_type{input_particle_type}, rest_mass{input_mass}, charge{input_charge} 
{
  momentum_vector = new std::vector<double>(4); 
  set_momentum_magnitude(input_energy);
  set_momentum_x(input_momentum_x);
  set_momentum_y(input_momentum_y);
  set_momentum_z(input_momentum_z);
};

lepton::~lepton() 
{   
  if(particle_type == "")
    std::cout<<std::endl; 
  else 
    std::cout<<"\nCalling destructor for "<<particle_type<< "."<<std::endl;
  delete momentum_vector; 
}

double lepton::dotProduct(const lepton &particle) const 
{
  double dot_product = (*momentum_vector)[0] * (*particle.momentum_vector)[0] 
                      + (*momentum_vector)[1] * (*particle.momentum_vector)[1] 
                      + (*momentum_vector)[2] * (*particle.momentum_vector)[2] 
                      + (*momentum_vector)[3] * (*particle.momentum_vector)[3];
  return dot_product;
}

std::vector<double> lepton::operator+(const lepton &particle) const 
{
  std::vector<double> total_vector(4, 0.0); 
  total_vector[0] = (*momentum_vector)[0] + (*particle.momentum_vector)[0];
  total_vector[1] = (*momentum_vector)[1] + (*particle.momentum_vector)[1];
  total_vector[2] = (*momentum_vector)[2] + (*particle.momentum_vector)[2];
  total_vector[3] = (*momentum_vector)[3] + (*particle.momentum_vector)[3];
  return total_vector;
}

lepton::lepton(const lepton &particle_to_copy) 
{
  std::cout<<"\nCalling copy constructor. "<<std::endl;
  particle_type = particle_to_copy.particle_type; 
  rest_mass = particle_to_copy.rest_mass;         
  charge = particle_to_copy.charge;              
  momentum_vector = new std::vector<double>(*particle_to_copy.momentum_vector);
}

lepton& lepton::operator=(const lepton &particle_to_copy) 
{
  std::cout<<"\nCalling copy assignment operator. "<<std::endl;
  if(&particle_to_copy == this) 
    return *this;
  delete momentum_vector; 
  momentum_vector = nullptr;
  particle_type = particle_to_copy.particle_type;
  rest_mass = particle_to_copy.rest_mass;
  charge = particle_to_copy.charge;
  momentum_vector = new std::vector<double>(*particle_to_copy.momentum_vector);
  return *this;
}

lepton::lepton(lepton &&target_particle)
{
  std::cout<<"\nCalling move constructor. "<<std::endl; 
  particle_type = target_particle.particle_type;
  rest_mass = target_particle.rest_mass;
  charge = target_particle.charge;
  momentum_vector = target_particle.momentum_vector;
  target_particle.particle_type = "";
  target_particle.rest_mass = 0;
  target_particle.charge = 0; 
  target_particle.momentum_vector = nullptr;
}

lepton &lepton::operator=(lepton &&target_particle)
{
  std::cout<<"\nCalling move assignment operator. "<<std::endl; 
  std::swap(particle_type, target_particle.particle_type);
  std::swap(rest_mass, target_particle.rest_mass);
  std::swap(charge, target_particle.charge);
  std::swap(momentum_vector, target_particle.momentum_vector);
  return *this; 
}

void print_lepton_data(const lepton &input_lepton) 
{
  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"\nType of lepton: "<<input_lepton.particle_type<<std::endl;
  std::cout<<"Rest Mass (MeV): "<<input_lepton.rest_mass<<std::endl;
  std::cout<<"Charge: "<<input_lepton.charge<<std::endl;
  std::cout<<"E/c: "<<(*input_lepton.momentum_vector)[0]<<std::endl; 
  std::cout<<"Momentum in x direction: "<<(*input_lepton.momentum_vector)[1]<<std::endl; 
  std::cout<<"Momentum in y direction: "<<(*input_lepton.momentum_vector)[2]<<std::endl; 
  std::cout<<"Momentum in z direction: "<<(*input_lepton.momentum_vector)[3]<<std::endl; 
}

int main() 
{
  std::vector<lepton> lepton_information_vector = 
    {
      lepton("electron", 0.511, -1, 0, 0, 0, 0), 
      lepton("electron", 0.511, -1, 0, 0, 0, 0),
      lepton("muon", 105.7, -1, 0, 0, 0, 0),
      lepton("muon", 105.7, -1, 0, 0, 0, 0),
      lepton("muon", 105.7, -1, 0, 0, 0, 0),
      lepton("muon", 105.7, -1, 0, 0, 0, 0), 
      lepton("anti-electron", 0.511, 1, 0, 0, 0, 0),
      lepton("anti-muon", 105.7, 1, 0, 0, 0, 0)      
    };
  std::vector<double> total_momentum = lepton_information_vector[0] + lepton_information_vector[1];
  std::cout<<"\nSum of each component of the momentum vectors of two electrons - "<<std::endl;
  std::cout<<"Total E/c: "<<total_momentum[0]<<std::endl;
  std::cout<<"Momentum in the x-direction: "<<total_momentum[1]<<std::endl;
  std::cout<<"Momentum in the y-direction: "<<total_momentum[2]<<std::endl;
  std::cout<<"Momentum in the z-direction: "<<total_momentum[3]<<std::endl;
  double dot_product = lepton_information_vector[2].dotProduct(lepton_information_vector[3]);
  std::cout<<"\nDot product of the first two muons: "<<dot_product<<std::endl; 
  lepton new_muon (lepton_information_vector[2]);
  std::cout<<"Deep copying first muon data to new muon (via copy constructor) -  "<<std::endl;
  print_lepton_data(new_muon);
  lepton new_electron; 
  new_electron = lepton_information_vector[0]; 
  std::cout<<"Deep copying first electron data to a new electron (via assignment operator) -  "<<std::endl;
  print_lepton_data (new_electron); 
  lepton new_antielectron(std::move(lepton_information_vector[6]));
  std::cout<<"Moving antielectron data from old antielectron object to new antielectron object - "<<std::endl; 
  print_lepton_data(new_antielectron);
  lepton new_antimuon;
  new_antimuon = std::move (lepton_information_vector[7]);
  std::cout<<"Moving antimuon data from old antimuon object to new antimuon object - "<<std::endl; 
  print_lepton_data(new_antimuon);
  return 0;
}