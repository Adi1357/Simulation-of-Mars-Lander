// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include <iostream>
#include <fstream>
#include "lander.h"
//#include "lander_graphics.cpp"
#include <vector>
#include <math.h>
using namespace std;

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
  // INSERT YOUR CODE HERE
    double h,e, P_out,K_h,K_p,delta, desired_rate;
    static vector<double> h_list, v_e_r_list, desired_rate_list;
    static double prev_sim_time;

    delta = (GRAVITY * MARS_MASS * (UNLOADED_LANDER_MASS + fuel * FUEL_DENSITY * FUEL_CAPACITY) / position.abs2()) / MAX_THRUST;
    K_h = 0.02;
    K_p = 0.4;

    h = position.abs() - MARS_RADIUS;
    e = -(0.5 + K_h*h + velocity*position.norm());
    desired_rate = -(0.5 + K_h * h);

    P_out = K_p * e;

    if (P_out <= -delta){
        throttle = 0;
    }
    else if ((-delta < P_out) && (P_out < 1 - delta)){
        throttle = delta + P_out;
    }
    else {
        throttle = 1;
    }

    h_list.push_back(h);
    v_e_r_list.push_back(velocity * position.norm());
    desired_rate_list.push_back(desired_rate);


    

    if (h <= 10) {
        ofstream fout;
        fout.open("verlet_mars.txt");
        if (fout) { // file opened successfully
            for (int i = 0; i < h_list.size(); i = i + 1) {
                fout << h_list[i] << ' ' << v_e_r_list[i] << ' '<< desired_rate_list[i] << endl;
                
            }
        }
        else { // file did not open successfully
            cout << "Could not open trajectory file for writing" << endl;
        }
    }

    prev_sim_time = simulation_time;
}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{
  // INSERT YOUR CODE HERE
  static vector3d previous_position;
  vector3d next_position, acceleration, thrust, drag_lander, drag, weight, net_force, position_current;
  vector<double> t_list;
  //vector<vector3d> pos_list, vel_list, accel_list;
  double m, G, M, r0_mag, dt, t_max, t, density;
  
  // mass, spring constant, initial position and velocity
  m = UNLOADED_LANDER_MASS + fuel*FUEL_CAPACITY* FUEL_DENSITY;
  G = GRAVITY;
  M = MARS_MASS;
  r0_mag = MARS_RADIUS;
  dt = delta_t;
  t_max = simulation_time;
  t = 0;
  

  thrust = thrust_wrt_world();
  density = atmospheric_density(position);
  drag_lander = -0.5 * density * DRAG_COEF_LANDER * M_PI* (LANDER_SIZE * LANDER_SIZE) * (velocity.abs2()) * velocity.norm();

  if (parachute_status == DEPLOYED) {
      drag = -0.5 * density * DRAG_COEF_CHUTE * (5*(2*LANDER_SIZE)* (2 * LANDER_SIZE)) * (velocity.abs2()) * velocity.norm() + drag_lander;
  }
  else drag = drag_lander;

  weight = -1 * G * M * m * position.norm() / (position.abs2());

  net_force = weight + drag + thrust;


  acceleration = net_force / m;
  if (simulation_time == 0) {
      previous_position = -(velocity * dt - position);
  }

  else {
      previous_position = previous_position;
  }

  position_current = position;
  position = 2 * position - previous_position + dt * dt * acceleration;
  velocity = (position - position_current) / dt;
  previous_position = position_current;


/*
  if (position.abs() > MARS_RADIUS) {
      acceleration = net_force / m;
      if (simulation_time == 0) {
          previous_position = -(velocity * dt - position);
      }
      
      else {
          previous_position = previous_position;
      }

      position_current = position;
      position = 2 * position - previous_position + dt * dt * acceleration;
      velocity = (position - position_current) / dt;
      previous_position = position_current;

  }

  else {
      previous_position = position;
      velocity.x  = 0.0;
      velocity.y = 0.0;
      velocity.z = 0.0;
  }


      // Verlet integration
  for (int i = 0; i <= t_max / dt; i++, t += dt) {

      // append current state to trajectories
      t_list.push_back(t);
      pos_list.push_back(position);
      vel_list.push_back(velocity);
      accel_list.push_back(acceleration);


      thrust = thrust_wrt_world();
      density = atmospheric_density(position);
      drag_lander = -0.5 * density * DRAG_COEF_LANDER * (LANDER_SIZE * LANDER_SIZE) * (velocity.abs2()) * velocity.norm();

      if (parachute_status == DEPLOYED) {
          drag = -0.5 * density * DRAG_COEF_CHUTE * (LANDER_SIZE * LANDER_SIZE) * (velocity.abs2()) * velocity.norm() + drag_lander;
      }
      else drag = drag_lander;

      weight = -1 * G * M * m * position / (position.abs2());

      net_force = weight + drag + thrust;


      position_current = position;

      if (position.abs() > MARS_RADIUS) {
          acceleration = net_force / m;
          position = 2 * position - previous_position + dt * dt * acceleration;
          velocity = (position - position_current) / dt;

          previous_position = position_current;
      }

      else {
          previous_position = position;
          position_current = position;
      }
      }

*/
      // Here we can apply an autopilot to adjust the thrust, parachute and attitude
      if (autopilot_enabled) autopilot();

      // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
      if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 6:
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;

  }
}
