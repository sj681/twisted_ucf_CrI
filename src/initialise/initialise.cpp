#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "exchange.hpp"
#include "positions.hpp"
#include "initialise.hpp"

double system_size_x;
double system_size_y;
double system_size_z;
int number_of_unit_cells_z;

double twist_angle;
double twist_loction;

std::vector < spin > atom(num_atoms);
std::vector < spin > nm_atom(num_nm_atoms);
std::vector < spin > below_twist;
std::vector < spin > above_twist;
std::vector < spin > all_nm_atoms;


std::ofstream outfile4 ("interactions.ucf");

void resize_arrays(std::vector < std::vector < double > > &A, int sizex, int sizey){
   A.resize(sizex);
      for (int i = 0; i < sizex; ++i)
          A[i].resize(sizey, 0.0);
}

void initialise_variables(){

   number_of_unit_cells_x = system_size_x/a0x;
   number_of_unit_cells_y = system_size_y/a0y;
   system_size_z = number_of_unit_cells_z*c0;
   twist_angle = twist_angle*3.14159265359/180.0;
   resize_arrays(Jint, 101,101);
   resize_arrays(Dx_inter, 201,201);
   resize_arrays(Dy_inter, 201,201);
   resize_arrays(Dz_inter, 201,201);
   resize_arrays(Dx_intra, 201,201);
   resize_arrays(Dy_intra, 201,201);
   resize_arrays(Dz_intra, 201,201);
}
