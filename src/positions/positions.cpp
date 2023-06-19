#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "initialise.hpp"



double a0x =7.0789999961853027;
double a0y= 12.260999679565430;
double c0 = 13.204;

int num_atoms = 8;
int num_nm_atoms = 24;

int number_of_unit_cells_x;
int number_of_unit_cells_y;

int num_above_atoms =0;
int num_below_atoms =0;


int total_atoms = 0;
int total_nm_atoms = 0;

void print_header(){

   std::ofstream outfile1 ("header.ucf");

   outfile1 << " #unit cell size " << std::endl;
   outfile1 << system_size_x << '\t' << system_size_y << '\t' << system_size_z << std::endl;
   outfile1 << " #unit cell vectors" << std::endl;
   outfile1 << " 1     0	   0" << std::endl;
   outfile1 << " 0     1	   0" << std::endl;
   outfile1 << " 0     0	   1" << std::endl;
   outfile1 << " #Atoms" << std::endl;
   outfile1 << total_atoms << '\t'	<< 2 << std::endl;

}

bool inside_system(double x, double y){

   if (x >=0 && x < system_size_x-0.001 && y >=0 && y< system_size_y -0.001) return true;
   else return false;
}


void create_atom_list(int n, bool magnetic, std::vector < spin > atom_list, std::vector < spin > &above, std::vector < spin > &below, int &total, std::string filename){

   double normalise_x = 100.0/(a0x*3.0/2.0);
   double normalise_y = 100.0/(a0x*sqrt(3/2.0));
   std::ofstream outfile2 (filename);
   for (int i = -number_of_unit_cells_x; i < 2*number_of_unit_cells_x; i++){
         for (int j = -number_of_unit_cells_y; j < 2*number_of_unit_cells_y; j++){
            for (int k = 0; k < number_of_unit_cells_z; k++){
               for (int atom_i = 0; atom_i < n; atom_i ++){

                  double x_j = atom_list[atom_i].x*a0x + i*a0x;
                  double y_j = atom_list[atom_i].y*a0y + j*a0y;
                  double z_j = atom_list[atom_i].z*c0 + k*c0;
                  if ( z_j > twist_loction){
                     double x_new = x_j*cos(twist_angle) - y_j*sin(twist_angle);
                     double y_new = y_j*cos(twist_angle) + x_j*sin(twist_angle);
                     if (inside_system(x_new, y_new)){
                        spin new_atom;
                        new_atom.x = x_new;
                        new_atom.y = y_new;
                        new_atom.z = z_j;
                        new_atom.S = 1;
                        new_atom.id = total;
                        if (magnetic){
                           double changex = fabs(x_new - x_j);
                           double changey = fabs(y_new - y_j);
                           //  std::cout << x_new << "\t" << y_new << "\t" << changex << "\t" << changey << '\t' << Jint[val_x2][val_y2] << std::endl;

                           int val_x2 = sqrt(changex*changex)*normalise_x;
                           int val_y2 = sqrt(changey*changey)*normalise_y;

                           while (val_x2 > 99){
                              val_x2 = val_x2 - 100;
                           }
                           while (val_y2 > 99){
                              val_y2 = val_y2 - 100;
                           }
                        //  std::cout << x_new << "\t" << y_new << "\t" << changex << "\t" << changey << std::endl;
                           new_atom.dx = val_x2;
                           new_atom.dy = val_y2;
                        }
                        above.push_back(new_atom);
                        outfile2 << total << "\t" << x_new/system_size_x << '\t' <<  y_new/system_size_y <<  "\t" << z_j/system_size_z << "\t" << 1 << "\t" << 0 << "\t" << 0 << std::endl;
                        total ++;
                        num_above_atoms++;
                     }
                  }
                  else if (inside_system(x_j, y_j)){
                     spin new_atom;
                     new_atom.x = x_j;
                     new_atom.y = y_j;
                     new_atom.z = z_j;
                     new_atom.id = total;
                     new_atom.S = 0;
                     below.push_back(new_atom);
                     outfile2 << total << "\t" << x_j/system_size_x << '\t' <<  y_j/system_size_y <<  "\t" << z_j/system_size_z << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
                     total ++;
                     num_below_atoms++;

              }
            }
         }
      }
   }
   std::cout << above.size() << "\t" << below.size() << "\t" << num_above_atoms << "\t" << num_below_atoms << std::endl;
}
