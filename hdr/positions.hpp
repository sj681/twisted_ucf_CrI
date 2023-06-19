#ifndef POSITION_HPP
#define POSITION_HPP

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>


   extern double twist_angle;

   extern int number_of_unit_cells_x;
   extern int number_of_unit_cells_y;

   extern double a0x;
   extern double a0y;
   extern double c0;

   extern int num_atoms;
   extern int num_nm_atoms;
   extern double twist_loction;
   extern int total_atoms;
   extern int total_nm_atoms;
   extern int num_above_atoms;
   extern int num_below_atoms;

   class spin {
      public:
         double x;
         double y;
         double z;
         double S;
         int id;
         double dx;
         double dy;
   };

   extern std::vector < spin > atom;
   extern std::vector < spin > nm_atom;
   extern std::vector < spin > below_twist;
   extern std::vector < spin > above_twist;
   extern std::vector < spin > all_nm_atoms;

   void read_in_atoms(std::string filename, int n_atoms, std::vector <spin > &atom2);
   void read_in_exchange(std::string filename);
   void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz);

   void print_header();
   void create_atom_list(int n, bool magnetic, std::vector < spin > uc_atom_list, std::vector < spin > &above, std::vector < spin > &below, int &total, std::string filename);

#endif
