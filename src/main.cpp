#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "exchange.hpp"
#include "initialise.hpp"


int main(){

   //in A
   system_size_x = 1000;
   system_size_y = 1000;
   //in number of unit cells
   number_of_unit_cells_z = 1;
   //in degrees
   twist_angle =1.41;

   initialise_variables();

   twist_loction = 0.25*system_size_z-0.01;;

   read_in_atoms("files/atom_list_aa", num_atoms, atom);
   read_in_atoms("files/nm_atoms", num_nm_atoms, nm_atom);
   read_in_exchange("files/interpolated_array");
   read_in_dmi("files/interpolated_Dij_intra", Dx_intra, Dy_intra, Dz_intra);
   read_in_dmi("files/interpolated_Dij_inter", Dx_inter, Dy_inter, Dz_inter);

   create_atom_list(num_atoms, true, atom, above_twist, below_twist, total_atoms, "atom_positions.ucf");
   create_atom_list(num_nm_atoms, false, nm_atom, all_nm_atoms, all_nm_atoms, total_nm_atoms, "nm_atom_positions.ucf");

    print_header();
    calc_in_plane_exchange(above_twist);
    calc_in_plane_exchange(below_twist);
    calc_out_of_plane_exchange(above_twist,below_twist);
    print_interaction_header();
 }
