#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "initialise.hpp"
#include "exchange.hpp"

//set nearest neighbour  distances (in plane nn 1,2,3)
double nn_dist_1 = a0x*pow(3,0.5)/3;
double nn_dist_2 = a0x;
double nn_dist_3 = a0x*pow(1.3333333333,0.5);

//Set exchange interaction values and associated constants
double eVtoJ = 1.602176634e-19;
double J_constant = 1.0*2.25*eVtoJ/1000;
double J_intra_1=2.68*J_constant;
double J_intra_2=0.71*J_constant;
double J_intra_3=-0.01*J_constant;

//set the initial jumber of interactions to zero for counter
int number_of_interactions = 0;

//initialise arrays to store exchage interactions
std::vector < std::vector < double > > Jint;
std::vector < std::vector < double > > Dx_inter;
std::vector < std::vector < double > > Dy_inter;
std::vector < std::vector < double > > Dz_inter;
std::vector < std::vector < double > > Dx_intra;
std::vector < std::vector < double > > Dy_intra;
std::vector < std::vector < double > > Dz_intra;


std::vector <double > crossProduct(std::vector <double >A, std::vector <double > B){
   std::vector <double > P(3,0.0);
    P[0] = A[1] * B[2] - A[2] * B[1];
    P[1] = A[2] * B[0] - A[0] * B[2];
    P[2] = A[0] * B[1] - A[1] * B[0];
    return P;
}

void print_interaction_header(){
   std::ofstream outfile3 ("header_interactions.ucf");
   outfile3 << number_of_interactions <<  "\t" << "tensorial"<< std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
}


std::vector <double >  calculate_dmi_vector(double x_i, double y_i, double z_i,double x_j, double y_j, double z_j, double dx, double dy, double dz){

   std::vector <double > Duv(3,0.0);
    for (int nm = 0; nm < total_nm_atoms; nm ++){
       double x_nm = all_nm_atoms[nm].x;
       double y_nm = all_nm_atoms[nm].y;
       double z_nm = all_nm_atoms[nm].z;
       double dl1 = (x_nm - x_i)*(x_nm - x_i) + (y_nm - y_i)*(y_nm - y_i) + (z_nm - z_i)*(z_nm - z_i);
       double dl2 = (x_nm - x_j)*(x_nm - x_j) + (y_nm - y_j)*(y_nm - y_j) + (z_nm - z_j)*(z_nm - z_j);
      // std::cout << dl1 << '\t' << dl2 << "\t" << nn_dist_1*nn_dist_1 << std::endl;
       if (dl1  <nn_dist_2*nn_dist_2 && dl2 < nn_dist_2*nn_dist_2){
          std::vector <double > uij(3,0.0);
          std::vector <double > zij(3,0.0);
          std::vector <double > Dij(3,0.0);
         // std::cout << "a" << std::endl;
         double norm = sqrt(dx*dx + dy*dy + dz*dz);
         uij[0] = dx/norm;
         uij[1] = dy/norm;
         uij[2] = dz/norm;
         double midx = ((x_i + x_j)/2);
         double midy = ((y_i + y_j)/2);
         double midz = ((z_i + z_j)/2);
         //std::cout << atom << "\t" << atomi << "\t"<< midx << "\t" <<  midy << "\t" <<  midz/c0 << "\t" <<  x_i << "\t" <<  y_i << "\t" <<  z_i << "\t" << x_j << "\t" <<  y_j << "\t" <<  z_j<< "\t" <<  std::endl;

         zij[0] = midx - x_nm;
         zij[1] = midy - y_nm;
         zij[2] = midz - z_nm;
         double norm2 = sqrt(zij[0]*zij[0] + zij[1]*zij[1] + zij[2]*zij[2]);
         zij[0] = zij[0]/norm2;
         zij[1] = zij[1]/norm2;
         zij[2] = zij[2]/norm2;
         Dij = crossProduct(uij,zij);
         Duv[0] = Duv[0] + Dij[0];
         Duv[1] = Duv[1] + Dij[1];
         Duv[2] = Duv[2] + Dij[2];
         //std::cout << uij[0] << "\t" <<  uij[1] << "\t" <<  uij[2] << "\t" <<  zij[0] << "\t" <<  zij[1] << "\t" <<  zij[2] << "\t" <<  Dij[0] << "\t" <<  Dij[1] << "\t" <<  Dij[2] << "\t" <<  std::endl;
         //std::cout << x_i << "\t" << atomi << "\t" << nm <<  "\t" << Dij[0] << '\t' << Dij[1] << "\t" << Dij[2] << std::endl;
         //std::cout << atom << "\t" << atomi << "\t" << nm <<  "\t" << sqrt(dl1) << '\t' << sqrt(dl2) << std::endl;
      }
   }
   double Dl = sqrt(Duv[0]*Duv[0] + Duv[1]*Duv[1] + Duv[2]*Duv[2]);
   Duv[0] = Duv[0]/Dl;
   Duv[1] = Duv[1]/Dl;
   Duv[2] = Duv[2]/Dl;
      return Duv;

}

void calc_in_plane_exchange(std::vector < spin > atom_list_1){

   int num_atoms_1 = atom_list_1.size();
   std::cout << num_atoms_1 <<std::endl;

   for (int atom_i = 0; atom_i < num_atoms_1; atom_i ++){

      double x_i = atom_list_1[atom_i].x;
      double y_i = atom_list_1[atom_i].y;
      double z_i = atom_list_1[atom_i].z;
      int id_i = atom_list_1[atom_i].id;
      //std::cout << x_i << '\t' << y_i << "\t" << z_i << "\t" << id_i << std::endl;
      int nn = 0;
      int nnn = 0;
      int nnnn = 0;
      for (int atom_j = 0; atom_j < num_atoms_1; atom_j ++){
         int id_j = atom_list_1[atom_j].id;

         if (id_i != id_j){
            double x_j = atom_list_1[atom_j].x;
            double y_j = atom_list_1[atom_j].y;
            double z_j = atom_list_1[atom_j].z;
            double dx = x_i - x_j;
            double dy = y_i - y_j;
            double dz = z_i - z_j;
            double dL = sqrt(dx*dx + dy*dy+dz*dz);

            if (dL < nn_dist_1 + 0.01){
            //            std::cout << "JIJ" << "\t" <<x_j << '\t' << y_j << "\t" << z_j << "\t" << id_j << std::endl;
               outfile4 << number_of_interactions <<  "\t" << id_i << '\t' << id_j << " 0 0 0 "<<J_intra_1 << "\t" << 0 << "\t" << 0<< "\t" << 0 << "\t" << J_intra_1 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << J_intra_1 <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
               nn++;
               number_of_interactions++;
            }
            else if (dL < nn_dist_2+0.01){
               outfile4 << number_of_interactions <<  "\t" << id_i << '\t' << id_j << " 0 0 0 "<<J_intra_2 << "\t" << 0 << "\t" << 0<< "\t" << 0 << "\t" << J_intra_2 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << J_intra_2 <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
               nnn++ ;
               number_of_interactions++;
            }
            else if (dL < nn_dist_3 + 0.01){
               outfile4 << number_of_interactions <<  "\t" << id_i << '\t' << id_j << " 0 0 0 "<<J_intra_3 << "\t" << 0 << "\t" << 0<< "\t" << 0 << "\t" << J_intra_3 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << J_intra_3 <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
               nnnn++;
               number_of_interactions++;
            }
         }
      }
   //   std::cout << x_i << '\t' << y_i << '\t' << z_i << "\t" << nn << '\t' << nnn << '\t' << nnnn << std::endl;
   }
}


void calc_out_of_plane_exchange(std::vector < spin > atom_list_1,std::vector < spin > atom_list_2){

   int num_atoms_1 = atom_list_1.size();
   int num_atoms_2 = atom_list_2.size();

   for (int atom_i = 0; atom_i < num_atoms_1; atom_i ++){

      double x_i = atom_list_1[atom_i].x;
      double y_i = atom_list_1[atom_i].y;
      double z_i = atom_list_1[atom_i].z;
      double id_i = atom_list_1[atom_i].id;
      double changex = atom_list_1[atom_i].dx;
      double changey = atom_list_1[atom_i].dy;
      //std::cout << x_i << "\t" << y_i << "\t" << changex << '\t' << changey << '\t' << Jint[changex][changey] << std::endl;
      double Jj = Jint[changex][changey]*J_constant;
      double Dx = Dx_inter[changex][changey]*J_constant;
      double Dy = Dy_inter[changex][changey]*J_constant;
      double Dz = Dz_inter[changex][changey]*J_constant;

      int nn = 0;
      for (int atom_j = 0; atom_j < num_atoms_2; atom_j ++){

         double x_j = atom_list_2[atom_j].x;
         double y_j = atom_list_2[atom_j].y;
         double z_j = atom_list_2[atom_j].z;
         double id_j = atom_list_2[atom_j].id;

         double dx = x_i - x_j;
         double dy = y_i - y_j;
         double dz = z_i - z_j;
         double dL = sqrt(dx*dx + dy*dy+dz*dz);
         //std::cout << dL << '\t' << nn_dist_3 <<std::endl;
         if (dL < nn_dist_3){
            nn++;
            //std::cout << id_i << std::endl;
            //std::cout << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx*Dij[0] << "\t" <<  Dy*Dij[1] << "\t" <<  Dz*Dij[2]<< "\t" <<  Jj << std::endl;
            //std::vector <double > Dij = calculate_dmi_vector(x_i, y_i,z_i, x_j,y_j,z_j,dx,dy,dz);
         //   std::cout << id_i << "\t" << id_j << '\t' << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx*Dij[0] << "\t" <<  Dy*Dij[1] << "\t" <<  Dz*Dij[2]<< "\t" <<  Jj << std::endl;//"\t" << "\t" << Dij[0] << "\t" << Dij[1] << "\t" << Dij[2] << std::endl;
            std::vector <double > Dij(3,0.0);
         // std::cout << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx << "\t" <<  Dy << "\t" <<  Dz<< "\t" <<  Dij[0] << '\t' << Dij[1] << '\t' << Dij[2]  << std::endl;


            outfile4 << number_of_interactions <<  "\t" << id_i << '\t' << id_j <<" 0 0 0 "<<Jj << "\t" << Dz*Dij[2] << "\t" << -Dy*Dij[1]<< "\t" << -Dz*Dij[2] << "\t" << Jj << "\t" << Dx*Dij[0] << "\t" << Dy*Dij[1] << "\t" << -Dx*Dij[0] << "\t" << Jj <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
            number_of_interactions++;

            //Dij = calculate_dmi_vector(x_j, y_j,z_j, x_i,y_i,z_i,-dx,-dy,-dz);
            outfile4 << number_of_interactions <<  "\t" << id_j << '\t' << id_i <<" 0 0 0 "<<Jj << "\t" << Dz*Dij[2] << "\t" << -Dy*Dij[1]<< "\t" << -Dz*Dij[2] << "\t" << Jj << "\t" << Dx*Dij[0] << "\t" << Dy*Dij[1] << "\t" << -Dx*Dij[0] << "\t" << Jj <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
            //std::cout << id_i << "\t" << id_j << '\t' << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx*Dij[0] << "\t" <<  Dy*Dij[1] << "\t" <<  Dz*Dij[2]<< "\t" <<  Jj << std::endl;//"\t" << "\t" << Dij[0] << "\t" << Dij[1] << "\t" << Dij[2] << std::endl;

            //outfile4 << number_of_interactions <<  "\t" << id_j << '\t' << id_i <<" 0 0 0 "<<Jj << "\t" << Dz*Dij[2] << "\t" << -Dy*Dij[1]<< "\t" << -Dz*Dij[2] << "\t" << Jj << "\t" << Dx*Dij[0] << "\t" << Dy*Dij[1] << "\t" << -Dx*Dij[0] << "\t" << Jj <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
            number_of_interactions++;
         }
      }
      //std::cout << nn << std::endl;
   }
}
