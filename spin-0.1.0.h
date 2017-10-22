#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <cstdlib>

#define NSITES 243
#define ETA 0.04
#define NITR 2000000
#define KB 1
#define TEMP 0.0000115

class SpinConfig{
	public:
		int randsite;
		string spinfilename;
		FILE * spinfile;
		double crossproduct[3],spinset[NSITES][3],oldspin[3],newspin[3];
		double divide(int,int);
		double dot(double (*), double (*));
		void cross(double (*), double (*));
		double create_rand(double,double);
		void create_initial_spinset();
		void import_spinset();
		void set_old_spin();
		void set_new_spin();
		void update_spinset();
		void print_spinset_to_file();
		void print_spinset_to_file2();
};
