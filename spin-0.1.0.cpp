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

using namespace std;

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

//Cast integer division to a double
double SpinConfig::divide(int x, int y)
{
	double quotient;
	quotient=static_cast<double>(x)/y;
	return quotient;
}

//Dot product of two vectors
double SpinConfig::dot(double firstvec[3],double secondvec[3])
{
	double product=0;
	for(int k=0;k<3;k++)
		product+=firstvec[k]*secondvec[k];
	return product;
}

//Define the cross product of two vectors
void SpinConfig::cross(double firstvec[3],double secondvec[3])
{
	crossproduct[0] = firstvec[1]*secondvec[2]-firstvec[2]*secondvec[1];
	crossproduct[1] = firstvec[2]*secondvec[0]-firstvec[0]*secondvec[2];
	crossproduct[2] = firstvec[0]*secondvec[1]-firstvec[1]*secondvec[0];
}

//Create random number in a given range
double SpinConfig::create_rand(double lowerbound,double upperbound)
{
	double randnum;
	double range;
	range=upperbound-lowerbound;
	randnum=range*divide(rand(),RAND_MAX)+lowerbound;
	return randnum;
}

//Create a random initial spin configuration
void SpinConfig::create_initial_spinset()
{
	double theta,phi;
	for(int i=0;i<NSITES;i++)
	{
		theta=acos(create_rand(0,1));
		phi=create_rand(0,2*M_PI);
		spinset[i][0]=sin(theta)*cos(phi);
		spinset[i][1]=sin(theta)*sin(phi);
		spinset[i][2]=cos(theta);
	}
}

//This method is used to import a spinset from file so that we can do simulated annealing
void SpinConfig::import_spinset()
{
	FILE * in_file;
        double d;

        in_file = fopen("/home/solidangle/FeCrAs_Model/Results/Monte_Carlo/Results/12x12/spinset-at-point5-run.dat","r");

        fscanf(in_file,"{");

        for(int i=0;i<NSITES-1;i++)
        {
                fscanf(in_file,"{");
                for(int j=0;j<2;j++)
                {
                        fscanf(in_file,"%lf,",&d);
                        spinset[i][j]=d;
                }
                fscanf(in_file,"%lf},",&d);
                spinset[i][2]=d;
        }

        fscanf(in_file,"{");
        for(int j=0;j<2;j++)
        {
                fscanf(in_file,"%lf,",&d);
                spinset[NSITES-1][j]=d;
        }
        fscanf(in_file,"%lf",&d);
        spinset[NSITES-1][2]=d;

        fclose(in_file);
}

//Use the randomly chosen index to set the "old" spin for this iteration
void SpinConfig::set_old_spin()
{
	randsite=rand()%NSITES;

	for(int k=0;k<3;k++)
		oldspin[k]=spinset[randsite][k];
}

//Flip the "old" spin into a new direction
void SpinConfig::set_new_spin()
{
	double theta,phi,alpha,khat[3];

	theta=acos(create_rand(-1,1));
	phi=create_rand(0,2*M_PI);
	alpha=ETA*create_rand(0,1);

	khat[0] = sin(theta)*cos(phi);
	khat[1] = sin(theta)*sin(phi);
	khat[2] = cos(theta);

	cross(khat,oldspin);

	for(int k=0;k<3;k++)
		newspin[k] = (cos(alpha)*oldspin[k])+(crossproduct[k]*sin(alpha))+(khat[k]*dot(khat,oldspin)*(1-cos(alpha)));
}

//Update the spinset (only happens when Monte Carlo condition is met)
void SpinConfig::update_spinset()
{
	for(int k=0;k<3;k++)
		spinset[randsite][k]=newspin[k];
}

//Print final spin configuration to file
void SpinConfig::print_spinset_to_file()
{
	//Output file for final spin configuration
	spinfile=fopen("Results/9x9/spinset-238-run.dat","w");

	//Format the output to be read into Mathematica
	fprintf(spinfile,"{");

	for(int n=0;n<NSITES-1;n++)
		fprintf(spinfile,"{%f,%f,%f},",spinset[n][0],spinset[n][1],spinset[n][2]);

	fprintf(spinfile,"{%f,%f,%f}",spinset[NSITES-1][0],spinset[NSITES-1][1],spinset[NSITES-1][2]);

	fprintf(spinfile,"}");

	fclose(spinfile);
}
