#include "spin-0.1.0.cpp"

class Observables {
	public:
		int count;
		double jij[NSITES][NSITES],oldenergy,deltaenergy,newenergy,magnetization[3],acceptance_rate, average_energy;
		double divide(int x,int y);
		void create_fake_jij();
		void load_real_jij();
		void initial_energy_calculation(double (*)[3]);
		void initial_magnetization_calculation(double (*)[3]);
		void staggered_magnetization_calc(double (*)[3]);
		void print_deltaenergy(FILE *);
		void update_energy_calc(int,int,double (*)[3],double (*),double (*));
		void update_magnetization_calc(int, double(*),double(*));
		void print_energy_to_file(FILE *);
		void print_energy_to_file_gnuplot(FILE *);
		void increment_count();
		void print_acceptance_to_file(FILE *,int);
		void print_acceptance_to_file_gnuplot(FILE *,int);
		void print_magnetization_to_file_gnuplot(FILE *,int);
		void print_average_energy_to_file(FILE *,FILE *,int);
}
;

//create fake j_ij (energy) matrix for testing the Monte Carlo
void Observables::create_fake_jij()
{
	for(int i=0;i<NSITES;i++)
		for(int j=0;j<NSITES;j++)
			jij[i][j] = -.05;
}

//Load j_ij (energy) matrix from file
void Observables::load_real_jij()
{
	FILE * in_file;

	//Variable used to load j_ij values with fscanf
	double d;

	in_file = fopen("9x9/fill-250.dat","r");

	fscanf(in_file,"{");

	//The following for loops take j_ij matrix data and load it into memory (2D array)
	for(int i=0;i<NSITES-1;i++)
	{
		fscanf(in_file,"{");
		for(int j=0;j<NSITES-1;j++)
		{
			fscanf(in_file,"%lf, ",&d);
			jij[i][j]=d;
		}
		fscanf(in_file,"%lf }, ",&d);
		jij[i][NSITES-1]=d;
	}

	fscanf(in_file,"{");
	for(int j=0;j<NSITES-1;j++)
	{
		fscanf(in_file,"%lf,",&d);
		jij[NSITES-1][j]=d;
	}
	fscanf(in_file,"%lf",&d);
	jij[NSITES-1][NSITES-1]=d;

	fclose(in_file);
}

//Cast integer division to a double
double Observables::divide(int x, int y)
{
	double quotient;
	quotient = static_cast<double>(x)/y;
	return quotient;
}

//first energy calculation
void Observables::initial_energy_calculation(double spinset[NSITES][3])
{
	oldenergy = 0;
	for(int i = 0;i<NSITES;i++)
		for(int j= 0;j<NSITES;j++)
			for(int k = 0 ; k<3; k++)
				oldenergy += .5*jij[i][j]*spinset[i][k]*spinset[j][k];
}

//calculation for magnetization
void Observables::initial_magnetization_calculation(double spinset[NSITES][3])
{
	for(int k = 0; k< 3 ;k++)
		for(int i = 0; i < NSITES ;i++)
			magnetization[k]+=spinset[i][k];
}

//calculation for staggered magnetization (fix me!)
void Observables::staggered_magnetization_calc(double spinset[NSITES][3])
{
	for(int k = 0; k< 3 ;k++)
		for(int i = 0; i < NSITES ;i++)
			magnetization[k]+=spinset[i][k];
}

//update energy after a spin is flipped
void Observables::update_energy_calc(int randsite, int nitr,double spinset[NSITES][3],double oldspin[3], double newspin[3])
{
	deltaenergy=0;
	for(int j = 0; j < randsite; j++)
		for(int k=0; k<3; k++)
			deltaenergy += jij[randsite][j]*spinset[j][k]*(newspin[k]-oldspin[k]);

	for(int j=randsite+1;j<NSITES;j++)
		for(int k=0;k<3;k++)
			deltaenergy += jij[randsite][j]*spinset[j][k]*(newspin[k]-oldspin[k]);

	newenergy = oldenergy + deltaenergy;
}

//Print change in energy to file
void Observables::print_deltaenergy(FILE * deltaenergy_out_file)
{
	fprintf(deltaenergy_out_file,"%f \n",deltaenergy);
}

//update magnetization value
void Observables::update_magnetization_calc(int nitr, double oldspin[3],double newspin[3])
{
	for(int k=0;k<3;k++)
	{
		magnetization[k] = magnetization[k]+newspin[k]-oldspin[k];
	}
}

void Observables::print_energy_to_file(FILE * energy_out_file) {

	fprintf(energy_out_file, "%f,",oldenergy);

}

void Observables::print_energy_to_file_gnuplot(FILE * energy_out_file) {

	fprintf(energy_out_file, "%f \n",oldenergy);

}

void Observables::increment_count() {
	count++;
}

void Observables::print_acceptance_to_file(FILE * acceptance_out_file,int nitr)
{
	acceptance_rate = divide(count,nitr);
	fprintf(acceptance_out_file, "%f,",acceptance_rate);
}

void Observables::print_acceptance_to_file_gnuplot(FILE * acceptance_out_file,int nitr)
{
	acceptance_rate = divide(count,nitr);
	fprintf(acceptance_out_file, "%f \n",acceptance_rate);
}

void Observables::print_magnetization_to_file_gnuplot(FILE * magnetization_out_file,int nitr)
{
	fprintf(magnetization_out_file, "%f   %f   %f",magnetization[0],magnetization[1],magnetization[2]);
}

void Observables::print_average_energy_to_file(FILE * energy_in_file, FILE * average_energy_out_file,int nitr)
{
	fscanf(energy_in_file,"{");
	double total_energy=0,d;
	for(int i=0; i<nitr; i++)
	{
		fscanf(energy_in_file,"%lf,",&d);
		total_energy += d;
	}
}
