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
};
