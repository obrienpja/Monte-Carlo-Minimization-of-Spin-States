#include "observables-0.1.0.cpp"

int main()
{
	//Declare variable that will hold a random number used in the Monte Carlo condition
	double randnum;

	//Declare a spin configuration object as well as an observables object. The observables act on the spin configuration object.
	SpinConfig sc;
	Observables obs;

	//Create initial spinset that will (hopefully) change to the ground state as we iterate
	sc.create_initial_spinset();

	//Load j_ij matrix from file so that we can find the ground state spin configuration for said j_ij matrix
	obs.load_real_jij();

	//initial energy calculation
	obs.initial_energy_calculation(sc.spinset);

	// Create file pointers and open files for writing data
	FILE * energy_out_file;
	FILE * acceptance_out_file;
	FILE * deltaenergy_out_file;

	energy_out_file = fopen("Results/energydata.dat","w");
	acceptance_out_file = fopen("Results/acceptancedata.dat","w");
	deltaenergy_out_file = fopen("Results/deltaenergy.dat","w");

	//monte carlo starts here

	obs.count=0;

	for(int nitr=0;nitr<NITR;nitr++)
	{
		//Track the energy
		obs.print_energy_to_file_gnuplot(energy_out_file);

		//Randomly choose a spin index of a spin we will flip
		randnum = sc.create_rand(0,1);

		//set new and old spin for this iteration
		sc.set_old_spin();
		sc.set_new_spin();

		//as soon as a new spin is chosen, run update energy calculation, then add to energy and mag arrays
		obs.update_energy_calc(sc.randsite,nitr,sc.spinset,sc.oldspin,sc.newspin);
		obs.print_deltaenergy(deltaenergy_out_file);

		obs.update_magnetization_calc(nitr,sc.oldspin,sc.newspin);

		// This is the condition for accepting the spin flip
		if(obs.deltaenergy<0||exp(-obs.deltaenergy/(KB*TEMP)) > randnum)
		{
			sc.update_spinset();
			obs.oldenergy = obs.newenergy;
			obs.increment_count();
		}

		obs.print_acceptance_to_file_gnuplot(acceptance_out_file,(nitr+1));
	}

	//Print results to file
	sc.print_spinset_to_file();

	//Close output files
	fclose(energy_out_file);
	fclose(acceptance_out_file);
	fclose(deltaenergy_out_file);

	return 0;
}
