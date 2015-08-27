#include <cstdlib>
#include <cstring>
#include <iostream>

#include "BD.h"

using std::atof;
using std::atoi;
using std::cout;
using std::endl;
using std::exit;
using std::strncmp;

/******************************************************************/
/******************************************************************//**
 * Main for BD simulation of Barnase/Barstar
 * \param a floating point number of salt concentration
 * \param fout an output file
 * \param ind an integer for temporary file index
 ******************************************************************/
int main1(int argc, char **argv) {
	if (argc != 7) {
		cout << "Correct input format: " << endl;
		cout << " ./exec sim [protein 1] [protein 2] [Salt conc] [outfile]"
				" [temp file #]" << endl;
		exit(0);
	}
	return run_sim(argv[2], argv[3], atof(argv[4]), argv[5], atoi(argv[6]));
} // end main1

/******************************************************************/
/******************************************************************//**
 * Main for computing energies, forces and torques on a replication
 * of a single molecule places equidistantly in a system.
 *   \param ifname an input file name
 *   \param num an integer number of molecules to include in the system.
 *              Values are: 2, 4, 6 or 8
 *   \param dist a floating point distance for molecule placement
 ******************************************************************/
int main2(int argc, char **argv) {
	if (argc != 5) {
		cout << "Correct input format: " << endl;
		cout << " ./exec slv [PQR file] [2, 4, 6 or 8 molecules]"
				" [distance between molecules]" << endl;
		exit(0);
	}

	cout << "SOLVE" << endl;
	return run_slv(argv[2], atoi(argv[3]), atof(argv[4]));
} // end main2

/******************************************************************/
/******************************************************************//**
 * Main for perturbation run.
 *   \param ifname a character pointer holding input file name
 *   \param num an integer describing the number of iterations to
 *              run force calculations
 *   \param dist an int describing the distance for molecule placement
 *   \param Dtr a floating point number containing the translational
 *              diffusion coefficient of input molecule
 *   \param Dr a floating point number containing the rotational
 *             diffusion coefficient of input molecule
 *   \param a string containing the perturbation file name
 ******************************************************************/
int main3(int argc, char **argv) {
	if (argc != 8) {
		cout << "Correct input format: " << endl;
		cout << " ./exec per [PQR file] [2, 4, 6 or 8 molecules]"
				" [distance between molecules] [translational diff. coeff]"
				" [rotational diff coeff] [runname]" << endl;
		exit(0);
	}
	cout << "PERTURB" << endl;

	return run_per(argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]),
			atof(argv[6]), argv[7]);
} // end main3

/******************************************************************/
/******************************************************************//**
 * Main for Computing the polarization forces, used for comparing the
 * effect of mutual polarization on force and torque.  The output is
 * a file named polar_[force/torque]_name_nmol_dist.txt.  The first line
 * indicates the force or torque computed per molecule in the absence of
 * mutual polarization.  The next line includes mutual polarization.
 * These two lines repeat for 1000 iterations.
 *   \param ifname is a character pointer for input file name
 *   \param num is an int of the number of molecules in the system
 *   \param dist is a floating point number indicating the
 *               distance for molecule placement
 *	\param a character string to describe the system
 ******************************************************************/
int main4(int argc, char **argv) {
	if (argc != 6) {
		cout << "Correct input format: " << endl;
		cout << " ./exec pol [PQR file] [2, 4, 6 or 8 molecules]"
				" [distance between molecules] [runname]" << endl;
		exit(0);
	}
	cout << "POL" << endl;

	return run_pol(argv[2], atoi(argv[3]), atoi(argv[4]), argv[5]);
} // end main4

/******************************************************************/
/******************************************************************//**
 * Main for computing the potential on a grid, given a number of num
 * identical molecules placed dist apart in a salt solution.
 *   \param ifname a character string containing an input file name
 *   \param num an integer number of molecules to introduce into system
 *   \param dist a floating point number for distance for molecule placement
 ******************************************************************/
int main5(int argc, char **argv) {
	if (argc != 5) {
		cout << "Correct input format: " << endl;
		cout << " ./exec dif [PQR file] [2, 4, 6 or 8 molecules]"
				" [distance between molecules]" << endl;
		exit(0);
	}
	cout << "GDIFF" << endl;

	return run_dif(argv[2], atoi(argv[3]), atof(argv[4]));
} // end main5

/******************************************************************/
/******************************************************************//**
 * Main for computing self-polarization and potential grid for a
 * single coarse-grained molecule positioned at (0,0,0).  Similar to
 * the option dif, or main 5, but only for a single molecule and
 *   \param ifname a character string describing the input file name
 *   \param fact a floating point scaling factor for the MPE radius
 ******************************************************************/
int main6(int argc, char **argv) {
	if (argc != 4) {
		cout << "Correct input format: " << endl;
		cout << " ./exec rad [PQR file] [scaling factor]" << endl;
		exit(0);
	}
	cout << "RAD" << endl;

	return run_rad(argv[2], atof(argv[3]));
} // end of main6

/******************************************************************/
/******************************************************************//**
 * Main for computing the effect of a third molecule on the interaction
 * free energy of two molecules, as a function of separation distance.
 * The output file generated, [runname]_dist*10.txt has 900 lines, each line
 * representing the interaction free energy between two molecules.
 * The first value in each line represents the system of just two molecules,
 * the next 30 values represents including a 3rd molecule in the system at 30
 * different orientations.
 *   \param ifname a character string describing input file name
 *   \param dist a floating point of the distance between two molecules
 *   \param a string for output for runname
 ******************************************************************/
int main7(int argc, char **argv) {
	if (argc != 5) {
		cout << "Correct input format: " << endl;
		cout << " ./exec cng [PQR file] [distance between molecles] [runname]"
				<< endl;
		exit(0);
	}
	cout << "CHANGE" << endl;

	return run_cng(argv[2], atof(argv[3]), argv[4]);
} // main7

/******************************************************************/
/******************************************************************//**
 * Main for making an infinte grid.  It prints out the potential, the force
 * and the torque for 8 molecules in a lattice with a given number of lattice
 * layers
 *   \param ifname a character string with input file name
 *   \param layer an integer describing the number of neighbors to consider
 *                when performing energy, force and torque calculations
 *   \param stretch a floating point number describing the scaling of
 *                  the CG spheres
 ******************************************************************/
int main8(int argc, char **argv) {
	if (argc != 5) {
		cout << "Correct input format: " << endl;
		cout << " ./exec inf [PQR file] [# layers] [stretch factor]" << endl;
		exit(0);
	}
	cout << "INFINITE GRID" << endl;

	return run_inf(argv[2], atoi(argv[3]), atof(argv[4]));
} // end main8

/******************************************************************//**
 * Main!
 ******************************************************************/
int main(int argc, char **argv) {
	// For running 2 molecule BD simulation
	if (strncmp(argv[1], "sim", 3) == 0)
		return main1(argc, argv);

	// For computing energies, torques and forces of many of the same
	// molecules in solution
	else if (strncmp(argv[1], "slv", 3) == 0)
		return main2(argc, argv);

	// For computing the energy of many molecules in solution as their
	// rotations and locations are perturbed
	else if (strncmp(argv[1], "per", 3) == 0)
		return main3(argc, argv);

	// For computing the forces/torques of many molecules  in solution
	// with and without mutual polarization
	else if (strncmp(argv[1], "pol", 3) == 0)
		return main4(argc, argv);

	// For computing a grid of potentials due to many molecules fixed
	// in solution
	else if (strncmp(argv[1], "dif", 3) == 0)
		return main5(argc, argv);

	// For computing a potential grid for a single molecule in solution
	else if (strncmp(argv[1], "rad", 3) == 0)
		return main6(argc, argv);

	// For computing the effect of a 3rd molecule of the total free energy
	// of a system at distance dist (in A)
	else if (strncmp(argv[1], "cng", 3) == 0)
		return main7(argc, argv);

	else if (strncmp(argv[1], "inf", 3) == 0)
		return main8(argc, argv);

	// else, bad option
	else {
		cout << "bad option!!! " << argv[1] << endl;
		return 1;
	}
}
