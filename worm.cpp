/////////////////////////////
// Author: Yu Fu
// Created: 2023-12-12
// Description: C++ program for simulating the hardcore boson dynamics on 3+1D lattice via worm algorithm
/////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <numeric>

struct Parameters {
    int Ns; //spatial site number
    int Nt; //temporal site number
    double Epsilon; //time step length
    double Mu; //chemical potential

    int NsCubeNt; //total site number in "4D space-time volume"
    int NsCube;  // total site number in "3D spatial volume"
    int NsSquare; // total site number in "2D spatial volume"

    int EqCount;  //times of worm-update to "approach to equilibrium" 
    int MeasCount; //times of measuring physicsal observables
    int SweepCount; //times of worm-update bwtween two measurements

    double PF; //time forward probability
    double PB; //time backward probability
    double Psh; // spatial hopping probability

    Parameters(int ns, int nt, double epsilon, double mu, double meas, double sweep)
            : Ns(ns), Nt(nt), Epsilon(epsilon), Mu(mu),
              NsCubeNt(Ns * Ns * Ns * Nt), NsCube(Ns * Ns * Ns), NsSquare(Ns * Ns),
              EqCount(5000), MeasCount(meas), SweepCount(sweep),
              PF(std::exp(mu * epsilon)), PB(std::exp(-1.0 * mu * epsilon)), Psh(6.0 * epsilon) {}
};


//information from lattice configuration after each worm-update
struct Configinfo{
	double Number;
	double Energy;
    double S_hopNum;
	double T_hopNum;
	double Chi_o;
	double Chi_omega;
	double Wx;
	double Wy;
	double Wz;
	double Wsquare;
	double arrow_fwdsite;
};

//coding lattice index
struct Lattice {
    std::vector<std::array<int, 7>> fwdNeighbor;
    Parameters para_inLat;

    Lattice(const Parameters& p) : para_inLat(p), fwdNeighbor(p.NsCubeNt) {}

    void lattice_index() {
        int idx_1D;

        for (int t = 0; t < para_inLat.Nt; t++) {
            for (int z = 0; z < para_inLat.Ns; z++) {
                for (int y = 0; y < para_inLat.Ns; y++) {
                    for (int x = 0; x < para_inLat.Ns; x++) {
                        idx_1D = ((((t * para_inLat.Ns + z) * para_inLat.Ns + y) * para_inLat.Ns) + x);
                        for (int i = 0; i < 7; i++) {
                            fwdNeighbor[idx_1D][i] =
                            ((((t + 1) % para_inLat.Nt) * para_inLat.Ns +
                              ((z + (i == 1) - (i == 4) + para_inLat.Ns) % para_inLat.Ns)) * para_inLat.Ns +
                             ((y + (i == 2) - (i == 5) + para_inLat.Ns) % para_inLat.Ns)) * para_inLat.Ns +
                            ((x + (i == 3) - (i == 6) + para_inLat.Ns) % para_inLat.Ns);
                        }
                    }
                }
            }
        }
    }
};

class MCsimulation {
public:

    Parameters params;
    Lattice lattice;
    Configinfo config;

    std::vector<int> SiteState;
    std::vector<std::array<int, 2>> Linkto;

    std::vector<double> Num_list;
    std::vector<double> Energy_list;
    std::vector<double> Chi_o_list;
    std::vector<double> Chi_omega_list;

    MCsimulation(const Parameters& p);

    void forward(std::mt19937& rng, int& worm, int& moveto, int& pick, int& WormSector, int& arrow);
    void backward(std::mt19937& rng,int& worm, int& moveto, int& pick, int& WormSector, int& arrow);
    void ToEq(std::mt19937& rng);
    void run(std::mt19937& rng);

    Configinfo WormUpdate(std::mt19937& rng);
};

MCsimulation::MCsimulation(const Parameters& p) : params(p), lattice(p) {
    SiteState.resize(params.NsCubeNt, 0);
    Linkto.resize(params.NsCubeNt, {{-1, -1}});
    Num_list.resize(params.MeasCount, 0);
    Energy_list.resize(params.MeasCount, 0);
    Chi_o_list.resize(params.MeasCount, 0);
    Chi_omega_list.resize(params.MeasCount, 0);
}

void MCsimulation::ToEq(std::mt19937& rng) {
    for (int i = 0; i < params.EqCount; i++) {
        WormUpdate(rng);
    }
}


void MCsimulation::backward(std::mt19937& rng, int& worm, int& moveto, int& pick, int& WormSector, int& arrow) {

	if (Linkto[worm][1] != -1) {
		moveto = Linkto[worm][1];
		Linkto[worm][1] = -1;
                Linkto[moveto][0] = -1;
		if (Linkto[worm][0] == -1) {
			SiteState[worm] = 0;
		}
		worm = moveto;
	}
	if(worm == pick){
		if(std::uniform_real_distribution<double>(0, 1)(rng) < params.PB){
			WormSector = 0;
	                SiteState[worm] = 0;
		}else{
			arrow = 1;
		}
	}
}

void MCsimulation::forward(std::mt19937& rng, int& worm, int& moveto, int& pick, int& WormSector, int& arrow) {

	if (std::uniform_real_distribution<double>(0, 1)(rng) < params.Psh) {
		int dice = std::uniform_int_distribution<int>(1, 6)(rng);
		moveto = lattice.fwdNeighbor[worm][dice];
	} else {
		moveto = lattice.fwdNeighbor[worm][0];
	}
	Linkto[worm][0] = moveto;
	if (Linkto[moveto][1] == -1) {
		Linkto[moveto][1] = worm;
		worm = moveto;
		SiteState[worm] = 1;
		if(worm == pick){
			WormSector = 0;
		}
	} else {
		int flag = Linkto[moveto][1];
		Linkto[moveto][1] = worm;
                Linkto[flag][0] = -1;
                worm = flag;
		arrow = 0;
		config.arrow_fwdsite ++;

		if(worm == pick){
			if(std::uniform_real_distribution<double>(0, 1)(rng) < params.PB ){
				WormSector = 0;
				SiteState[worm] = 0;
			}else{
				arrow = 1;
			}
		}
	}
}

Configinfo MCsimulation::WormUpdate(std::mt19937& rng) {

    int pick = std::uniform_int_distribution<int>(0, params.NsCubeNt - 1)(rng);
    int worm = pick;
    int moveto;
    int WormSector = 1;
    int arrow = (SiteState[pick]) ? 0 : 1;
    config.arrow_fwdsite = 0;

    if(arrow == 0){
	    moveto = Linkto[worm][1];
            Linkto[worm][1] = -1;
            Linkto[moveto][0] = -1;
	    if (Linkto[worm][0] == -1) {
		    SiteState[worm] = 0;
	    }
	    worm = moveto;
    }
    else{
	    SiteState[pick] = 1;
	    if (std::uniform_real_distribution<double>(0, 1)(rng) < params.Psh) {
		    int dice = std::uniform_int_distribution<int>(1, 6)(rng);
		    moveto = lattice.fwdNeighbor[worm][dice];
	    } else {
		    moveto = lattice.fwdNeighbor[worm][0];
	    }
	    Linkto[worm][0] = moveto;
	    if (Linkto[moveto][1] == -1) {
		    Linkto[moveto][1] = worm;
		    worm = moveto;
		    SiteState[worm] = 1;
	    } else{
		    int flag = Linkto[moveto][1];
                    Linkto[moveto][1] = worm;
                    Linkto[flag][0] = -1;
                    worm = flag;
		    arrow = 0;
		    config.arrow_fwdsite ++;
	    }
    }

    config.arrow_fwdsite ++;

    WormSector = (worm == pick) ? 0 : 1;
    while (WormSector) {

	    if (arrow == 0){
		    if(std::uniform_real_distribution<double>(0, 1)(rng) < params.PB ){
			    arrow = 0;
		    }
		    else{
			    arrow =1;
			    config.arrow_fwdsite ++;
		    }
	    }
	    else{
		    arrow = 1;
		    config.arrow_fwdsite ++;
	    }
    
	    if (arrow == 0) {
		    backward(rng, worm, moveto, pick, WormSector, arrow);
		    config.arrow_fwdsite ++;
	    }else{
		    forward(rng, worm, moveto, pick, WormSector, arrow);
	    }
    }
//    std::cout << "fwdarrow: " << arrow_fwdsite << std::endl;

     config.S_hopNum = 0.0;
     config.T_hopNum= 0.0;
     config.Energy = 0.0;
     config.Number = 0.0;
     config.Chi_o =0.0;
     config.Chi_omega =0.0;
     config.Wx = 0.0;
     config.Wy = 0.0;
     config.Wz = 0.0;
     config.Wsquare = 0.0;

     for(int i=0; i < params.NsCube; i++){
		config.Number = config.Number + SiteState[i];
	}

     for ( int i=0; i < params.NsCubeNt; i++){
		if (Linkto[i][0] != -1){
			if (Linkto[i][0] == lattice.fwdNeighbor[i][0]){
				config.T_hopNum  +=1;
			}
			else{
				config.S_hopNum  +=1;
			}
		}
	}

for (int i = 0; i < params.NsCubeNt; i++) {
    if (SiteState[i]) {
        int hop_flag = -1;
        for (int j = 0; j < 7; j++) {
            if (lattice.fwdNeighbor[i][j] == Linkto[i][0]) {
                hop_flag = j;
                break; 
            }
        }
        switch (hop_flag) {
            case 1: config.Wz++; break;
            case 2: config.Wy++; break;
            case 3: config.Wx++; break;
            case 4: config.Wz--; break;
            case 5: config.Wy--; break;
            case 6: config.Wx--; break;
            default: break; 
        }
    }
}
//      std::cout << "Temp hop num" << config.T_hopNum<< " " <<  "Spat hop num" << config.S_hopNum <<  std::endl;

        config.Energy = - 1.0/static_cast<double>(params.Nt) * 1.0 /params.Epsilon * config.S_hopNum + 1.0/static_cast<double>(params.Nt) * (6/(1-6*params.Epsilon)) * config.T_hopNum;
	config.Wsquare = (config.Wx/8.0 *config.Wx/8.0);
	config.Chi_o = static_cast<double>(config.arrow_fwdsite) * params.Epsilon;

        return config;
}

void MCsimulation::run(std::mt19937& rng) {
	lattice.lattice_index();
        ToEq(rng);

        for (int i = 0; i < params.MeasCount; i++) {
		double num_avg_per_meas =0.0;
		double energy_avg_per_meas =0.0;
		double chi_o_avg_per_meas =0.0 ;
		double chi_omega_avg_per_meas = 0.0;

		for (int j = 0; j < params.SweepCount; j++) {
			WormUpdate(rng);
			//Measurement after each worm-update
			num_avg_per_meas += config.Number;
			energy_avg_per_meas += config.Energy;
			chi_omega_avg_per_meas += config.Wsquare/static_cast<double>(params.Ns) /params.Epsilon / static_cast<double>(params.Nt);
			chi_o_avg_per_meas += config.Chi_o;
		}
		//Avweage over all updatea in one sweep
		num_avg_per_meas = num_avg_per_meas/static_cast<double>(params.SweepCount);
		energy_avg_per_meas = energy_avg_per_meas/static_cast<double>(params.SweepCount);
		chi_omega_avg_per_meas = chi_omega_avg_per_meas/static_cast<double>(params.SweepCount);
		chi_o_avg_per_meas = chi_o_avg_per_meas/static_cast<double>(params.SweepCount);

		Num_list[i] = num_avg_per_meas;
		Energy_list[i] = energy_avg_per_meas;
		Chi_omega_list[i] =chi_omega_avg_per_meas;
		Chi_o_list[i] =chi_o_avg_per_meas;
	}

	//Average over all measurements
	double num_avg_over_meas = 0.0;
	double energy_avg_over_meas = 0.0;
	double chi_o_avg_over_meas = 0.0;
	double chi_omega_avg_over_meas = 0.0;
	for(int i=0; i < params.MeasCount; i++){

		num_avg_over_meas += Num_list[i];
		energy_avg_over_meas += Energy_list[i];
		chi_omega_avg_over_meas += Chi_omega_list[i];
		chi_o_avg_over_meas += Chi_o_list[i];
	}
	num_avg_over_meas = num_avg_over_meas / static_cast<double>(params.MeasCount);
	energy_avg_over_meas = energy_avg_over_meas / static_cast<double>(params.MeasCount);
	chi_omega_avg_over_meas = chi_omega_avg_over_meas / static_cast<double>(params.MeasCount);
	chi_o_avg_over_meas = chi_o_avg_over_meas / static_cast<double>(params.MeasCount);

        //calculate errors
	double num_error = 0.0;
	double energy_error = 0.0;
	double chi_omega_error = 0.0;
	double chi_o_error =0.0;

	for(int i=0; i < params.MeasCount; i++){
		num_error += (Num_list[i]-num_avg_over_meas) * (Num_list[i]-num_avg_over_meas);
		energy_error += (Energy_list[i]- energy_avg_over_meas) * (Energy_list[i]-energy_avg_over_meas);
                chi_omega_error += (Chi_omega_list[i]- chi_omega_avg_over_meas) * (Chi_omega_list[i]- chi_omega_avg_over_meas);
		chi_o_error += (Chi_o_list[i]- chi_o_avg_over_meas) * (Chi_o_list[i]- chi_o_avg_over_meas);
	}

	num_error = std::pow(num_error,0.5) / static_cast<double>(params.MeasCount);
	energy_error = std::pow(energy_error,0.5) / static_cast<double>(params.MeasCount);
	chi_omega_error = std::pow(chi_omega_error,0.5) / static_cast<double>(params.MeasCount);
	chi_o_error = std::pow(chi_o_error,0.5) / static_cast<double>(params.MeasCount);

	// output to file
	std::ofstream f("result.dat");
	f << "num_avg_over_all_meas: " <<  num_avg_over_meas << " ± " <<  num_error << std::endl;
	f << "energy_avg_over_all_meas: " << energy_avg_over_meas << " ± " << energy_error << std::endl;
	f << "chi_omega_avg_over_all_meas: " << chi_omega_avg_over_meas << " ± " << chi_omega_error << std::endl;
	f << "chi_o_avg_over_all_meas: " << chi_o_avg_over_meas << " ± " << chi_o_error << std::endl;

	f << "----------------------------------------------------" <<std::endl;

	for( int i=0; i < params.MeasCount; i++){
		f << Num_list[i] << " " << Energy_list[i] << " "  << Chi_omega_list[i] << " " << Chi_o_list[i] << std::endl;
	}
}


int main() {
    //Parameters params(32, 110, 0.01, 1.0, 1, 5000); //ns,nt,epsilon,mu, meas, sweep
	Parameters params(8, 110, 0.01, 1.0, 20, 5000); //ns,nt,epsilon,mu, meas, sweep
	//Parameters params(2, 400, 0.03, 1.4, 5, 5000); //ns,nt,epsilon,mu, meas, sweep
    Lattice lat(params);
	MCsimulation simulation(params);
	std::mt19937 rng(std::random_device{}());  // Initialize RNG
	simulation.run(rng);
	return 0;
}
