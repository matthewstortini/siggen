#include <iostream>
#include <fstream>

#include "Siggen.h"
// #include "ICPC.h"
#include <vector>

using namespace std;
using namespace Siggen;

int main(int argc, char** argv)
{

  // Setup setup("config_files/icpc.config");
  // Detector<ICPC>* icpc = new Detector<ICPC>(setup);
  // icpc->field_setup();
  // SignalGenerator<ICPC> siggen(icpc, setup);

  Setup setup("config_files/ortec4b.config");
  Detector<GEM>* gem = new Detector<GEM>(setup);
  gem->field_setup();
  SignalGenerator<GEM> siggen(gem, setup);

  int nsteps = siggen.get_output_length();
  int nsegs = siggen.get_nsegments();

  std::cout << "output steps: "<< nsteps << std::endl;

  vector<float> output_arr(nsteps*nsegs,0.);

  point pt;
  pt.x=20;
  pt.y=0;
  pt.z=20;

  int flag = siggen.get_signal(pt, &output_arr[0]);

  ofstream myfile;
  myfile.open ("output.txt");

  for(int i = 0; i< nsteps; i++) {
      for (int j = 0; j< nsegs; j++){
          myfile << output_arr[j*nsteps+i];
          if (j < nsegs-1) myfile <<",";
      }
      myfile << std::endl;
  }
  myfile.close();

	// start<StraightLine>(argc, argv);
	return 0;
}
