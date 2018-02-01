#include "Setup.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include <cmath>

#include "Utils.h"

// using namespace std;

namespace Siggen
{

Setup::Setup(std::string filename)
{
	read_config(filename);
}

int Setup::read_config(std::string config_file_name) {
  /* reads and parses configuration file of name config_file_name
     returns 0 on success, 1 otherwise
  */

  int n=0;
  std::string line, key, value, conf_set, conf_path;
  std::map<std::string,std::string>* currentMap;

  std::ifstream file(config_file_name.c_str());
  if(!file){
    std::cout <<"\nERROR: config file " << config_file_name << " does not exist?\n";
    return 1;
  }

  //Figure out the path of the config file
  size_t found;
  found=config_file_name.find_last_of("/\\");
  conf_path = config_file_name.substr(0,found+1);

  std::cout << "\nReading values from config file " << config_file_name << "\n";

  while(getline(file, line))
  {
    n++;

    // ignore comments and blank lines
    if (line.length() <3 || line.find_first_not_of(' ') == line.npos || line[line.find_first_not_of(' ')] == '#') continue;

    //look for lines that start with "[", which denote a new group
    int first_idx = line.find_first_not_of(" ");
    if (line[first_idx] == '['){
      conf_set = line.substr(first_idx+1, line.find_first_of("]") - first_idx-1);
      if (conf_set.compare("general")==0) continue;
      else if (conf_set.compare("geometry")==0) currentMap = &geometry_map;
      // else if (conf_set.compare("fieldgen")==0) currentMap = &fieldgen_map;
      else if (conf_set.compare("siggen")==0)   currentMap = &siggen_map;
      else if (conf_set.compare("detector")==0) currentMap = &detector_map;
      else{
        std::cout << "\nERROR: group name "<< conf_set << " not valid (line number " << n << ").\n";
        std::cout << "--> full line: " << line << std::endl;
        return 1;
      }
      continue;
    }

    //parse the key, value pair (first two things separated by white space)
    //turn to stringstream for easy parsing
    std::stringstream line_stream(line);
    line_stream >> key >> value;
    if (line_stream.fail()){
      std::cout << "\nERROR: couldn't parse line number " << n <<": " << line << "\n";
      return 1;
    }

    std::stringstream valstream(value);
    if (key == "verbosity_level"){
      valstream >> verbosity;
      verbosity_setting = (verbosity_level) verbosity;
      continue;
    }else if (key == "xtal_radius"){
      valstream >> xtal_radius;
    }else if (key == "xtal_length"){
      valstream >> xtal_length;
    }else if (key == "xtal_grid"){
      valstream >> xtal_grid;
    }else if (key == "field_name"){
      valstream >> field_name;
      field_name = conf_path + field_name;
    }else if (key == "wp_name"){
      valstream >> wp_name;
      wp_name = conf_path + wp_name;
    }else if (key == "imp_min"){
      valstream >> imp_min;
    }else if (key == "imp_max"){
      valstream >> imp_max;
    }else if (key == "imp_num"){
      valstream >> imp_num;
    }else if (key == "grad_min"){
      valstream >> grad_min;
    }else if (key == "grad_max"){
      valstream >> grad_max;
    }else if (key == "grad_num"){
      valstream >> grad_num;
    }

    (*currentMap)[key] = value;

  }

  std::cout << "Done reading file\n";
  return 0;
}

    //



  	// } else if (strstr(key_word[i], "top_bullet_radius")) {
  	//   top_bullet_radius = fi;
  	// } else if (strstr(key_word[i], "bottom_bullet_radius")) {
  	//   bottom_bullet_radius = fi;
  	// } else if (strstr(key_word[i], "pc_length")) {
  	//   pc_length = fi;
  	// } else if (strstr(key_word[i], "pc_radius")) {
  	//   pc_radius = fi;
  	// } else if (strstr(key_word[i], "bulletize_PC")) {
  	//   bulletize_PC = ii;
  	// } else if (strstr(key_word[i], "taper_length")) {
  	//   taper_length = fi;
  	// } else if (strstr(key_word[i], "wrap_around_radius")) {
  	//   wrap_around_radius = fi;
  	// } else if (strstr(key_word[i], "ditch_depth")) {
  	//   ditch_depth = fi;
  	// } else if (strstr(key_word[i], "ditch_thickness")) {
  	//   ditch_thickness = fi;
    // } else if (strstr(key_word[i], "central_hole_l")) {
    //   central_hole_l = fi;
    // } else if (strstr(key_word[i], "central_hole_r")) {
    //   central_hole_r = fi;
    // } else if (strstr(key_word[i], "taper_z")) {
    //   taper_z = fi;
    // } else if (strstr(key_word[i], "taper_dr")) {
    //   taper_dr = fi;



  	// } else {
  	//   printf("WARNING: unrecognized keyword %s\n", key_word[i]);
  	//   return 1;
  	// }
    //
  	// if (verbosity >= CHATTY) {
  	//   // printf("%s", line);
  	//   if (iint) {
  	//     printf("%s: %d\n", key_word[i], ii);
  	//   } else if (strlen(name) > 0) {
  	//     printf("%s: %s\n", key_word[i], name);
  	//   } else {
  	//     printf("%s: %f\n", key_word[i], fi);
  	//   }
  	// }



// int parse_value()

}//namespace Siggen
