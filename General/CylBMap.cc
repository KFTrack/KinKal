#include "KinKal/General/CylBMap.hh"

#include <fstream>
#include <sstream>
#include <math.h>

namespace KinKal{
  
  double parse_param_double(const std::string& line, const std::string& substring) {
    // split up grid parameter expected to be a double
    std::size_t i = line.find(substring)+substring.size();
    std::size_t j = line.find(' ', i);
    std::size_t d = j-i;
    return std::stod(line.substr(i, d));
  }

  int parse_param_int(const std::string& line, const std::string& substring) {
    // split up grid parameter expected to be an int
    std::size_t i = line.find(substring)+substring.size();
    std::size_t j = line.find(' ', i);
    std::size_t d = j-i;
    return std::stoi(line.substr(i, d));
  }

  CylBMap::CylBMap(const std::string& file) {
    // read in data from a file
    ntot_ = 0; // track number of lines
    std::string line;
    std::ifstream stream(file);
    while (std::getline(stream, line))
      ++ntot_;
    stream.close();
    // loop through  again, saving data
    // first sift through header
    stream.open(file);
    std::string const data_flag="data";
    std::string const grid_flag="grid";
    double R0=0., Z0=0., dR=0., dZ=0.;
    int nR=0, nZ=0;
    bool grid_parsed=false;
    while (line.substr(0, 4) != data_flag)
    {
      // get the next line
      std::getline(stream, line);
      // parse grid, if parameters in line
      if (line.substr(0, 4) == grid_flag)
      {
        if (grid_parsed) throw std::invalid_argument("multiple grid parameter lines found!");
        R0 = parse_param_double(line, "R0=");
        Z0 = parse_param_double(line, "Z0=");
        nR = parse_param_int(line, "nR=");
        nZ = parse_param_int(line, "nZ=");
        dR = parse_param_double(line, "dR=");
        dZ = parse_param_double(line, "dZ=");
        grid_parsed=true;
      }
      // decrease data count
      --ntot_;
    }
    // throw if no grid parameters found
    if (!grid_parsed) throw std::invalid_argument("no grid parameters found in data file!");
    // construct r and z vectors from grid parameters
    m_ = nR;
    n_ = nZ;
    int ntot_exp = m_*n_;
    if (ntot_ != ntot_exp) throw std::invalid_argument("number of data does not match expected, based on grid parameters");
    r_.resize(m_);
    z_.resize(n_);
    // set r and z from grid parameters
    for (int i=0; i<m_; i++){
      r_[i] = R0 + i * dR;
    }
    for (int j=0; j<n_; j++){
      z_[j] = Z0 + j * dZ;
    }
    // allocate matrices
    // resize outer vectors
    Br_.resize(m_);
    Bz_.resize(m_);
    // resize inner vectors
    for (int i=0; i<m_; i++){
      Br_[i].resize(n_);
      Bz_[i].resize(n_);
    }
    // loop through data lines
    double temp, delta_pos, z_temp, r_temp;
    int step=0, i=0, j=0;
    while(std::getline(stream, line)) {
      i = step / n_;
      j = step % n_;
      std::istringstream iss(line);
      // grab r & z values
      iss >> temp;
      r_temp = temp;
      iss >> temp;
      z_temp = temp;
      // grab field values
      iss >> temp;
      Br_[i][j] = temp;
      iss >> temp;
      Bz_[i][j] = temp;
      // check that r and z grid values are correct
      delta_pos = pow(pow(r_temp-r_[i], 2) + pow(z_temp-z_[j], 2), 0.5);
      if(std::abs(delta_pos) > 1e-5) throw std::invalid_argument("field spacing isn't uniform! dist="+
        std::to_string(delta_pos)+" for line i="+ std::to_string(i)+", j="+std::to_string(j)+
        ", step="+std::to_string(step));
      // increment
      step++;
    }
    // Shrink to fit -- fixes issues with passing in values to InterpBilinear
    r_.shrink_to_fit(), z_.shrink_to_fit(); 
    Br_.shrink_to_fit(), Bz_.shrink_to_fit();
    // shrink inner vectors
    for (int i=0; i<m_; i++){
      Br_[i].shrink_to_fit();
      Bz_[i].shrink_to_fit();
    }
  }
}
