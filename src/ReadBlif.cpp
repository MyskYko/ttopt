#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>

void ReadBlifHeader(std::ifstream &f, std::string &modulename, std::vector<std::string> &inputs, std::vector<std::string> &outputs) {
  std::string str;
  while(std::getline(f, str)) {
    if(str.length() > 6 && str.substr(0, 6) == ".model") {
      break;
    }
  }
  modulename = str.substr(7);
  while(std::getline(f, str)) {
    if(str.length() > 7 && str.substr(0, 7) == ".inputs") {
      break;
    }
  }
  {
    std::stringstream ss(str);
    std::string input;
    std::getline(ss, input, ' ');
    while(std::getline(ss, input, ' ')) {
      inputs.push_back(input);
    }
  }
  while(std::getline(f, str)) {
    if(str.length() > 8 && str.substr(0, 8) == ".outputs") {
      break;
    }
  }
  {
    std::stringstream ss(str);
    std::string output;
    std::getline(ss, output, ' ');
    while(std::getline(ss, output, ' ')) {
      outputs.push_back(output);
    }
  }
}

int ReadBlifFuncs(std::ifstream &f, int nGroupSize, std::vector<std::string> &LUTInputs, std::vector<std::string> &LUTOutputs, std::vector<std::vector<int> > &onsets) {
  LUTInputs.clear();
  LUTOutputs.clear();
  onsets.clear();
  std::string str;
  for(int i = 0; i < nGroupSize; i++) {
    while(std::getline(f, str)) {
      if(str.length() > 6 && str.substr(0, 6) == ".names") {
        break;
      }
    }
    if(!f) {
      return 0;
    }
    std::vector<std::string> LUTInputs_;
    {
      std::stringstream ss(str);
      std::string input;
      std::getline(ss, input, ' ');
      while(std::getline(ss, input, ' ')) {
        LUTInputs_.push_back(input);
      }
      std::string LUTOutput = LUTInputs_.back();
      LUTOutputs.push_back(LUTOutput);
      LUTInputs_.pop_back();
    }
    if(i == 0) {
      LUTInputs = LUTInputs_;
    } else {
      assert(LUTInputs == LUTInputs_);
    }
    std::vector<int> onset;
    int pos;
    while(pos = f.tellg(), std::getline(f, str)) {
      if(str.length() == 0 || (str[0] != '0' && str[0] != '1')) {
        break;
      }
      str = str.substr(0, LUTInputs.size());
      int pat = std::stoi(str, 0, 2);
      onset.push_back(pat);
    }
    f.seekg(pos);
    onsets.push_back(onset);
  }
  return 1;
}
