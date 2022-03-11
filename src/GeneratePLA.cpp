#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>

void GeneratePLA(std::string filename, std::vector<std::vector<int> > onsets, std::vector<char *> pValues, int nBPatterns, int rarity) {
  int LUTSize = pValues.size();
  int nOutputs = onsets.size();
  std::vector<int> careset;
  {
    std::vector<int> count(1 << LUTSize);
    for(int i = 0; i < nBPatterns; i++) {
      for(int j = 0; j < 8; j++) {
        int value = 0;
        for(auto pValue: pValues) {
          value = ((pValue[i] >> j) & 1) + (value << 1);
        }
        count[value]++;
        if(count[value] == rarity) {
          careset.push_back(value);
        }
      }
    }
  }
  
  std::vector<int> count(1 << LUTSize);
  for(int i = 0; i < nOutputs; i++) {
    int pow10 = 1;
    for(int j = i; j < nOutputs-1; j++) {
      pow10 = pow10 * 10;
    }
    for(int value: onsets[i]) {
      count[value] += pow10;
    }
  }

  std::ofstream f(filename);
  f << ".i " << LUTSize << std::endl;
  f << ".o " << nOutputs << std::endl;
  f << ".type fr" << std::endl;
  for(int value: careset) {
    std::string str;
    for(int i = 0; i < LUTSize; i++) {
      str += ((value >> i) & 1) + '0';
    }
    std::reverse(str.begin(), str.end());
    f << str << " " << std::setfill('0') << std::setw(nOutputs) << count[value] << std::endl;
  }

  f.close();
}
