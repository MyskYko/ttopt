#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <cassert>

void ReadSim(std::string filename, int nInputs, std::vector<char *> &vpBPats, int &nBPats) {
  for(auto pBPats: vpBPats) {
    free(pBPats);
  }
  vpBPats.clear();
  FILE *pFile;
  pFile = fopen (filename.c_str(), "rb");
  fseek(pFile, 0, SEEK_END);
  int size = ftell(pFile);
  rewind(pFile);
  nBPats = size / nInputs;
  vpBPats.resize(nInputs);
  for(int i = 0; i < nInputs; i++) {
    vpBPats[i] = (char *)malloc(nBPats);
    int r = fread(vpBPats[i], 1, nBPats, pFile);
    assert(r == nBPats);
  }
  fclose(pFile);
}
