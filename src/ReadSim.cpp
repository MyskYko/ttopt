#include <cstdio>
#include <cstdlib>
#include <string>
#include <cassert>

void ReadSim(std::string filename, int nInputs, char **pValues, int &nBPatterns) {
  FILE *pFile;
  pFile = fopen (filename.c_str(), "rb");
  fseek(pFile, 0, SEEK_END);
  int size = ftell(pFile);
  rewind(pFile);
  nBPatterns = size / nInputs;
  for(int i = 0; i < nInputs; i++) {
    pValues[i] = (char *)malloc(nBPatterns);
    int r = fread(pValues[i], 1, nBPatterns, pFile);
    assert(r == nBPatterns);
  }
  fclose(pFile);
}
