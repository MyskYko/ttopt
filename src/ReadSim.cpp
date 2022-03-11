#include <cstdio>
#include <cstdlib>

void ReadSim(char *filename, int nInputs, char **pValues, int &nBPatterns) {
  FILE *pFile;
  pFile = fopen (filename, "rb");
  fseek(pFile, 0, SEEK_END);
  int size = ftell(pFile);
  rewind(pFile);
  nBPatterns = size / nInputs;
  printf("%d bytes, %d patterns\n", size, nBPatterns * 8);

  for(int i = 0; i < nInputs; i++) {
    pValues[i] = (char *)malloc(nBPatterns);
    fread(pValues[i], 1, nBPatterns, pFile);
  }

  fclose(pFile);
}
