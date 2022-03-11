#include <string>
#include <algorithm>

std::string BinaryToString(int bin, int size) {
  std::string str;
  for(int i = 0; i < size; i++) {
    str += (bin & 1) + '0';
    bin = bin >> 1;
  }
  std::reverse(str.begin(), str.end());
  return str;
}
