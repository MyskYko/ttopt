#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <random>
#include <cassert>
#include <map>

extern std::string BinaryToString(int bin, int size);

class TT {
public:
  int nInputs;
  int nSize;
  int nOutputs;
  std::vector<uint> t;

  std::vector<std::vector<int> > vvIndices;
  std::vector<int> vLevels;

  std::mt19937 rng;
  static const uint ones[];
  static const uint sieve[];

  TT(std::vector<std::vector<int> > const &onsets, int nInputs): nInputs(nInputs) {
    assert(nInputs >= 5);
    nSize = 1 << (nInputs - 5);
    nOutputs = onsets.size();
    t.resize(nSize * nOutputs);
    for(int i = 0; i < nOutputs; i++) {
      for(int pat: onsets[i]) {
        int index = pat / 32;
        int pos = pat % 32;
        t[index + nSize * i] |= 1 << pos;
      }
    }
    vLevels.resize(nInputs);
    std::iota(vLevels.begin(), vLevels.end(), 0);
  }

  void GeneratePla(std::string filename) {
    std::ofstream f(filename);
    f << ".i " << nInputs << std::endl;
    f << ".o " << nOutputs << std::endl;
    for(int index = 0; index < nSize; index++) {
      for(int pos = 0; pos < 32; pos++) {
        int pat = (index << 5) + pos;
        f << BinaryToString(pat, nInputs) << " ";
        for(int i = 0; i < nOutputs; i++) {
          f << (char)('0' + ((t[index + nSize * i] >> pos) & 1));
        }
        f << std::endl;
      }
    }
  }

  uint GetValue(int index_lev, int lev) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= 5);
    int index = index_lev >> (lev + 5 - nInputs);
    int pos = (index_lev % (1 << (lev + 5 - nInputs))) << (nInputs - lev);
    return (t[index] >> pos) & ones[nInputs - lev];
  }

  void SetValue(int index_lev, int lev, uint value) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= 5);
    int index = index_lev >> (lev + 5 - nInputs);
    int pos = (index_lev % (1 << (lev + 5 - nInputs))) << (nInputs - lev);
    t[index] &= ~(ones[nInputs - lev] << pos);
    t[index] ^= value << pos;
  }

  int BDDFind(int index, int lev) {
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      bool fZero = true;
      bool fOne = true;
      uint one = ones[5];
      for(int i = 0; i < nScope; i++) {
        uint value = t[nScope * index + i];
        fZero &= value == 0;
        fOne &= value == one;
        if(!fZero && !fOne) {
          break;
        }
      }
      if(fZero) {
        return -2;
      }
      if(fOne) {
        return -1;
      }
      for(int index2: vvIndices[lev]) {
        bool fEq = true;
        bool fCompl = true;
        for(int i = 0; i < nScope; i++) {
          uint value = t[nScope * index + i];
          uint value2 = t[nScope * index2 + i];
          fEq &= value == value2;
          fCompl &= value == ~value2;
          if(!fEq && !fCompl) {
            break;
          }
        }
        if(fEq) {
          return index2 << 1;
        }
        if(fCompl) {
          return (index2 << 1) ^ 1;
        }
      }
    } else {
      uint one = ones[nInputs - lev];
      uint value = GetValue(index, lev);
      if(value == 0) {
        return -2;
      }
      if(value == one) {
        return -1;
      }
      for(int index2: vvIndices[lev]) {
        uint value2 = GetValue(index2, lev);
        if(value2 == value) {
          return index2 << 1;
        }
        if((value2 ^ one) == value) {
          return (index2 << 1) ^ 1;
        }
      }
    }
    return -3;
  }

  int CountBDDNodesOne(int index, int lev) {
    int r = BDDFind(index, lev);
    if(r >= -2) {
      return r;
    }
    vvIndices[lev].push_back(index);
    return index << 1;
  }
  
  int CountBDDNodes() {
    vvIndices.clear();
    vvIndices.resize(nInputs);
    std::vector<std::vector<int> > vvIndicesRedundant(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      CountBDDNodesOne(i, 0);
    }
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0 = CountBDDNodesOne(index << 1, i);
        int cof1 = CountBDDNodesOne((index << 1) ^ 1, i);
        if(cof0 == cof1) {
          vvIndicesRedundant[i-1].push_back(index);
        }
      }
    }
    int count = 1; // const node
    for(int i = 0; i < nInputs; i++) {
      auto it = vvIndicesRedundant[i].begin();
      std::vector<int> vIndicesNew;
      for(int j: vvIndices[i]) {
        if(it == vvIndicesRedundant[i].end() || j != *it) {
          vIndicesNew.push_back(j);
        } else {
          it++;
        }
      }
      vvIndices[i] = vIndicesNew;
      count += vvIndices[i].size();
    }
    return count;
  }

  bool Imply(int index1, int index2, int lev) {
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      for(int i = 0; i < nScope; i++) {
        if(t[nScope * index1 + i] & ~t[nScope * index2 + i]) {
          return false;
        }
      }
      return true;
    }
    return !(GetValue(index1, lev) & (GetValue(index2, lev) ^ ones[nInputs - lev]));
  }

  int GenerateBDDBlifRec(std::vector<std::vector<int> > &vvNodes, int &nNodes, int index, int lev, std::ofstream &f, std::string const &prefix) {
    int r = BDDFind(index, lev);
    if(r >= 0) {
      auto it = std::lower_bound(vvIndices[lev].begin(), vvIndices[lev].end(), r >> 1);
      int i = it - vvIndices[lev].begin();
      return (vvNodes[lev][i] << 1) ^ (r & 1);
    }
    if(r >= -2) {
      return r + 2;
    }
    int cof0 = GenerateBDDBlifRec(vvNodes, nNodes, index << 1, lev + 1, f, prefix);
    int cof1 = GenerateBDDBlifRec(vvNodes, nNodes, (index << 1) ^ 1, lev + 1, f, prefix);
    if(cof0 == cof1) {
      return cof0;
    }
    int cof0id = cof0 >> 1;
    int cof1id = cof1 >> 1;
    bool cof0c = cof0 & 1;
    bool cof1c = cof1 & 1;
    bool imp01 = Imply(index << 1, (index << 1) ^ 1, lev + 1);
    bool imp10 = Imply((index << 1) ^ 1, index << 1, lev + 1);
    f << ".names " << prefix << "v" << lev << " " << prefix << "n" << cof0id << " " << prefix << "n" << cof1id << " " << prefix << "n" << nNodes << std::endl;
    f << (imp01? "-" : "0") << !cof0c << "- 1" << std::endl;
    f << (imp10? "--" : "1-") << !cof1c << " 1" << std::endl;
    vvIndices[lev].push_back(index);
    vvNodes[lev].push_back(nNodes);
    return (nNodes++) << 1;
  }

  void GenerateBDDBlif(std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
    std::string prefix = outputs.front();
    int nNodes = 1; // const node
    vvIndices.clear();
    vvIndices.resize(nInputs);
    std::vector<std::vector<int> > vvNodes(nInputs);
    std::vector<int> vOutputs;
    f << ".names " << prefix << "n0" << std::endl;
    for(int i = 0; i < nInputs; i++) {
      f << ".names " << inputs[i] << " " << prefix << "v" << vLevels[i] << std::endl;
      f << "1 1" << std::endl;
    }
    for(int i = 0; i < nOutputs; i++) {
      int node = GenerateBDDBlifRec(vvNodes, nNodes, i, 0, f, prefix);
      int id = node >> 1;
      bool c = node & 1;
      f << ".names " << prefix << "n" << id << " " << outputs[i] << std::endl;
      f << !c << " 1" << std::endl;
    }
  }

  void SwapLevel(int lev) {
    auto it0 = std::find(vLevels.begin(), vLevels.end(), lev);
    auto it1 = std::find(vLevels.begin(), vLevels.end(), lev + 1);
    std::swap(*it0, *it1);
    if(nInputs - lev > 6) {
      int nScope = 1 << (nInputs - lev - 5 - 2);
      for(int i = nScope; i < nSize * nOutputs; i += (nScope << 2)) {
        for(int j = 0; j < nScope; j++) {
          std::swap(t[i + j], t[i + nScope + j]);
        }
      }
    } else if(nInputs - lev == 6) {
      for(int i = 0; i < nSize * nOutputs; i += 2) {
        t[i+1] ^= t[i] >> 16;
        t[i] ^= t[i+1] << 16;
        t[i+1] ^= t[i] >> 16;
      }
    } else {
      for(int i = 0; i < nSize * nOutputs; i++) {
        int d = (nInputs - lev - 2);
        int shamt = 1 << d;
        t[i] ^= (t[i] >> shamt) & sieve[d];
        t[i] ^= (t[i] & sieve[d]) << shamt;
        t[i] ^= (t[i] >> shamt) & sieve[d];
      }
    }
  }

  int SiftReo() {
    int best = CountBDDNodes();
    std::list<int> vars(nInputs);
    std::iota(vars.begin(), vars.end(), 0);
    while(!vars.empty()) {
      int maxvar = -1;
      uint maxnodes = 0;
      std::list<int>::iterator maxit;
      for(auto it = vars.begin(); it != vars.end(); it++) {
        if(vvIndices[vLevels[maxvar]].size() > maxnodes) {
          maxnodes = vvIndices[vLevels[maxvar]].size();
          maxvar = *it;
          maxit = it;
        }
      }
      if(maxvar == -1) {
        break;
      }
      vars.erase(maxit);
      auto bestt = t;
      auto vLevelsBest = vLevels;
      auto oldt = t;
      auto vLevelsOld = vLevels;
      for(int i = vLevels[maxvar]; i < nInputs - 1; i++) {
        SwapLevel(i);
        int count = CountBDDNodes();
        if(best > count) {
          best = count;
          bestt = t;
          vLevelsBest = vLevels;
        }
      }
      t = oldt;
      vLevels = vLevelsOld;
      for(int i = vLevels[maxvar]-1; i >= 0; i--) {
        SwapLevel(i);
        int count = CountBDDNodes();
        if(best > count) {
          best = count;
          bestt = t;
          vLevelsBest = vLevels;
        }
      }
      t = bestt;
      vLevels = vLevelsBest;
    }
    return best;
  }

  void Reo(std::vector<int> vLevelsNew) {
    for(int i = 0; i < nInputs; i++) {
      int var = std::find(vLevelsNew.begin(), vLevelsNew.end(), i) - vLevelsNew.begin();
      int lev = vLevels[var];
      if(lev < i) {
        for(int j = lev; j < i; j++) {
          SwapLevel(j);
        }
      } else if(lev > i) {
        for(int j = lev-1; j >= i; j--) {
          SwapLevel(j);
        }
      }
    }
    assert(vLevels == vLevelsNew);
  }

  int RandomSiftReo(int nRound) {
    int best = SiftReo();
    auto bestt = t;
    auto vLevelsBest = vLevels;
    for(int i = 0; i < nRound; i++) {
      std::vector<int> vLevelsNew(nInputs);
      std::iota(vLevelsNew.begin(), vLevelsNew.end(), 0);
      std::shuffle(vLevelsNew.begin(), vLevelsNew.end(), rng);
      Reo(vLevelsNew);
      int r = SiftReo();
      if(best > r) {
        best = r;
        bestt = t;
        vLevelsBest = vLevels;
      }
    }
    t = bestt;
    vLevels = vLevelsBest;
    return best;
  }

  void ShowIndices() {
    for(uint i = 0; i < vvIndices.size(); i++) {
      std::cout << "var " << i << ":" << std::endl;
      for(int index: vvIndices[i]) {
        std::cout << index << ", ";
      }
      std::cout << std::endl;
    }
  }
};

const uint TT::ones[] = {0x00000001,
                         0x00000003,
                         0x0000000f,
                         0x000000ff,
                         0x0000ffff,
                         0xffffffff};

const uint TT::sieve[] = {0x22222222,
                          0x0c0c0c0c,
                          0x00f000f0,
                          0x0000ff00};

class TTDC : public TT{
public:
  std::vector<uint> caret;

  TTDC(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TT(onsets, nInputs) {
    std::vector<uint> care(nSize);
    std::vector<int> count(1 << nInputs);
    for(int i = 0; i < nBPats; i++) {
      for(int j = 0; j < 8; j++) {
        int pat = 0;
        for(auto pBPat: pBPats) {
          pat <<= 1;
          pat |= ((pBPat[i] >> j) & 1);
        }
        count[pat]++;
        if(count[pat] == rarity) {
          int index = pat / 32;
          int pos = pat % 32;
          care[index] |= 1 << pos;
        }
      }
    }
    for(int i = 0; i < nOutputs; i++) {
      caret.insert(caret.end(), care.begin(), care.end());
    }
  }

  uint GetCare(int index_lev, int lev) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= 5);
    int index = (index_lev >> (lev + 5 - nInputs));
    int pos = (index_lev % (1 << (lev + 5 - nInputs))) << (nInputs - lev);
    return (caret[index] >> pos) & ones[nInputs - lev];
  }

  bool CheckDC(int index, int lev) {
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      for(int i = 0; i < nScope; i++) {
        uint value = caret[nScope * index + i];
        if(value != 0) {
          return false;
        }
      }
      return true;
    } else {
      uint value = GetCare(index, lev);
      if(value != 0) {
        return false;
      }
      return true;
    }
  }

  void MergeCare(int index1, int index2, int lev) {
    assert(index2 >= 0);
    if(index1 < 0) {
      return;
    }
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      for(int i = 0; i < nScope; i++) {
        caret[nScope * index1 + i] |= caret[nScope * index2 + i];
      }
    } else {
      uint value = GetCare(index2, lev);
      int index = index1 >> (lev + 5 - nInputs);
      int pos = (index1 % (1 << (lev + 5 - nInputs))) << (nInputs - lev);
      caret[index] |= value << pos;
    }
  }

  int CountBDDNodesOSM() {
    auto oldcaret = caret;
    vvIndices.clear();
    vvIndices.resize(nInputs);
    std::vector<std::vector<int> > children(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(CheckDC(i, 0)) {
        continue;
      }
      int r = CountBDDNodesOne(i, 0);
      if(r != i << 1) {
        MergeCare(r >> 1, i, 0);
      }
    }
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0 = -3, cof1 = -3;
        bool cof0dc = CheckDC(index << 1, i);
        if(!cof0dc) {
          cof0 = CountBDDNodesOne(index << 1, i);
          if(cof0 != index << 2) {
            MergeCare(cof0 >> 1, index << 1, i);
          }
          children[i-1].push_back(cof0);
        }
        bool cof1dc = CheckDC((index << 1) ^ 1, i);
        if(!cof1dc) {
          cof1 = CountBDDNodesOne((index << 1) ^ 1, i);
          if(cof1 != ((index << 2) ^ 2)) {
            MergeCare(cof1 >> 1, (index << 1) ^ 1, i);
          }
          children[i-1].push_back(cof1);
        }
        assert(!cof0dc || !cof1dc);
        if(cof0dc) {
          children[i-1].push_back(cof1);
        }
        if(cof1dc) {
          children[i-1].push_back(cof0);
        }
      }
    }
    std::vector<std::vector<int> > vvIndicesRedundant(nInputs);
    std::map<int, std::pair<int, int> > skipped;
    for(int i = nInputs - 2; i >= 0; i--) {
      std::map<int, std::pair<int, int> > nextskipped;
      std::map<std::pair<std::pair<int, int>, std::pair<int, int> >, int> unique;
      for(uint j = 0; j < vvIndices[i].size(); j++) {
        std::pair<int, int> cof0 = {i+1, children[i][j+j]};
        std::pair<int, int> cof1 = {i+1, children[i][j+j+1]};
        int cof0ind = children[i][j+j] >> 1;
        int cof1ind = children[i][j+j+1] >> 1;
        bool cof0c = children[i][j+j] & 1;
        bool cof1c = children[i][j+j+1] & 1;
        if(cof0ind < 0) {
          cof0.first = nInputs;
        } else if(skipped.count(cof0ind)) {
          cof0 = skipped[cof0ind];
          cof0.second ^= cof0c;
        }
        if(cof1ind < 0) {
          cof1.first = nInputs;
        } else if(skipped.count(cof1ind)) {
          cof1 = skipped[cof1ind];
          cof1.second ^= cof1c;
        }
        if(cof0 == cof1) {
          nextskipped[vvIndices[i][j]] = cof0;
          vvIndicesRedundant[i].push_back(vvIndices[i][j]);
          continue;
        }
        bool neg = cof0.second & 1;
        if(neg) {
          cof0.second ^= 1;
          cof1.second ^= 1;
        }
        if(unique.count({cof0, cof1})) {
          nextskipped[vvIndices[i][j]] = {i, unique[{cof0, cof1}] ^ neg};
          vvIndicesRedundant[i].push_back(vvIndices[i][j]);
          continue;
        }
        unique[{cof0, cof1}] = (vvIndices[i][j] << 1) ^ neg;
      }
      skipped = nextskipped;
    }
    int count = 1; // const node
    for(int i = 0; i < nInputs; i++) {
      auto it = vvIndicesRedundant[i].begin();
      std::vector<int> vIndicesNew;
      for(int j: vvIndices[i]) {
        if(it == vvIndicesRedundant[i].end() || j != *it) {
          vIndicesNew.push_back(j);
        } else {
          it++;
        }
      }
      vvIndices[i] = vIndicesNew;
      count += vvIndices[i].size();
    }
    caret = oldcaret;
    return count;
  }

  void CopyFunc(int index1, int index2, bool neg, int lev) {
    assert(index1 >= 0);
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      uint one = ones[5];
      for(int i = 0; i < nScope; i++) {
        if(index2 < 0) {
          t[nScope * index1 + i] = neg? one: 0;
        } else {
          t[nScope * index1 + i] = neg? ~t[nScope * index2 + i]: t[nScope * index2 + i];
        }
      }
    } else {
      uint value;
      uint one = ones[nInputs - lev];
      if(index2 < 0) {
        value = neg? one: 0;
      } else {
        value = GetValue(index2, lev);
        if(neg) {
          value ^= one;
        }
      }
      SetValue(index1, lev, value);
    }
  }

  void OSM() {
    std::vector<std::vector<std::pair<int, int> > > merged(nInputs);
    vvIndices.clear();
    vvIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(CheckDC(i, 0)) {
        for(int j = 0; j < nSize; j++) {
          t[j + nSize * i] = 0;
        }
        continue;
      }
      int r = CountBDDNodesOne(i, 0);
      if(r != i << 1) {
        MergeCare(r >> 1, i, 0);
        merged[0].push_back({r, i});
      }
    }
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0, cof1;
        bool cof0dc = CheckDC(index << 1, i);
        if(!cof0dc) {
          cof0 = CountBDDNodesOne(index << 1, i);
          if(cof0 != index << 2) {
            MergeCare(cof0 >> 1, index << 1, i);
            merged[i].push_back({cof0, index << 1});
          }
        }
        bool cof1dc = CheckDC((index << 1) ^ 1, i);
        if(!cof1dc) {
          cof1 = CountBDDNodesOne((index << 1) ^ 1, i);
          if(cof1 != ((index << 2) ^ 2)) {
            MergeCare(cof1 >> 1, (index << 1) ^ 1, i);
            merged[i].push_back({cof1, (index << 1) ^ 1});
          }
        }
        assert(!cof0dc || !cof1dc);
        if(cof0dc) {
          merged[i].push_back({(index << 2) ^ 2, index << 1});
        }
        if(cof1dc) {
          merged[i].push_back({index << 2, (index << 1) ^ 1});
        }
      }
    }
    for(int i = nInputs - 1; i >= 0; i--) {
      for(auto p: merged[i]) {
        CopyFunc(p.second, p.first >> 1, p.first & 1, i);
      }
    }
  }
};

void TTTest(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
  int nInputs = inputs.size();
  TT tt(onsets, nInputs);
  std::cout << tt.CountBDDNodes() << std::endl;
  TTDC ttdc(onsets, nInputs, pBPats, nBPats, rarity);
  std::cout << ttdc.CountBDDNodes() << std::endl;
  std::cout << ttdc.CountBDDNodesOSM() << std::endl;
  ttdc.OSM();
  std::cout << ttdc.CountBDDNodes() << std::endl;
  ttdc.GenerateBDDBlif(inputs, outputs, f);
}
