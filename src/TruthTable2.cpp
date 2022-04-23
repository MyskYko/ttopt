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
#include <bitset>
#include <unordered_map>
#include "HashPair.hpp"

using namespace std;
namespace test {
  typedef unsigned long long wrd; // 64 bit type
  vector<wrd> t;     // truth table
  int nWrd;    // number of words in truth table
  int nInput;  // number of inputs
  int nOutput; // number of outputs

  const wrd ones[] = {0x0000000000000001ull,
                      0x0000000000000003ull,
                      0x000000000000000full,
                      0x00000000000000ffull,
                      0x000000000000ffffull,
                      0x00000000ffffffffull,
                      0xffffffffffffffffull};

  wrd GetValue(int idx, int nRowPerSeg) {
    int nSegPerWrd = 64 / nRowPerSeg;
    wrd v = t[idx / nSegPerWrd];
    v >>= (idx % nSegPerWrd) * nRowPerSeg;
    return v & ones[(int)log2(nRowPerSeg)];
  }

  bool IsConst0(int idx, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if(t[nWrdPerSeg * idx + i] != 0)
          return false;
    } else if(GetValue(idx, nRowPerSeg) != 0)
      return false;
    return true;
  }

  bool IsEq(int idx1, int idx2, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if(t[nWrdPerSeg * idx1 + i] !=
           t[nWrdPerSeg * idx2 + i])
          return false;
    } else if(GetValue(idx1, nRowPerSeg) !=
              GetValue(idx2, nRowPerSeg))
      return false;
    return true;
  }

  bool IsConst1(int idx, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if(t[nWrdPerSeg * idx + i] != ones[6])
          return false;
    } else if(GetValue(idx, nRowPerSeg) != ones[(int)log2(nRowPerSeg)])
      return false;
    return true;
  }
  
  bool IsComplEq(int idx1, int idx2, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if(~t[nWrdPerSeg * idx1 + i] !=
           t[nWrdPerSeg * idx2 + i])
          return false;
    } else if((GetValue(idx1, nRowPerSeg) ^ ones[(int)log2(nRowPerSeg)]) !=
              GetValue(idx2, nRowPerSeg))
      return false;
    return true;
  }

  
  // vector of useful indices for each level
  vector<vector<int>> vvIdx;
  // vector of redundant indices for each level
  vector<vector<int>> vvRedIdx;
  
  // hash table type
  typedef unordered_map<pair<int, int>, int> ht;
  // vector of cofactors for each level
  vector<vector<int>> vvCof;

  int FindOrAdd(int idx, int lev);
  /*
  int FindOrAdd(int idx, int lev) {
    if(IsConst0(idx, lev))
      return -2;
    if(IsConst1(idx, lev))
      return -1;
    for(int i = 0; i < vvIdx[lev].size(); i++) {
      if(IsEq(idx, vvIdx[lev][i], lev))
        return i << 1;
      if(IsComplEq(idx, vvIdx[lev][i], lev))
        return (i << 1) | 1;
    }
    vvIdx[lev].push_back(idx);
    return (vvIdx[lev].size() - 1) << 1;
  }
  */

  int CountNodes() {
    vvIdx.assign(nInput, vector<int>());
    vvRedIdx.assign(nInput, vector<int>());
    for(int i = 0; i < nOutput; i++)
      FindOrAdd(i, 0);
    for(int lev = 1; lev < nInput; lev++)
      for(int idx: vvIdx[lev - 1]) {
        int cof0 = FindOrAdd(idx * 2, lev);
        int cof1 = FindOrAdd(idx * 2 + 1, lev);
        if(cof0 == cof1)
          vvRedIdx[lev - 1].push_back(idx);
      }
    int count = 1; // const node
    for(int lev = 0; lev < nInput; lev++)
      count += vvIdx[lev].size() -
        vvRedIdx[lev].size();
    return count;
  }

  /*
  const wrd swapmask[] = {0x2222222222222222ull,
                          0x0c0c0c0c0c0c0c0cull,
                          0x00f000f000f000f0ull,
                          0x0000ff000000ff00ull,
                          0x00000000ffff0000ull};

  void SwapTable(int lev) {
    int nRowPerSeg = 1 << (nInput - lev - 2);
    if(nRowPerSeg >= 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrd; i += nWrdPerSeg * 4)
        for(int j = 0; j < nWrdPerSeg; j++)
          swap(t[i + nWrdPerSeg + j],
               t[i + nWrdPerSeg * 2 + j]);
    } else if(nRowPerSeg == 32)
      for(int i = 0; i < nWrd; i += 2) {
        t[i + 1] ^= t[i] >> 32;
        t[i] ^= t[i + 1] << 32;
        t[i + 1] ^= t[i] >> 32;
      }
    else
      for(int i = 0; i < nWrd; i++) {
        wrd mask = swapmask[(int)log2(nRowPerSeg)];
        t[i] ^= (t[i] >> nRowPerSeg) & mask;
        t[i] ^= (t[i] & mask) << nRowPerSeg;
        t[i] ^= (t[i] >> nRowPerSeg) & mask;
      }
  }
  */
  
  /*
  int FindOrAddUt(int index, int cof0, int cof1,
                  int lev, ht &ut,
                  vector<int> &vCofLow) {
    if(cof0 < 0 && cof0 == cof1)
      return cof0;
    bool fCompl = cof0 & 1;
    if(fCompl)
      cof0 ^= 1, cof1 ^= 1;
    if(ut.count({cof0, cof1}))
      return (ut[{cof0, cof1}] << 1) ^ fCompl;
    vvIdx[lev].push_back(index);
    int id = vvIdx[lev].size() - 1;
    ut[{cof0, cof1}] = id;
    vCofLow.push_back(cof0);
    vCofLow.push_back(cof1);
    if(cof0 == cof1)
      vvRedIdx[lev].push_back(index);
    return (id << 1) ^ fCompl;
  }

  void SwapBDD(int lev) {
    ht ut(2 * vvIdx[lev + 1].size()); // unique table
    vector<int> vCofHigh, vCofLow;
    vvIdx[lev+1].clear();
    vvRedIdx[lev].clear();
    vvRedIdx[lev+1].clear();
    for(int i = 0; i < vvIdx[lev].size(); i++) {
      int idx = vvIdx[lev][i];
      int cof0id = vvCof[lev][i * 2] >> 1;
      bool cof0c = vvCof[lev][i * 2] & 1;
      int cof1id = vvCof[lev][i * 2 + 1] >> 1;
      bool cof1c = vvCof[lev][i * 2 + 1] & 1;
      int cof00, cof01, cof10, cof11;
      if(cof0id < 0)
        cof00 = cof01 = -2 ^ cof0c;
      else
        cof00 = vvCof[lev + 1][cof0id * 2] ^ cof0c,
                cof01 = vvCof[lev + 1][cof0id * 2 + 1] ^
          cof0c;
      if(cof1id < 0)
        cof10 = cof11 = -2 ^ cof1c;
      else
        cof10 = vvCof[lev + 1][cof1id * 2] ^ cof1c,
                cof11 = vvCof[lev + 1][cof1id * 2 + 1] ^
          cof1c;
      int cof0 = FindOrAddUt(idx * 2, cof00, cof10,
                             lev + 1, ut,
                             vCofLow);
      int cof1 = FindOrAddUt(idx * 2 + 1, cof01,
                             cof11, lev + 1, ut,
                             vCofLow);
      vCofHigh.push_back(cof0);
      vCofHigh.push_back(cof1);
      if(cof0 == cof1)
        vvRedIdx[lev].push_back(idx);
    }
    vvCof[lev] = vCofHigh, vvCof[lev + 1] = vCofLow;
  }


  int CountNodes() {
    vvIdx.assign(nInput, vector<int>());
    vvRedIdx.assign(nInput, vector<int>());
    vvCof.assign(nInput, vector<int>());
    for(int i = 0; i < nOutput; i++)
      FindOrAdd(i, 0);
    for(int lev = 1; lev < nInput + 1; lev++)
      for(int idx: vvIdx[lev - 1]) {
        int cof0 = FindOrAdd(idx * 2, lev);
        int cof1 = FindOrAdd(idx * 2 + 1, lev);
        vvCof[lev - 1].push_back(cof0);
        vvCof[lev - 1].push_back(cof1);
        if(cof0 == cof1)
          vvRedIdx[lev - 1].push_back(idx);
      }
    int count = 1; // const node
    for(int lev = 0; lev < nInput; lev++)
      count += vvIdx[lev].size() -
        vvRedIdx[lev].size();
    return count;
  }
  */
  
  vector<wrd> c; // care set

  // update the function and care set at idx2
  void Merge(int idx2, int idx1, int lev, bool fCompl = false);
  wrd GetCare(int index_lev, int lev);
  
  bool TSM(int idx1, int idx2, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if(c[nWrdPerSeg * idx1 + i] &
           c[nWrdPerSeg * idx2 + i] &
           (t[nWrdPerSeg * idx1 + i] ^
            t[nWrdPerSeg * idx2 + i]))
          return false;
    } else if(GetCare(idx1, nRowPerSeg) &
              GetCare(idx2, nRowPerSeg) &
              (GetValue(idx1, nRowPerSeg) ^
               GetValue(idx2, nRowPerSeg)))
      return false;
    Merge(idx2, idx1, lev);
    return true;
  }





  bool ComplTSM(int idx1, int idx2, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if(c[nWrdPerSeg * idx1 + i] &
           c[nWrdPerSeg * idx2 + i] &
           (~t[nWrdPerSeg * idx1 + i] ^
            t[nWrdPerSeg * idx2 + i]))
          return false;
    } else if(GetCare(idx1, nRowPerSeg) &
              GetCare(idx2, nRowPerSeg) &
              (GetValue(idx1, nRowPerSeg) ^
               GetValue(idx2, nRowPerSeg) ^
              ones[(int)log2(nRowPerSeg)]))
      return false;
    Merge(idx2, idx1, lev, true);
    return true;
  }

  bool IsConst0DC(int idx, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if((t[nWrdPerSeg * idx + i] & c[nWrdPerSeg * idx + i]) != 0)
          return false;
    } else if((GetValue(idx, nRowPerSeg) & GetCare(idx, nRowPerSeg)) != 0)
      return false;
    return true;
  }

  bool IsConst1DC(int idx, int lev) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++)
        if((t[nWrdPerSeg * idx + i] | ~c[nWrdPerSeg * idx + i]) != ones[6])
          return false;
    } else if((GetValue(idx, nRowPerSeg) | ~GetCare(idx, nRowPerSeg)) != ones[6])
      return false;
    return true;
  }

  int FindOrAdd(int idx, int lev) {
    if(IsConst0DC(idx, lev))
      return -2;
    if(IsConst1DC(idx, lev))
      return -1;
    for(int i = 0; i < vvIdx[lev].size(); i++) {
      if(TSM(idx, vvIdx[lev][i], lev))
        return i << 1;
      if(ComplTSM(idx, vvIdx[lev][i], lev))
        return (i << 1) | 1;
    }
    vvIdx[lev].push_back(idx);
    return (vvIdx[lev].size() - 1) << 1;
  }

  wrd GetCare(int idx, int nRowPerSeg) {
    int nSegPerWrd = 64 / nRowPerSeg;
    wrd v = c[idx / nSegPerWrd];
    v >>= (idx % nSegPerWrd) * nRowPerSeg;
    return v & ones[(int)log2(nRowPerSeg)];
  }
  
  void SetValue(int index_lev, int lev, wrd value) {
    int logwidth = nInput - lev;
    int index = index_lev >> (6 - logwidth);
    int pos = (index_lev % (1 << (6 - logwidth))) << logwidth;
    t[index] &= ~(ones[logwidth] << pos);
    t[index] ^= value << pos;
  }
  void SetCare(int index_lev, int lev, wrd value) {
    int logwidth = nInput - lev;
    int index = index_lev >> (6 - logwidth);
    int pos = (index_lev % (1 << (6 - logwidth))) << logwidth;
    c[index] &= ~(ones[logwidth] << pos);
    c[index] ^= value << pos;
  }
  
  void Merge(int idx1, int idx2, int lev, bool fCompl) {
    int nRowPerSeg = 1 << (nInput - lev);
    if(nRowPerSeg > 64) {
      int nWrdPerSeg = nRowPerSeg / 64;
      for(int i = 0; i < nWrdPerSeg; i++) {
        if(fCompl)
          t[nWrdPerSeg * idx1 + i] =
            (t[nWrdPerSeg * idx1 + i] &
             ~c[nWrdPerSeg * idx2 + i]) |
            (~t[nWrdPerSeg * idx2 + i] &
             c[nWrdPerSeg * idx2 + i]);
        else
          t[nWrdPerSeg * idx1 + i] =
            (t[nWrdPerSeg * idx1 + i] &
             ~c[nWrdPerSeg * idx2 + i]) |
            (t[nWrdPerSeg * idx2 + i] &
             c[nWrdPerSeg * idx2 + i]);
        c[nWrdPerSeg * idx1 + i] |=
          c[nWrdPerSeg * idx2 + i];
      }
    } else {
      wrd v1 = GetValue(idx1, nRowPerSeg);
      wrd v2 = GetValue(idx2, nRowPerSeg);
      if(fCompl)
        v2 ^= ones[(int)log2(nRowPerSeg)];
      wrd c1 = GetCare(idx1, nRowPerSeg);
      wrd c2 = GetCare(idx2, nRowPerSeg);
      v1 &= c2 ^ ones[(int)log2(nRowPerSeg)];
      v1 |= c2 & v2;
      SetValue(idx1, lev, v1);
      SetCare(idx1, lev, c1 | c2);
    }
  }


  vector<int> vLevels;
/*  
  void SwapIdx(int lev) {
    for(int lev2 = lev + 2; lev2 < nInput; lev2++) {
      int nSubseg = 1 << (lev2 - (lev + 2));
      for(int i = 0; i < vvIdx[lev2].size(); i++) {
        int idx = vvIdx[lev2][i] / nSubseg;
        if(idx % 4 == 1)
          vvIdx[lev2][i] += nSubseg;
        else if(idx % 4 == 2)
          vvIdx[lev2][i] -= nSubseg;
      }
    }
  }

  int RecountNodes(int lev) {
    vvIdx[lev].clear();
    vvIdx[lev+1].clear();
    for(int i = lev; i < lev + 2; i++) {
      if(!i) {
        for(int j = 0; j < nOutput; j++) {
          FindOrAdd(j, 0);
        }
      } else {
        vvRedIdx[i-1].clear();
        for(int idx: vvIdx[i - 1]) {
          int cof0 = FindOrAdd(idx * 2, i);
          int cof1 = FindOrAdd(idx * 2 + 1, i);
          if(cof0 == cof1)
            vvRedIdx[i - 1].push_back(idx);
        }
      }
    }
    if(lev < nInput - 2) {
      vvRedIdx[lev+1].clear();
      for(int index: vvIdx[lev+1]) {
        if(IsEq(index << 1, (index << 1) ^ 1, lev + 2)) {
          vvRedIdx[lev+1].push_back(index);
        }
      }
    }
    int count = 1; // const node
    for(int lev = 0; lev < nInput; lev++)
      count += vvIdx[lev].size() -
        vvRedIdx[lev].size();
    return count;
  }
  
  int Swap(int lev) {
    auto it0 = std::find(vLevels.begin(), vLevels.end(), lev);
    auto it1 = std::find(vLevels.begin(), vLevels.end(), lev + 1);
    std::swap(*it0, *it1);
    // SwapTable(lev);
    // int a = RecountNodes(lev);
    // SwapIdx(lev);
    // return a;
    SwapBDD(lev);
    int count = 1; // const node
    for(int lev = 0; lev < nInput; lev++)
      count += vvIdx[lev].size() -
        vvRedIdx[lev].size();
    return count;
  }

  vector<vector<wrd> > savedt;
  vector<vector<vector<int> > > vvIndicesSaved;
  vector<vector<vector<int> > > vvRedundantIndicesSaved;
  vector<vector<int> > vLevelsSaved;
  vector<vector<vector<int> > > vvChildrenSaved;
  
  // void Save(uint i) {
  //   if(savedt.size() < i + 1) {
  //     savedt.resize(i + 1);
  //     vLevelsSaved.resize(i + 1);
  //   }
  //   savedt[i] = t;
  //   vLevelsSaved[i] = vLevels;
  // }

  // void Load(uint i) {
  //   assert(i < savedt.size());
  //   t = savedt[i];
  //   vLevels = vLevelsSaved[i];
  // }

  // void SaveIndices(uint i) {
  //   if(vvIndicesSaved.size() < i + 1) {
  //     vvIndicesSaved.resize(i + 1);
  //     vvRedundantIndicesSaved.resize(i + 1);
  //   }
  //   vvIndicesSaved[i] = vvIdx;
  //   vvRedundantIndicesSaved[i] = vvRedIdx;
  // }

  // void LoadIndices(uint i) {
  //   vvIdx = vvIndicesSaved[i];
  //   vvRedIdx = vvRedundantIndicesSaved[i];
  // }

  void Save(uint i) {
    if(vLevelsSaved.size() < i + 1) {
      vLevelsSaved.resize(i + 1);
    }
    vLevelsSaved[i] = vLevels;
  }

  void Load(uint i) {
    assert(i < vLevelsSaved.size());
    vLevels = vLevelsSaved[i];
  }

  void SaveIndices(uint i) {
    if(vvIndicesSaved.size() < i + 1) {
      vvIndicesSaved.resize(i + 1);
      vvRedundantIndicesSaved.resize(i + 1);
      vvChildrenSaved.resize(i + 1);
    }
    vvIndicesSaved[i] = vvIdx;
    vvRedundantIndicesSaved[i] = vvRedIdx;
    vvChildrenSaved[i] = vvCof;
  }

  void LoadIndices(uint i) {
    vvIdx = vvIndicesSaved[i];
    vvRedIdx = vvRedundantIndicesSaved[i];
    vvCof = vvChildrenSaved[i];
  }
  
  int BDDNodeCountLevel(int lev) {
    return vvIdx[lev].size() - vvRedIdx[lev].size();
  }
  
  int SiftReo() {
    int best = CountNodes();
    Save(0);
    SaveIndices(0);
    std::vector<int> vars(nInput);
    std::iota(vars.begin(), vars.end(), 0);
    std::sort(vars.begin(), vars.end(), [&](int i1, int i2) {return BDDNodeCountLevel(vLevels[i1]) > BDDNodeCountLevel(vLevels[i2]);});
    bool turn = true;
    for(int var: vars) {
      bool updated = false;
      int lev = vLevels[var];
      for(int i = lev; i < nInput - 1; i++) {
        int count = Swap(i);
        if(best > count) {
          best = count;
          updated = true;
          Save(turn);
          SaveIndices(turn);
        }
      }
      if(lev) {
        Load(!turn);
        LoadIndices(!turn);
        for(int i = lev - 1; i >= 0; i--) {
          int count = Swap(i);
          if(best > count) {
            best = count;
            updated = true;
            Save(turn);
            SaveIndices(turn);
          }
        }
      }
      turn ^= updated;
      Load(!turn);
      LoadIndices(!turn);
    }
    return best;
  }
*/
  
  void test(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
    nInput = inputs.size();
    nOutput = onsets.size();
    int nSize = 0;
    if(nInput >= 6) {
      nSize = 1 << (nInput - 6);
      nWrd = nSize * nOutput;
      t.resize(nWrd);
      for(int i = 0; i < nOutput; i++) {
        for(int pat: onsets[i]) {
          int index = pat / 64;
          int pos = pat % 64;
          t[nSize * i + index] |= 1ull << pos;
        }
      }
    } else {
      nWrd = ((1 << nInput) * nOutput + 64 - 1) / 64;
      t.resize(nWrd);
      for(int i = 0; i < nOutput; i++) {
        int padding = i * (1 << nInput);
        for(int pat: onsets[i]) {
          int pos = (padding + pat) % 64;
          t[padding / 64] |= 1ull << pos;
        }
      }
    }
    vLevels.resize(nInput);
    std::iota(vLevels.begin(), vLevels.end(), 0);

    vector<wrd> care;
    if(nSize) {
      care.resize(nSize);
    } else {
      care.resize(1);
    }
    std::vector<int> count(1 << nInput);
    for(int i = 0; i < nBPats; i++) {
      for(int j = 0; j < 8; j++) {
        int pat = 0;
        for(auto pBPat: pBPats) {
          pat <<= 1;
          pat |= (pBPat[i] >> j) & 1;
        }
        count[pat]++;
        if(count[pat] == rarity) {
          int index = pat / 64;
          int pos = pat % 64;
          care[index] |= 1ull << pos;
        }
      }
    }
    if(nSize) {
      for(int i = 0; i < nOutput; i++) {
        c.insert(c.end(), care.begin(), care.end());
      }
    } else {
      c.resize(nWrd);
      for(int i = 0; i < nOutput; i++) {
        int padding = i * (1 << nInput);
        c[padding / 64] |= care[0] << (padding % 64);
      }
    }

    std::cout << CountNodes() << std::endl;
    //    std::cout << SiftReo() << std::endl;
  }
}
