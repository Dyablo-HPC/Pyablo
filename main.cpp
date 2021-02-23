#include <bits/stdc++.h>

#include "Snapshot.h"
#include "XdmfReader.h"

using namespace dyablo;

int main(int argc, char **argv) {
  std::string filename {argv[1]};
  XdmfReader reader;
  Snapshot s = reader.readSnapshot(filename);
  s.print();

  auto [min, max] = s.getCellBoundingBox(0);
  std::cout << "Bounding box of cell #0 : (" << min[0] << "; " << min[1] 
    << ") - (" << max[0] << "; " << max[1] << ")" << std::endl;

  Vec pos{0.51, 0.51};
  int iCell = s.getCellFromPosition(pos);
  std::cout << "Cell holding position : (0.5, 0.5) is Cell #" << iCell << std::endl;
  
  std::vector<float> y;
  std::vector<float> rho;
  Vec v {0.50001, 0.0, 0.0};
  float dx = 0.001;
  while (v[1] < 1.0) {
    float lrho = s.probeDensity(v);
    y.push_back(v[1]);
    rho.push_back(lrho);
    v[1] += dx;
  }

  for (auto vy: y)
    std::cout << vy << " ";
  std::cout << std::endl;
  for (auto r: rho)
    std::cout << r << " ";
  std::cout << std::endl;

  return 0;
}