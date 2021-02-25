#include <bits/stdc++.h>

#include <omp.h>
#include "Snapshot.h"
#include "XdmfReader.h"

using namespace dyablo;

int main(int argc, char **argv) {
  std::string filename {argv[1]};
  XdmfReader reader;
  Snapshot s = reader.readSnapshot(filename);
  s.print();

  std::vector<float> y;
  Vec v {0.50001, 0.0, 0.0};
  std::vector<Vec> positions;
  float dx = 0.001;
  while (v[1] < 1.0) {
    positions.push_back(v);
    v[1] += dx;
  }

  std::cout << "Probing positions" << std::endl;
  std::vector<float> rho = s.probeDensity(positions);
  std::cout << "DONE !" << std::endl;

  for (auto v: positions)
    std::cout << v[1] << " ";
  std::cout << std::endl;
  for (auto r: rho)
    std::cout << r << " ";
  std::cout << std::endl;

  return 0;
}