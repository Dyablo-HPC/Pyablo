#pragma once

#include <bits/stdc++.h>

#include "pugixml.hpp"

#include "Snapshot.h"

namespace dyablo
{

class XdmfReader {
public: 
  XdmfReader() = default;
  Snapshot readSnapshot(std::string filename);
};

}
