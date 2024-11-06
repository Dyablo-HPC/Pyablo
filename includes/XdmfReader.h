#pragma once

#include <bits/stdc++.h>

#include "pugixml.hpp"

#include "Snapshot.h"

namespace dyablo
{
  
using TimeSeries = std::vector<std::string>;

class XdmfReader {
public: 
  XdmfReader() = default;
  Snapshot<Geometry_xmf>   readSnapshot(std::string filename);
  Snapshot<Geometry_restart>   readRestart(std::string filename);
  TimeSeries readTimeSeries(std::string filename);
};

}
