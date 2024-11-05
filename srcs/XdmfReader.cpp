#include "XdmfReader.h"

namespace dyablo {

// Local helpers
auto splitH5Filename(std::string xpath) {
  std::string file {""};
  std::string path {""};
  int mode = 0;
  for (auto c : xpath) {
    if (c == ' ' || c == '\n')
      continue;
    if (c == ':')
      mode++;
    else {
      if (mode == 0)
        file += c;
      else
        path += c;
    }
  }
  return std::make_pair(file, path);
}

inline int extractDimensions(std::string dims) {
  std::string d1 = dims.substr(0, dims.find(" "));
  return std::atoi(d1.c_str());
}

/**
 * Reads a snapshot from an xdmf file
 * @param filename the path to the file to read
 * @return A Snapshot object
 **/
Snapshot XdmfReader::readSnapshot(std::string restart_filename) {
  Snapshot snap;

  std::map<std::string, hid_t> h5_handles;

  snap.setTime(0.2);
  snap.setName("test");
  snap.setBlockSize(4,4,4);
  snap.setDomainBoundingBox(std::make_pair<Vec, Vec>({0,0,0},{1,1,1}));


  /**
   * Connectivity info
   **/
  snap.addH5Handle("restart_file", restart_filename);
  snap.setOctCoordinates("restart_file", "/Octree");

  /**
   * Attributes info
   **/
  for (auto field_name: {"rho","e_tot","rho_vx"}) {    
    std::string xpath = std::string("/fields/") + field_name;
    snap.addAttribute("restart_file", xpath, field_name, "Double", "?");  
  }
  return snap;
}

/**
 * Reads a time-series from an xdmf file
 * @param filename the path to the file to read
 * @return A TimeSeries object
 **/
TimeSeries XdmfReader::readTimeSeries(std::string filename) {
  // Extracting path to have it when needed for the hdf5 files
  std::string path = filename.substr(0, filename.rfind("/")+1);
  if (path == "")
    path = "./";

  // Opening xmf file
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(filename.c_str());

  TimeSeries series;

  auto xdmf     = doc.child("Xdmf");
  auto domain   = xdmf.child("Domain");
  auto grid     = domain.child("Grid");
  for (auto includes : grid.children("xi:include")) {
    std::filesystem::path ref = includes.attribute("href").value();
    series.push_back(ref);
  }
  return series;
}

}