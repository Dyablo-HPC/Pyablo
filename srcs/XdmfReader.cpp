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
Snapshot XdmfReader::readSnapshot(std::string filename) {
  // Extracting path to have it when needed for the hdf5 files
  std::string path = filename.substr(0, filename.rfind("/")+1);
  if (path == "")
    path = "./";
  
  // Opening xmf file
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(filename.c_str());

  Snapshot snap;

  std::map<std::string, hid_t> h5_handles;
  
  auto xdmf     = doc.child("Xdmf");
  auto domain   = xdmf.child("Domain");
  auto grid     = domain.child("Grid");
  auto time     = grid.child("Time");
  auto topology = grid.child("Topology");
  auto geometry = grid.child("Geometry");

  std::string time_str = time.attribute("Value").value();
  float time_value = std::strtod(time_str.c_str(), nullptr);
  snap.setTime(time_value);
  std::string grid_name = grid.attribute("Name").value();
  snap.setName(grid_name);

  /**
   * Connectivity info
   **/
  auto connectivity = topology.child("DataItem");
  int nCells {std::atoi(topology.attribute("NumberOfElements").value())};
  std::string top_type {topology.attribute("TopologyType").value()};
  if (top_type == "Quadrilateral")
    snap.setNDim(2);
  else if (top_type == "Hexahedron")
    snap.setNDim(3);
  else {
    std::cerr << "ERROR: Topology type not handled: " << top_type << std::endl;
    std::exit(1);
  }

  auto [conn_handle, conn_path] = splitH5Filename(connectivity.child_value());
  std::string conn_filename = path + conn_handle;
  snap.addH5Handle(conn_handle, conn_filename);
  snap.setConnectivity(conn_handle, conn_path, nCells);

  /**
   * Geometry info
   **/
  auto coordinates  = geometry.child("DataItem");
  std::string geom_type = geometry.attribute("GeometryType").value();
  if (geom_type != "XYZ") {
    std::cerr << "ERROR: Unknown geometry type: " << geom_type << std::endl;
    std::exit(1);
  }

  int nVertices = extractDimensions(coordinates.attribute("Dimensions").value());
  auto [coord_handle, coord_path] = splitH5Filename(coordinates.child_value());
  std::string coord_filename = path + coord_handle;
  snap.addH5Handle(coord_handle, coord_filename);
  snap.setCoordinates(coord_handle, coord_path, nVertices);

  /**
   * Attributes info
   **/
  for (auto att: grid.children("Attribute")) {
    std::string name   = att.attribute("Name").value();
    std::string center = att.attribute("Center").value();

    auto di = att.child("DataItem");
    std::string type = di.attribute("NumberType").value();
    
    auto [att_handle, att_path] = splitH5Filename(di.child_value());
    std::string att_filename = path + att_handle;
    snap.addH5Handle(att_handle, att_filename);
    snap.addAttribute(att_handle, att_path, name, type, center);  
  }
  return snap;
}

}