#include <bits/stdc++.h>

#include <omp.h>
#include "Snapshot.h"
#include "XdmfReader.h"

using namespace dyablo;

enum ToolMode {
  TM_Invalid,
  TM_TimeSeries,
  TM_Slice,
  TM_Profile
};

using RealSeries = std::vector<double>;

struct Case {
  ToolMode mode;
  std::filesystem::path path;
  
  Case() : mode(ToolMode::TM_Invalid) {};

  void init(int argc, char **argv) {
    std::string tool_type{argv[1]};
    path = argv[2];
    if (tool_type == "time-series")
      initTimeseries(argc, argv);
    else if (tool_type == "slice")
      initSlice(argc, argv);
    else if (tool_type == "profile")
      initProfile(argc, argv);
  }

  void initTimeseries(int argc, char **argv) {
    if (path.extension() != ".xmf") {
      std::cout << "ERROR ! Given file is not an .xmf file !" << std::endl;
      return;
    }

    if (argc < 4) {
      std::cout << "ERROR no filename given :" << std::endl;
      std::cout << "USAGE " << argv[0] << " time-series INPUT OUTPUT" << std::endl;
      return;
    }

    XdmfReader reader;
    TimeSeries ts = reader.readTimeSeries(path);

    RealSeries t, E, Ek, e, M;

    for (auto filename : ts) {
      auto file_path = path.replace_filename(filename); 
      std::cout << "Opening file " << file_path.string() << std::endl;
      auto snap = reader.readSnapshot(file_path);

      t.push_back(snap.getTime());
      E.push_back(snap.getTotalEnergy());
      Ek.push_back(snap.getTotalKineticEnergy());
      e.push_back(snap.getTotalInternalEnergy(1.666666666667));
      M.push_back(snap.getTotalMass());

      snap.close();
    }

    std::ofstream f_out(argv[3]);

    size_t nsnaps = ts.size();
    for (auto isnap=0; isnap < nsnaps; ++isnap)
      f_out << t[isnap] << " " << M[isnap] << " " << Ek[isnap] << " " << e[isnap] << " " << E[isnap] << std::endl;

    f_out.close();

    mode = ToolMode::TM_TimeSeries;
  }

  void initSlice(int argc, char **argv) {

    mode = ToolMode::TM_Slice;
  }

  void initProfile(int argc, char **argv) {

    mode = ToolMode::TM_Profile;
  }
};

void usage(int argc, char **argv) {
  std::cout << "USAGE: " << argv[0] << " MODE PATH [OPTIONS]" << std::endl;
  std::cout << "  MODE: time-series slice or profile" << std::endl;
  std::cout << "  PATH: path to the file or folder containing files" << std::endl;
}

Case parse_arguments(int argc, char **argv) {
  Case result;
  
  if (argc < 3)
    usage(argc, argv);
  else 
    result.init(argc, argv);

  return result;
}

int main(int argc, char **argv) {
  // Reading case from arguments
  Case cur_case = parse_arguments(argc, argv);
  if (cur_case.mode == ToolMode::TM_Invalid)
    return 1;





  return 0;
}