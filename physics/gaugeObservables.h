#include "./lattices/gaugefield.hpp"
#include <fstream>
#include <cmath>

namespace physics{

  class gaugeObservables{
  public:
    
    /**
     * Measures plaquette and polyakov loop and writes them to file
     */
    void measurePlaqAndPoly(const physics::lattices::Gaugefield& gf, int iter, const std::string& filename, meta::Inputparameters params);
    void measurePlaqAndPoly(const physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params);

  private:
    void writePlaqAndPolyToFile(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter,  const std::string& filename);
    
  };
}
