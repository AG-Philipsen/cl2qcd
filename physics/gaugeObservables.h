#include "./lattices/gaugefield.hpp"
#include <fstream>
#include <cmath>
#include "../meta/inputparameters.hpp"
#include "algorithms/kappa_clover.hpp"

namespace physics{

  class gaugeObservables{
  public:

    /**
     * Measures all gauge observables according to parameter settings
     */
    void measureGaugeObservables(physics::lattices::Gaugefield& gf, int iteration, meta::Inputparameters params);

    /**
     * Measures plaquette and polyakov loop and writes them to file
     */
    void measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, const std::string& filename, meta::Inputparameters params);
    void measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params);

    void measureTransportcoefficientKappa(physics::lattices::Gaugefield& gaugefield, int iteration, meta::Inputparameters parameters);
    
    void writeTransportcoefficientKappaToFile(hmc_float kappa, std::string filename, int iteration);
    
    void writeTransportcoefficientKappaToFileUsingOpenOutputStream(hmc_float kappa, int iteration);



  private:
    std::ofstream outputToFile;

    void writePlaqAndPolyToFile(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter,  const std::string& filename);
    
  };
}
