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

    /**
     * Measures rectangle products of link variables and writes it to file
     */
    void measureRectangles(physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params);

    /**
     * Measures transportcoefficient kappa and writes it to file
     */
    void measureTransportcoefficientKappa(physics::lattices::Gaugefield& gaugefield, int iteration, meta::Inputparameters parameters);

  private:
    double plaquette;
    double plaquette_temporal;
    double plaquette_spatial;
    double rectangles;
    double kappa;
    hmc_complex polyakov;
    std::ofstream outputToFile;

    void writePlaqAndPolyToFile(int iter,  const std::string& filename);

    void writeTransportcoefficientKappaToFile(std::string filename, int iteration);
    
    void writeTransportcoefficientKappaToFileUsingOpenOutputStream(int iteration);

    void writeRectanglesToFile(int iter, const std::string& filename);
  };
}
