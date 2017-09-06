/** @file
 * Declaration of the physics::lattices::Spinorfield class
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _PHYSICS_LATTICES_SPINORFIELD_
#define _PHYSICS_LATTICES_SPINORFIELD_

#include "../../hardware/system.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/buffers/plain.hpp"
#include "../prng.hpp"
#include "scalar.hpp"
#include "../../common_header_files/types_fermions.h"
#include "latticesInterfaces.hpp"
#include "../interfacesHandler.hpp"
#include "../../hardware/lattices/spinorfield.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

template <class Lattice, typename Basetype> void pseudo_randomize(const Lattice* to, int seed);

/**
 * Representation of a gaugefield.
 */
class Spinorfield {
    public:
        /**
         * Construct a gaugefield based on the input-files of the system
         */
        Spinorfield(const hardware::System&, const SpinorfieldParametersInterface&,const bool place_on_host = false);

        /**
         * Release resources
         */
        virtual ~Spinorfield();

        /*
         * Spinorfields cannot be copied
         */
        Spinorfield& operator=(const Spinorfield&) = delete;
        Spinorfield(const Spinorfield&) = delete;
        Spinorfield() = delete;

        /**
         * Get the buffers containing the gaugefield state on the devices.
         */
        const std::vector<const hardware::buffers::Plain<spinor> *> get_buffers() const noexcept;

        /**
         * Apply Gamma5 to the Spinorfield
         */
        void gamma5() const;

        /**
         * Set Spinorfield to zero
         */
        void zero() const;

        /**
         * Set Spinorfield to zero
         */
        void cold() const;

        /**
         * Set Spinorfield to be gaussian.
         */
        void gaussian(const physics::PRNG& prng) const;

        /**
         * Update the halos of the spinorfield buffers.
         */
        void update_halo() const;

        /**
         * Get the number of elements.
         */
        unsigned get_elements() const noexcept;

    protected:
        /**
         * Allow (re-)creation of buffers by children
         */
        void fill_buffers();

        /**
         * Allow destruction of buffers by children
         */
        void clear_buffers();

    private:
        hardware::System const& system;
        //TODO: turn the following pointer in a reference
        const SpinorfieldParametersInterface& spinorfieldParametersInterface;
    hardware::lattices::Spinorfield spinorfield;
        const bool place_on_host;
        void import(const spinor * const host) const;

        friend hmc_complex scalar_product(const Spinorfield& left, const Spinorfield& right);
        friend hmc_float squarenorm(const Spinorfield& field);
        friend void saxpy(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const Spinorfield& y);
        friend void sax(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x);
        friend void saxsbypz(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const hmc_complex beta, const Spinorfield& y, const Spinorfield& z);
        friend void pseudo_randomize<Spinorfield, spinor>(const Spinorfield* to, int seed);
};

/**
 * Create n spinorfields.
 *
 * \param n The number of spinorfields to create
 */
std::vector<Spinorfield *> create_spinorfields(const hardware::System& system, const size_t n, physics::InterfacesHandler& interfacesHandler, const bool place_on_host = false);

/**
 * Release the given spinorfields
 */
void release_spinorfields(const std::vector<Spinorfield *> fields);

/**
 * Calculate the scalar product of two spinorfields.
 */
hmc_complex scalar_product(const Spinorfield& left, const Spinorfield& right);
void scalar_product(const Scalar<hmc_complex>* res, const Spinorfield& left, const Spinorfield& right);

/**
 * Calculate the squarenorm of the spinorfield
 */
hmc_float squarenorm(const Spinorfield& field);
void squarenorm(const Scalar<hmc_float>* res, const Spinorfield& field);

/**
 * Perform the BLAS operation saxpy.
 *
 * out = alpha * x + y
 */
void saxpy(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const Spinorfield& y);
void saxpy(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x, const Spinorfield& y);

/**
 * Perform the BLAS operation sax.
 *
 * out = alpha * x
 */
void sax(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x);
void sax(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x);

/**
 * Perform the BLAS operation saxsbypz.
 *
 * out = alpha * x + beta * y + z
 */
void saxsbypz(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const hmc_complex beta, const Spinorfield& y, const Spinorfield& z);
void saxsbypz(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x, const Scalar<hmc_complex>& beta, const Spinorfield& y, const Spinorfield& z);

void log_squarenorm(const std::string& msg, const physics::lattices::Spinorfield& x);

/**
 * Fill the given field with a window of the other field.
 *
 * A window means all buffers of the first field will be filled with the same content of one buffer of the second field.
 */
void fill_window(const Spinorfield* out, const Spinorfield& src, const size_t idx);

}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_ */
