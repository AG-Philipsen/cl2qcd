/*
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

/**
 * @attention This file is under refactoring. Its code is being copied somewhere else and then
 *            it is doubled around in other files. This file will be deleted at the end of the
 *            refactoring!
 *
 *            This is also the reson why for the moment it is commented out.
 */

/*
 * 
 * #include "lattices/gaugefield.hpp"
 * #include "lattices/staggeredfield_eo.hpp"
 * #include "../meta/util.hpp"
 * #include "../common_header_files/operations_complex.h"
 * #include "../hardware/device.hpp"
 * #include "../hardware/code/fermions_staggered.hpp"
 * #include "../hardware/code/spinors.hpp"
 * #include "../interfaceImplementations/latticesParameters.hpp"
 * #include <vector>
 * 
 * class TestGaugefield_stagg {
 * 
 *  public:
 *  TestGaugefield_stagg(const hardware::System * system) : system(system), prngParameters( &system->get_inputparameters() ), params(&system->get_inputparameters()),
 *                                                          prng(*system, &prngParameters), gf(*system, params, prng)
 *  {
 *     BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
 *     const auto & inputfile = system->get_inputparameters();
 *     meta::print_info_hmc(inputfile);
 *   };
 * 
 *   const hardware::code::Fermions_staggered * get_device();
 *   const hardware::buffers::SU3 * get_gaugefield();
 *   const hardware::code::Gaugefield* get_gf_code();
 *   void save_conf();
 * 
 *  private:
 *   const hardware::System * const system;
 *   const physics::ParametersPrng_fromMetaInputparameters prngParameters;
 *   const physics::lattices::GaugefieldParametersImplementation params;
 *   physics::PRNG prng;
 *   physics::lattices::Gaugefield gf; //I changed this variable from const to not const to be able to save the conf to a lime file!     
 * };
 * 
 *  
 * const hardware::code::Fermions_staggered* TestGaugefield_stagg::get_device()   
 * {
 *   return system->get_devices()[0]->getFermionStaggeredCode(); 
 * }
 * 
 * const hardware::code::Gaugefield* TestGaugefield_stagg::get_gf_code()      
 * {
 *   return system->get_devices()[0]->getGaugefieldCode();                                                                             
 * } 
 * 
 * const hardware::buffers::SU3 * TestGaugefield_stagg::get_gaugefield()
 * {
 *   return gf.get_buffers().at(0);                                                                                                      
 * }
 * 
 * void TestGaugefield_stagg::save_conf()
 * { 
 *   gf.save("conf_lime",13);                                                                                                            
 * }                                     
 * 
 * /---------------------------------------------------------------------------------/
 * 
 * //
 * // Fuction that "convert" a matrix to a string with a proper structure to be 
 * // written to the text file that will be later used for the reference code 
 * //
 * std::string matrix_to_string(Matrixsu3 m)
 * {
 *   std::ostringstream os;
 *   os.precision(16);
 *   os << "(" << m.e00.re << "," << m.e00.im << ") (" << m.e01.re << "," << m.e01.im << ") (" << m.e02.re << "," << m.e02.im << ")\n";
 *   os << "(" << m.e10.re << "," << m.e10.im << ") (" << m.e11.re << "," << m.e11.im << ") (" << m.e12.re << "," << m.e12.im << ")\n";
 *   os << "(" << m.e20.re << "," << m.e20.im << ") (" << m.e21.re << "," << m.e21.im << ") (" << m.e22.re << "," << m.e22.im << ")\n\n";
 *   return os.str();
 * }
 * 
 * //
 * // Fuction that "convert" a su3vec to a string with a proper structure to be   
 * // written to the text file that will be later used for the reference code     
 * //
 * std::string su3vec_to_string(su3vec m)
 * {
 *   std::ostringstream os;
 *   os.precision(16);
 *   os << "(" << m.e0.re << "," << m.e0.im << ") (" << m.e1.re << "," << m.e1.im << ") (" << m.e2.re << "," << m.e2.im << ")\n\n";
 *   return os.str();
 * }
 * 
 * //Tool to be used in the function print_gaugefield_to_textfile 
 * void get_full_coord_from_site_idx(int site_idx, int &x, int &y, int &z, int &t, const int ns)
 * {
 *   int volspace=ns*ns*ns;
 *   int space=site_idx%volspace;
 *   t=site_idx/volspace;
 *   z=space/ns/ns;
 *   int acc=z;
 *   y=space/ns-ns*acc;
 *   acc=ns*acc+y;
 *   x=space-ns*acc;
 * }
 * 
 * Matrixsu3 unit_matrixsu3()
 * {
 *   Matrixsu3 out;
 *   out.e00.re = 1.;
 *   out.e00.im = 0.;
 *   out.e01.re = 0.;
 *   out.e01.im = 0.;
 *   out.e02.re = 0.;
 *   out.e02.im = 0.;
 * 
 *   out.e10.re = 0.;
 *   out.e10.im = 0.;
 *   out.e11.re = 1.;
 *   out.e11.im = 0.;
 *   out.e12.re = 0.;
 *   out.e12.im = 0.;
 * 
 *   out.e20.re = 0.;
 *   out.e20.im = 0.;
 *   out.e21.re = 0.;
 *   out.e21.im = 0.;
 *   out.e22.re = 1.;
 *   out.e22.im = 0.;
 * 
 *   return out;
 * }
 * 
 * 
 * inline Matrixsu3 multiply_matrixsu3_by_complex (Matrixsu3 in, hmc_complex factor)
 * {
 *   Matrixsu3 out;
 *   out.e00 = complexmult(in.e00, factor);
 *   out.e01 = complexmult(in.e01, factor);
 *   out.e02 = complexmult(in.e02, factor);
 *   out.e10 = complexmult(in.e10, factor);
 *   out.e11 = complexmult(in.e11, factor);
 *   out.e12 = complexmult(in.e12, factor);
 *   out.e20 = complexmult(in.e20, factor);
 *   out.e21 = complexmult(in.e21, factor);
 *   out.e22 = complexmult(in.e22, factor);
 *   return out;
 * }
 * 
 * //
 * //  In the reference code the lattice is reorganized in the following way:    
 * //
 * //  links used according to this scheme
 * //   0             size            size2           size3         no_links  
 * //   |-------|-------|-------|-------|-------|-------|-------|-------|
 * //      e        o       e       o       e       o       e       o  
 * //        x-dir         y-dir         z-dir         t-dir  
 * //
 * //  where e=even, o=odd, whereas size=VOL4D.  
 * //  Hence, in order to use the same random configuration in tests 
 * //  I have to print all links to a text file according this scheme.  
 * //
 * //  @note: In our program mu=0 is the TIME direction and mu=1,2,3 are the x,y,z direction!!!   
 * // 
 * //
 * #include "../hardware/code/gaugefield.hpp"
 * void print_gaugefield_to_textfile(std::string outputfile, TestGaugefield_stagg * cpu, const meta::Inputparameters & params)
 * {
 *   int nt=params.get_ntime();
 *   int ns=params.get_nspace();
 *   if(ns!=nt){
 *     logger.fatal() << "The lattice must be isotropic to call the function print_gaugefield_to_textfile(...)!";
 *     abort();
 *   }
 *   //conf_old is the Matrixsu3 array with the links in the standard order (standard for this code)                                     
 *   //conf_new is the Matrixsu3 array in the right order (ref. code scheme) to be written to the file                                   
 *   Matrixsu3 *conf_old=new Matrixsu3[ns*ns*ns*nt*4];
 *   Matrixsu3 *conf_new=new Matrixsu3[ns*ns*ns*nt*4];
 *   cpu->get_gf_code()->exportGaugefield(conf_old,cpu->get_gaugefield());
 *   //Now I have conf_old and I have to fill properly conf_new                                                                          
 *   int x,y,z,t,num,even,size;
 *   size=ns*ns*ns*nt;
 *   for(int i=0; i<ns*ns*ns*nt; i++){
 *     get_full_coord_from_site_idx(i,x,y,z,t,ns);
 *     even = (x+y+z+t)%2;
 *     // even=0 for even sites                                                                                                          
 *     // even=1 for odd sites                                                                                                           
 *     num = even*size/2 + (x+y*ns+z*ns*ns+t*ns*ns*ns)/2;
 *     // num is where, in conf_new, conf_old[...] is to be written                                                                      
 *     conf_new[num       ]=conf_old[4*i+1]; //x-dir                                                                                     
 *     conf_new[num+size  ]=conf_old[4*i+2]; //y-dir                                                                                     
 *     conf_new[num+size*2]=conf_old[4*i+3]; //z-dir                                                                                     
 *     conf_new[num+size*3]=conf_old[4*i  ]; //t-dir      
 *     //Only to check if the function works correctly ---> Comment out these line if all is fine                                        
 *     //
 *     //  conf_new[num].e01={y,0.};                                                                                                         
 *     //  conf_new[num].e02={z,0.};                                                                                                         
 *     //  conf_new[num].e10={t,0.};                                                                                                         
 *     //  conf_new[num].e11={0.,0.};                                                                                                        
 *     //  conf_new[num].e12={0.,0.};                                                                                                        
 *     //  conf_new[num].e20={0.,0.};                                                                                                        
 *     //  conf_new[num].e21={num,0.};                                                                                                       
 *     //  conf_new[num].e22={even,0.};                                                                                                      
 *     //  conf_new[num+size].e00={x,0.};                                                                                                    
 *     //  conf_new[num+size].e01={y,0.};                                                                                                    
 *     //  conf_new[num+size].e02={z,0.};                                                                                                    
 *     //  conf_new[num+size].e10={t,0.};                                                                                                    
 *     //  conf_new[num+size].e11={0.,0.};                                                                                                   
 *     //  conf_new[num+size].e12={0.,0.};                                                                                                   
 *     //  conf_new[num+size].e20={0.,0.};                                                                                                   
 *     //  conf_new[num+size].e21={num,0.};                                                                                                  
 *     //  conf_new[num+size].e22={even,0.};                                                                                                 
 *     //  conf_new[num+size*2].e00={x,0.};                                                                                                  
 *     //  conf_new[num+size*2].e01={y,0.};                                                                                                  
 *     //  conf_new[num+size*2].e02={z,0.};                                                                                                  
 *     //  conf_new[num+size*2].e10={t,0.};                                                                                                  
 *     //  conf_new[num+size*2].e11={0.,0.};                                                                                                 
 *     //  conf_new[num+size*2].e12={0.,0.};                                                                                                 
 *     //  conf_new[num+size*2].e20={0.,0.};                                                                                                 
 *     //  conf_new[num+size*2].e21={num,0.};                                                                                                
 *     //  conf_new[num+size*2].e22={even,0.};                                                                                               
 *     //  conf_new[num+size*3].e00={x,0.};                                                                                                  
 *     //  conf_new[num+size*3].e01={y,0.};                                                                                                  
 *     //  conf_new[num+size*3].e02={z,0.};                                                                                                  
 *     //  conf_new[num+size*3].e10={t,0.};                                                                                                  
 *     //  conf_new[num+size*3].e11={0.,0.};                                                                                                 
 *     //  conf_new[num+size*3].e12={0.,0.};                                                                                                 
 *     //  conf_new[num+size*3].e20={0.,0.};                                                                                                 
 *     //  conf_new[num+size*3].e21={num,0.};                                                                                                
 *     //  conf_new[num+size*3].e22={even,0.};
 *     
 *   }
 *   //Now we can write conf_new to the file                                                                                             
 *   std::ofstream file(outputfile.c_str());
 *   file << ns << " " << ns << " " << ns << " " << nt << " ";
 *   file << params.get_beta() << " " << params.get_mass() << " 12345" << std::endl;
 *   //The last number that I set to 12345 should be the hmc iteration; here it is not relevant                                          
 *   for(int i=0; i<ns*ns*ns*nt*4; i++)
 *     file << matrix_to_string(conf_new[i]);
 *   file.close();
 *
 *   // //To keep only links in one direction: it is usefull only for debugging
 *   // for(int i=0; i<ns*ns*ns*nt*4; i++){
 *   //   if(i!=0)
 *   //     conf_old[i]=unit_matrixsu3();
 *   //
 *   //   //I multiply each link by the BC-phase to run the Ref.Code with periodic-BC
 *   //   //if((i%4)==0)
 *   //   //  conf_old[i]=multiply_matrixsu3_by_complex(conf_old[i],{1./sqrt(2.),1./sqrt(2.)});
 *   //
 *   //   // conf_old[i].e00.re=0.;
 *   //   // conf_old[i].e00.im=0.;
 *   //   // conf_old[i].e01.re=0.;
 *   //   // conf_old[i].e01.im=0.;
 *   //   // conf_old[i].e02.re=0.;
 *   //   // conf_old[i].e02.im=0.;
 *   //   // conf_old[i].e10.re=0.;
 *   //   // conf_old[i].e10.im=0.;
 *   //   // conf_old[i].e11.re=0.;
 *   //   // conf_old[i].e11.im=0.;
 *   //   // conf_old[i].e12.re=0.;
 *   //   // conf_old[i].e12.im=0.;
 *   //   // conf_old[i].e20.re=0.;
 *   //   // conf_old[i].e20.im=0.;
 *   //   // conf_old[i].e21.re=0.;
 *   //
 *   // }
 *   // cpu->get_gf_code()->importGaugefield(cpu->get_gaugefield(),conf_old);
 *   // cpu->save_conf();
 * }
 * 
 * //
 * // In the reference code the lattice is reorganized in the following way: 
 * //
 * //  0             size
 * //  |-------|-------|   
 * //      e       o 
 * //
 * // where e=even, o=odd, whereas size=VOL4D. 
 * // Hence, in order to use the same random staggered field in tests  
 * // I have to print it to a text file according this scheme.
 * //
 * 
 * void print_staggeredfield_to_textfile(std::string outputfile, su3vec * sf, const meta::Inputparameters & params)
 * {
 *   int nt=params.get_ntime();
 *   int ns=params.get_nspace();
 *   if(ns!=nt){
 *     logger.fatal() << "The lattice must be isotropic to call the function print_staggeredfield_to_textfile(...)!";
 *     abort();
 *   }
 *   //sf     is the su3vec array ordered with the "superindex scheme"                                                                     
 *   //sf_new is the su3vec array in the right order (ref. code scheme) to be written to the file                                          
 *   su3vec *sf_new = new su3vec[ns*ns*ns*nt];
 *   //Now I have conf_old and I have to fill properly conf_new                                                                            
 *   int x,y,z,t,num,even,size;
 *   size=ns*ns*ns*nt;
 *   for(int i=0; i<ns*ns*ns*nt; i++){
 *     get_full_coord_from_site_idx(i,x,y,z,t,ns);
 * 
 *     // logger.warn() << "(" << x << "," << y << "," << z << "," << t << ") => "                                                         
 *     //            << "(" << sf[i].e0.re << "," << sf[i].e0.im << ") ("                                                                  
 *     //            << sf[i].e1.re << "," << sf[i].e1.im << ") ("                                                                         
 *     //            << sf[i].e2.re << "," << sf[i].e2.im << ")";                                                                          
 * 
 *     even = (x+y+z+t)%2;
 *     // even=0 for even sites                                                                                                            
 *     // even=1 for odd sites                                                                                                             
 *     num = even*size/2 + (x+y*ns+z*ns*ns+t*ns*ns*ns)/2;
 *     // num is where, in conf_new, conf_old[...] is to be written                                                                        
 *     sf_new[num]=sf[i];
 *   }
 *   //Now we can write sf_new to the file 
 *   std::ofstream file(outputfile.c_str());
 *   file << ns << " " << ns << " " << ns << " " << nt << std::endl;
 *   for(int i=0; i<ns*ns*ns*nt; i++){
 *     get_full_coord_from_site_idx(i,x,y,z,t,ns);
 *     file << su3vec_to_string(sf_new[i]);
 *   }
 *   file.close();
 * }
 * 
 * void print_staggeredfield_eo_to_textfile(std::string outputfile, su3vec * sf, const meta::Inputparameters &  params)
 * {
 *   int nt=params.get_ntime();
 *   int ns=params.get_nspace();
 *   if(ns!=nt){
 *     logger.fatal() << "The lattice must be isotropic to call the function print_staggeredfield_to_textfile(...)!";
 *     abort();
 *   }
 *   //sf     is the su3vec array ordered with the "even-odd superindex scheme"                                                                     
 *   //sf_new is the su3vec array in the right order (ref. code scheme) to be written to the file
 *   // ======> hence sf_new is in this case equal to sf that contain the values of the field only in
 *   //         even (or odd) sites
 *   //We can write sf directly to the file 
 *   std::ofstream file(outputfile.c_str());
 *   file << ns << " " << ns << " " << ns << " " << nt << std::endl;
 *   for(int i=0; i<ns*ns*ns*nt/2; i++)
 *     file << su3vec_to_string(sf[i]);
 *   file.close();
 * }
 * 
 * void print_staggeredfield_eo_to_textfile(std::string outputfile, const physics::lattices::Staggeredfield_eo* sf, const hardware::System& system)
 * {
 *   const auto & params = system.get_inputparameters();//deprecated
 *   su3vec * out_sf;
 *   size_t NUM_ELEMENTS_SF_EO = 0;//@todo: this does not work anyway because of the fct. above! hardware::code::get_eoprec_spinorfieldsize(params);
 *   out_sf = new su3vec[NUM_ELEMENTS_SF_EO];
 *   auto sf_bufs = sf->get_buffers();
 *   if(sf_bufs.size() > 1){
 *     logger.fatal() << "Print staggeredfield to textfile not implemented for multi device!";
 *     abort();
 *   }
 *   sf_bufs[0]->dump(out_sf);
 *   print_staggeredfield_eo_to_textfile(outputfile,out_sf,params);
 * }
 * 
 * 
 * //
 * // Function that returns a vector with the 6 real number contained in an su3vec
 * //
 * std::vector<hmc_float> reals_from_su3vec(su3vec v){
 *   std::vector<hmc_float> out;
 *   out.push_back(v.e0.re);
 *   out.push_back(v.e0.im);
 *   out.push_back(v.e1.re);
 *   out.push_back(v.e1.im);
 *   out.push_back(v.e2.re);
 *   out.push_back(v.e2.im);
 *   return out;
 * }
 * 
 *   
 *   
 */
