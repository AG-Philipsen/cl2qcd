#include "host_random.h"

//this should be avoidable by using the extern Random rnd
// Random zufall(seed);

//Zufallszahl 1,2,3 vom Typ int 
inline int random_123 ()
{
  //return zufall.int64() % 3 +1;
  return rnd.int64() % 3 +1;
}

//Gibt drei Zufallszahlen 1,2,3
void random_1_2_3 (int rand[3])
{
  rand[0] = random_123();
  do
    {rand[1] = random_123();}
  while (rand[1] == rand[0]);
  rand[2] = 6 - rand[1] - rand[0];
}

void init_random_seeds(Random random, hmc_ocl_ran * out, const int NUM, usetimer * timer){
  (*timer).reset();
  int dummy;
  //The initial values for the first the x, y and z components of the state should be > 128
  for(int i = 0; i<NUM; i++){
    do{
      dummy = random.int64();
    } while(dummy < 128);
    (out[i]).x = dummy;
    do{
      dummy = random.int64();
    } while(dummy < 128);
    (out[i]).y = dummy;
    do{
      dummy = random.int64();
    } while(dummy < 128);
    (out[i]).z = dummy;    
    dummy = random.int64();
    (out[i]).w = dummy;    
  }
  (*timer).add();
  return;
}

void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha)
{
  hmc_float delta;
  hmc_float a0 ;
  hmc_float eta ;
  do
  {
    delta = -log(rnd.doub())/alpha*pow(cos(2. * PI * rnd.doub()), 2.) -log(rnd.doub())/alpha;
    a0 = 1.-delta;
    eta = rnd.doub();
  }while ( (1.-0.5*delta) < eta*eta); 
  hmc_float phi = 2.*PI*rnd.doub();
  hmc_float theta = asin(2.*rnd.doub() - 1.);
  dst[0] = a0;
  dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
  dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
  dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}

