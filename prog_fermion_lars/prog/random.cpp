#include "random.h"

Random zufall(seed);

//Zufallszahl vom Typ hmc_float zwischen 0 und 1 
double random_01 ()
{
  return zufall.doub();
}

/*
double random_01_halboffen ()
{
  hmc_float temp;
  do
  {
    temp = random_01();
  }
  while
    (toBool(operator==(temp, hmc_float(0))));

  return temp;
};
*/

//Zufallszahl 1,2,3 vom Typ int 
int random_123 ()
{
  return zufall.int64() % 3 +1;
}


//Zufallszahl 0,1,2,3 vom Typ int
int random_0123 ()
{
  return zufall.int64() % 4;
}

/*
//Zufallszahl zwischen -1 und 1
double random_11 ()
{
  return hmc_float(2) * zufall.doub() - hmc_float(1);
};
*/

/*
//Zufälliger Gitterpunkt auf dem Gitter
multi1d<int> random_gitter(const multi1d<int>& latt_size)
{
  multi1d<int> temp (Nd);
  
  for (int i=0; i< temp.size(); i++)
  {
    temp [i] = zufall.int64() % latt_size[i];
  };

  return temp;
};
*/

//Zufallszahl 0,1 vom Typ int 
int random_0_1 ()
{
  return zufall.int64() % 2 ;
}

//Gibt Zufallszahl zwischen 0 und 2*Pi
double random_02pi ()
{
  hmc_float temp = hmc_float(2)*M_PI;
  return  temp - random_01()*temp;		//a -> a(0-z)+z
}


//Gibt drei Zufallszahlen 1,2,3
void random_1_2_3 (int rand[3])
{
  rand[0] = random_123();
  
  do
    {rand[1] = random_123();}
  while (rand[1] == rand[0]);

  if (rand[0]!=1 && rand[1]!=1)
    rand[2] = 1;
  else if (rand[0]!=2 && rand[1]!=2)
    rand[2] = 2;
  else
    rand[2] = 3;
}
