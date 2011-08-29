/** @file
 *
 * Everything required by heatbath's main()
 */
#ifndef _TKKAPPAH_
#define _TKKAPPAH_
//should only be included in main prog

#include "heatbath.h"
#include "gaugefield.h"
#include "gaugefield_heatbath.h"
#include "gaugefield_k_hybrid.h"

//couple of timers
usetimer totaltime;
usetimer inittime;
usetimer polytime;
usetimer plaqtime;
usetimer updatetime;
usetimer overrelaxtime;
usetimer copytime;

#endif /* _TKKAPPAH_ */
