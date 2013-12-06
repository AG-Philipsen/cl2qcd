/** @file
 *
 * Everything required by gaugeobservable's main()
 */

#ifndef GAUGEOBSERVABLESEXECUTABLE_H_
#define GAUGEOBSERVABLESEXECUTABLE_H_

#include "generalExecutable.h"

class gaugeobservablesExecutable : public multipleConfigurationExecutable
{
public:
	gaugeobservablesExecutable(int argc, const char* argv[]);
};


#endif /* GAUGEOBSERVABLESEXECUTABLE_H_ */
