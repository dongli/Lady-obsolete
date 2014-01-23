#ifndef __Lady_TracerSpeciesInfo__
#define __Lady_TracerSpeciesInfo__

#include "lady_commons.h"

namespace lady {

class TracerSpeciesInfo {
    string name;
    string units;
    string brief;
public:
    TracerSpeciesInfo();
    TracerSpeciesInfo(const string &name, const string &units,
                      const string &brief);
    ~TracerSpeciesInfo();

    const string &getName() { return name; }
    const string &getUnits() { return units; }
    const string &getBrief() { return brief; }
};

}

#endif