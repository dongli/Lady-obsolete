#include "geomtk.h"
#include "TracerManager.h"

int main(int argc, const char *argv[])
{
    LADY_DOMAIN domain(3);
    lady::TracerManager tracerManager;

    tracerManager.init(domain, 10);
    return 0;
}
