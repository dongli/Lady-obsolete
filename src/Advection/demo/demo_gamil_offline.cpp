#include "lady.h"

int main(int argc, char const *argv[])
{
    geomtk::TimeManager timeManager;
    lady::GamilAdaptor gamilAdaptor;
    lady::AdvectionManager advectionManager;
    geomtk::StampString o1("tracers.", ".nc"), o2("velocity.", ".nc");

//    gamilAdaptor.init("");

    return 0;
}