#define G__DICTIONARY

#include "HTTutilities/Jet2TauFakes/interface/IFunctionWrapper.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTFormula.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTGraph.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH1D.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH2F.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH2D.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH3D.h"
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"

#include <map>
#include <vector>
#include <string>


namespace 
{
    struct HTTutilities_Jet2TauFakes 
    {
        IFunctionWrapper ifctw;
        WrapperTGraph wtgr;
        WrapperTFormula wtfo;
        WrapperTH1D wth1d;
        WrapperTH2F wth2f;
        WrapperTH2D wth2d;
        WrapperTH3D wth3d;

        std::map<std::string, std::vector<size_t>> m1;
        std::map<std::string, std::vector<std::vector<size_t>>> m2;

        FakeFactor ff;
    };
}
