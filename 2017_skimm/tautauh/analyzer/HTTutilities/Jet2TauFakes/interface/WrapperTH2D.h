
#ifndef Jet2TauFakes_WrapperTH2D_h
#define Jet2TauFakes_WrapperTH2D_h

#include "HTTutilities/Jet2TauFakes/interface/IFunctionWrapper.h"

#include "TH2D.h"

class WrapperTH2D : public IFunctionWrapper
{

    public:
        WrapperTH2D():IFunctionWrapper() {};
        WrapperTH2D(const TH2D& h, const std::string& name):IFunctionWrapper(name),m_histo(h) {};
        virtual ~WrapperTH2D();

        double value(size_t size, const double* xs) override
        {
            if(size<2) return 0.;
            int bx = m_histo.GetXaxis()->FindBin(xs[0]);
            int by = m_histo.GetYaxis()->FindBin(xs[1]);
            // Don't use overflow bins
            if(bx>m_histo.GetNbinsX()) bx = m_histo.GetNbinsX();
            else if(bx==0) bx = 1;
            if(by>m_histo.GetNbinsY()) by = m_histo.GetNbinsY();
            else if(by==0) by = 1;
            return m_histo.GetBinContent(bx,by);
        }
        double value(const std::vector<double>& xs) override
        {
            return value(xs.size(), xs.data());
        }

    private:
        TH2D m_histo;


    //private:
        //ClassDef(WrapperTH2D,1)
};


#endif
