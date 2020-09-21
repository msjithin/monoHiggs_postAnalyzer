import ROOT

class FFApplicationTool():
    def __init__(self,theFFDirectory,channel):
        self.theFFDirectory = theFFDirectory
        
        self.theRawFile = ROOT.TFile(theFFDirectory+"uncorrected_fakefactors_"+channel+".root")
        if self.theRawFile.IsZombie():            
            raise RuntimeError("Problem loading the files!")
        self.ff_qcd_0jet = self.theRawFile.Get("rawFF_"+channel+"_qcd_0jet")
        self.ff_qcd_0jet_unc1_up = self.theRawFile.Get("rawFF_"+channel+"_qcd_0jet_unc1_up")
        self.ff_qcd_0jet_unc1_down = self.theRawFile.Get("rawFF_"+channel+"_qcd_0jet_unc1_down")
        self.ff_qcd_0jet_unc2_up = self.theRawFile.Get("rawFF_"+channel+"_qcd_0jet_unc2_up")
        self.ff_qcd_0jet_unc2_down = self.theRawFile.Get("rawFF_"+channel+"_qcd_0jet_unc2_down")

        self.ff_qcd_1jet = self.theRawFile.Get("rawFF_"+channel+"_qcd_1jet")
        self.ff_qcd_1jet_unc1_up = self.theRawFile.Get("rawFF_"+channel+"_qcd_1jet_unc1_up")
        self.ff_qcd_1jet_unc1_down = self.theRawFile.Get("rawFF_"+channel+"_qcd_1jet_unc1_down")
        self.ff_qcd_1jet_unc2_up = self.theRawFile.Get("rawFF_"+channel+"_qcd_1jet_unc2_up")
        self.ff_qcd_1jet_unc2_down = self.theRawFile.Get("rawFF_"+channel+"_qcd_1jet_unc2_down")

        self.ff_qcd_2jet = self.theRawFile.Get("rawFF_"+channel+"_qcd_2jet")
        self.ff_qcd_2jet_unc1_up = self.theRawFile.Get("rawFF_"+channel+"_qcd_2jet_unc1_up")
        self.ff_qcd_2jet_unc1_down = self.theRawFile.Get("rawFF_"+channel+"_qcd_2jet_unc1_down")
        self.ff_qcd_2jet_unc2_up = self.theRawFile.Get("rawFF_"+channel+"_qcd_2jet_unc2_up")
        self.ff_qcd_2jet_unc2_down = self.theRawFile.Get("rawFF_"+channel+"_qcd_2jet_unc2_down")

        self.ff_w_0jet = self.theRawFile.Get("rawFF_"+channel+"_w_0jet")
        self.ff_w_0jet_unc1_up = self.theRawFile.Get("rawFF_"+channel+"_w_0jet_unc1_up")
        self.ff_w_0jet_unc1_down = self.theRawFile.Get("rawFF_"+channel+"_w_0jet_unc1_down")
        self.ff_w_0jet_unc2_up = self.theRawFile.Get("rawFF_"+channel+"_w_0jet_unc2_up")
        self.ff_w_0jet_unc2_down = self.theRawFile.Get("rawFF_"+channel+"_w_0jet_unc2_down")

        self.ff_w_1jet = self.theRawFile.Get("rawFF_"+channel+"_w_1jet")
        self.ff_w_1jet_unc1_up = self.theRawFile.Get("rawFF_"+channel+"_w_1jet_unc1_up")
        self.ff_w_1jet_unc1_down = self.theRawFile.Get("rawFF_"+channel+"_w_1jet_unc1_down")
        self.ff_w_1jet_unc2_up = self.theRawFile.Get("rawFF_"+channel+"_w_1jet_unc2_up")
        self.ff_w_1jet_unc2_down = self.theRawFile.Get("rawFF_"+channel+"_w_1jet_unc2_down")

        self.ff_w_2jet = self.theRawFile.Get("rawFF_"+channel+"_w_2jet")
        self.ff_w_2jet_unc1_up = self.theRawFile.Get("rawFF_"+channel+"_w_2jet_unc1_up")
        self.ff_w_2jet_unc1_down = self.theRawFile.Get("rawFF_"+channel+"_w_2jet_unc1_down")
        self.ff_w_2jet_unc2_up = self.theRawFile.Get("rawFF_"+channel+"_w_2jet_unc2_up")
        self.ff_w_2jet_unc2_down = self.theRawFile.Get("rawFF_"+channel+"_w_2jet_unc2_down")

        self.ff_tt_0jet = self.theRawFile.Get("mc_rawFF_"+channel+"_tt")
        self.ff_tt_0jet_unc1_up = self.theRawFile.Get("mc_rawFF_"+channel+"_tt_unc1_up")
        self.ff_tt_0jet_unc1_down = self.theRawFile.Get("mc_rawFF_"+channel+"_tt_unc1_down")
        self.ff_tt_0jet_unc2_up = self.theRawFile.Get("mc_rawFF_"+channel+"_tt_unc2_up")
        self.ff_tt_0jet_unc2_down = self.theRawFile.Get("mc_rawFF_"+channel+"_tt_unc2_down")

        self.theFMvisFile = ROOT.TFile(theFFDirectory+"FF_corrections_1.root")
        if self.theFMvisFile.IsZombie():
            raise RuntimeError("Problem loading the files!")
        
        #
        self.mVisClosure_QCD_0jet = self.theFMvisFile.Get("closure_mvis_"+channel+"_0jet_qcd")
        self.mVisClosure_QCD_1jet = self.theFMvisFile.Get("closure_mvis_"+channel+"_1jet_qcd")
        self.mVisClosure_QCD_2jet = self.theFMvisFile.Get("closure_mvis_"+channel+"_2jet_qcd")
        self.mVisClosure_W_0jet = self.theFMvisFile.Get("closure_mvis_"+channel+"_0jet_w")
        self.mVisClosure_W_1jet = self.theFMvisFile.Get("closure_mvis_"+channel+"_1jet_w")
        self.mVisClosure_W_2jet = self.theFMvisFile.Get("closure_mvis_"+channel+"_2jet_w")
        self.mVisClosure_TT = self.theFMvisFile.Get("closure_mvis_"+channel+"_ttmc")        

        #MT and OSSS closure
        self.theFOSSSClosureFile = ROOT.TFile(theFFDirectory+"FF_QCDcorrectionOSSS.root")
        if self.theFOSSSClosureFile.IsZombie():
            raise RuntimeError("Problem loading the files!")
        
        self.OSSSClosure_QCD = self.theFOSSSClosureFile.Get("closure_OSSS_mvis_"+channel+"_qcd")
        self.OSSSClosure_QCD_unc1_up = self.theFOSSSClosureFile.Get("closure_OSSS_mvis_"+channel+"_qcd_unc1_up")
        self.OSSSClosure_QCD_unc1_down = self.theFOSSSClosureFile.Get("closure_OSSS_mvis_"+channel+"_qcd_unc1_down")
        self.OSSSClosure_QCD_unc2_up = self.theFOSSSClosureFile.Get("closure_OSSS_mvis_"+channel+"_qcd_unc2_up")
        self.OSSSClosure_QCD_unc2_down = self.theFOSSSClosureFile.Get("closure_OSSS_mvis_"+channel+"_qcd_unc2_down")

        self.MTClosure_W = self.theFOSSSClosureFile.Get("closure_mt_"+channel+"_w")
        self.MTClosure_W_unc1_up = self.theFOSSSClosureFile.Get("closure_mt_"+channel+"_w_unc1_up")
        self.MTClosure_W_unc1_down = self.theFOSSSClosureFile.Get("closure_mt_"+channel+"_w_unc1_down")
        self.MTClosure_W_unc2_up = self.theFOSSSClosureFile.Get("closure_mt_"+channel+"_w_unc2_up")
        self.MTClosure_W_unc2_down = self.theFOSSSClosureFile.Get("closure_mt_"+channel+"_w_unc2_down")

    def get_raw_FF(self,pt,fct):
        ff=1.0
        ff=fct.Eval(pt)    
        #if(pt>80):
        #    ff=fct.Eval(80)
        return ff

    def get_mvis_closure(self,mvis,fct):
        corr = 1.0
        corr = fct.Eval(mvis)
        if(mvis>350):
            corr=fct.Eval(350)
        #if (mvis<50):
        #    corr=fct.Eval(50)
        return corr

    def get_mt_closure(self,mt, fct):
        corr=1.0
        corr=fct.Eval(mt)
        #if (mt>120):
        #    corr=fct.Eval(120)
        return corr

    def get_ff(self, pt, mt, mvis, njets, frac_tt, frac_qcd, frac_w, unc='',upOrDown=''):
        ff_qcd = 1.0
        ff_w = 0
        ff_tt = 1.0
    
        #Raw ff
        if(njets==0):
            ff_qcd=self.get_raw_FF(pt,self.ff_qcd_0jet)
            if unc == 'ff_qcd_0jet_unc1':
                if upOrDown == 'up':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_0jet_unc1_up)
                elif upOrDown == 'down':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_0jet_unc1_down)
            elif unc == 'ff_qcd_0jet_unc2':
                if upOrDown == 'up':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_0jet_unc2_up)
                elif upOrDown == 'down':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_0jet_unc2_down)
            ff_w=self.get_raw_FF(pt,self.ff_w_0jet)
            if unc == 'ff_w_0jet_unc1':
                if upOrDown == 'up':
                    ff_w=self.get_raw_FF(pt,self.ff_w_0jet_unc1_up)
                elif upOrDown == 'down':
                    ff_w=self.get_raw_FF(pt,self.ff_w_0jet_unc1_down)
            elif unc == 'ff_w_0jet_unc2':
                if upOrDown == 'up':
                    ff_w=self.get_raw_FF(pt,self.ff_w_0jet_unc2_up)
                elif upOrDown == 'down':
                    ff_w = self.get_raw_FF(pt,self.ff_w_0jet_unc2_down)
        elif(njets==1):
            ff_qcd=self.get_raw_FF(pt,self.ff_qcd_1jet)
            if unc == 'ff_qcd_1jet_unc1':
                if upOrDown == 'up':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_1jet_unc1_up)
                elif upOrDown == 'down':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_1jet_unc1_down)
            elif unc == 'ff_qcd_1jet_unc2':
                if upOrDown == 'up':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_1jet_unc2_up)
                elif upOrDown == 'down':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_1jet_unc2_down)
            ff_w=self.get_raw_FF(pt,self.ff_w_1jet)
            if unc == 'ff_w_1jet_unc1':
                if upOrDown == 'up':
                    ff_w=self.get_raw_FF(pt,self.ff_w_1jet_unc1_up)
                elif upOrDown == 'down':
                    ff_w=self.get_raw_FF(pt,self.ff_w_1jet_unc1_down)
            elif unc == 'ff_w_1jet_unc2':
                if upOrDown == 'up':
                    ff_w=self.get_raw_FF(pt,self.ff_w_1jet_unc2_up)
                elif upOrDown == 'down':
                    ff_w = self.get_raw_FF(pt,self.ff_w_1jet_unc2_down)
        else:
            ff_qcd=self.get_raw_FF(pt,self.ff_qcd_2jet)
            if unc == 'ff_qcd_2jet_unc1':
                if upOrDown == 'up':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_2jet_unc1_up)
                elif upOrDown == 'down':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_2jet_unc1_down)
            elif unc == 'ff_qcd_2jet_unc2':
                if upOrDown == 'up':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_2jet_unc2_up)
                elif upOrDown == 'down':
                    ff_qcd = self.get_raw_FF(pt,self.ff_qcd_2jet_unc2_down)
            ff_w=self.get_raw_FF(pt,self.ff_w_2jet)
            if unc == 'ff_w_2jet_unc1':
                if upOrDown == 'up':
                    ff_w=self.get_raw_FF(pt,self.ff_w_2jet_unc1_up)
                elif upOrDown == 'down':
                    ff_w=self.get_raw_FF(pt,self.ff_w_2jet_unc1_down)
            elif unc == 'ff_w_2jet_unc2':
                if upOrDown == 'up':
                    ff_w=self.get_raw_FF(pt,self.ff_w_2jet_unc2_up)
                elif upOrDown == 'down':
                    ff_w = self.get_raw_FF(pt,self.ff_w_2jet_unc2_down)

        ff_tt=self.get_raw_FF(pt,self.ff_tt_0jet)
        if unc == 'ff_tt_0jet_unc1':
            if upOrDown == 'up':
                ff_tt=self.get_raw_FF(pt,self.ff_tt_0jet_unc1_up)
            elif upOrDown == 'down':
                ff_tt=self.get_raw_FF(pt,self.ff_tt_0jet_unc1_down)
        elif unc == 'ff_tt_0jet_unc2':
            if upOrDown == 'up':
                ff_tt=self.get_raw_FF(pt,self.ff_tt_0jet_unc2_up)
            elif upOrDown == 'down':
                ff_tt=self.get_raw_FF(pt,self.ff_tt_0jet_unc2_down)

        #mvis closure
        if njets == 0:
            if unc == 'mvisclosure_qcd_0jet':
                if upOrDown == 'up':
                    ff_qcd = ff_qcd*(1+2*(self.get_mvis_closure(mvis,self.mVisClosure_QCD_0jet)-1))
                elif upOrDown == 'down':
                    ff_qcd = ff_qcd
            else:
                ff_qcd = ff_qcd*self.get_mvis_closure(mvis,self.mVisClosure_QCD_0jet)
                
            if unc == 'mvisclosurce_w_0jet':
                if upOrDown == 'up':
                    ff_w = ff_w*(1+2*(self.get_mvis_closure(mvis,self.mVisClosure_W_0jet)-1))
                elif upOrDown == 'down':
                    ff_w = ff_w
            else:
                ff_w = self.get_mvis_closure(mvis,self.mVisClosure_W_0jet)
        elif njets == 1:
            if unc == 'mvisclosure_qcd_0jet':
                if upOrDown == 'up':
                    ff_qcd = ff_qcd*(1+2*(self.get_mvis_closure(mvis,self.mVisClosure_QCD_1jet)-1))
                elif upOrDown == 'down':
                    ff_qcd = ff_qcd
            else:
                ff_qcd = ff_qcd*self.get_mvis_closure(mvis,self.mVisClosure_QCD_1jet)
                
            if unc == 'mvisclosurce_w_0jet':
                if upOrDown == 'up':
                    ff_w = ff_w*(1+2*(self.get_mvis_closure(mvis,self.mVisClosure_W_1jet)-1))
                elif upOrDown == 'down':
                    ff_w = ff_w
            else:
                ff_w = self.get_mvis_closure(mvis,self.mVisClosure_W_1jet)

        else:
            if unc == 'mvisclosure_qcd_0jet':
                if upOrDown == 'up':
                    ff_qcd = ff_qcd*(1+2*(self.get_mvis_closure(mvis,self.mVisClosure_QCD_2jet)-1))
                elif upOrDown == 'down':
                    ff_qcd = ff_qcd
            else:
                ff_qcd = ff_qcd*self.get_mvis_closure(mvis,self.mVisClosure_QCD_2jet)
                
            if unc == 'mvisclosurce_w_0jet':
                if upOrDown == 'up':
                    ff_w = ff_w*(1+2*(self.get_mvis_closure(mvis,self.mVisClosure_W_2jet)-1))
                elif upOrDown == 'down':
                    ff_w = ff_w
            else:
                ff_w = self.get_mvis_closure(mvis,self.mVisClosure_W_2jet)

        if unc == 'mvisclosure_tt':
            if upOrDown == 'up':
                ff_tt = ff_tt*(1+2*(self.get_mvis_closure(mvis,self.mVisClosure_TT)-1))
            elif upOrDown == 'down':
                ff_tt = ff_tt
        else:
            ff_tt = self.get_mvis_closure(mvis,self.mVisClosure_TT)
        
        #MT and OSSS corrections
        if unc == 'mtclosure_w_unc1':
            if upOrDown == 'up':
                ff_w = ff_w*self.get_mt_closure(mt,self.MTClosure_W_unc1_up)
            elif upOrDown == 'down':
                ff_w = ff_w*self.get_mt_closure(mt,self.MTClosure_W_unc1_down)
        elif unc == 'mtclosure_w_unc2':
            if upOrDown == 'up':
                ff_w = ff_w*self.get_mt_closure(mt,self.MTClosure_W_unc2_up)
            elif upOrDown == 'down':
                ff_w = ff_w*self.get_mt_closure(mt,self.MTClosure_W_unc2_down)
        else:
            ff_w = ff_w*self.get_mt_closure(mt,self.MTClosure_W)
        if unc == 'osssclosure_qcd_unc1':
            if upOrDown == 'up':
                ff_qcd = ff_qcd*(1+2*(self.get_mvis_closure(mvis,self.OSSSClosure_QCD)-1))
            elif upOrDown == 'down':
                ff_qcd = ff_qcd        
        else:
            ff_qcd = ff_qcd*self.get_mvis_closure(mvis,self.OSSSClosure_QCD)
        
        ff_cmb = frac_tt*ff_tt + frac_qcd*ff_qcd + frac_w*ff_w
        return ff_cmb
