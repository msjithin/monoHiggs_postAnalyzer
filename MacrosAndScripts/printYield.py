import csv


# with open('eventYield.csv', mode='w') as yield_file:
#   yield_write = csv.writer(yield_file, delimiter=',', quotechar='"')
#   yield_write.writerow(['Name', 'Yield' ])

if histoname=='cutflow_n' :
  with open('eventYield.csv', mode='a') as yield_file:
    yield_write = csv.writer(yield_file, delimiter=',', quotechar='"')
    yield_write.writerow(['From Cutflow bin 8' ])
    yield_write.writerow(['Data_hist', Data_hist.GetBinContent(8) ])
    yield_write.writerow(['ZTT',  ZTT_hist.GetBinContent(8)  ])
    yield_write.writerow(['ZLL',  ZLL_hist.GetBinContent(8)  ])
    yield_write.writerow(['Fake', F_bkg.GetBinContent(8) ])
    yield_write.writerow(['TT',   TT_hist.GetBinContent(8) ])
    yield_write.writerow(['ggh',  GluGluH_hist.GetBinContent(8) ])
    yield_write.writerow(['VV' ,  VV_hist.GetBinContent(8) ])
if histoname=='tauCharge_6' :
  with open('eventYield.csv', mode='a') as yield_file:
    yield_write = csv.writer(yield_file, delimiter=',', quotechar='"')
    yield_write.writerow(['From tau charge, integral, with mT cut' ])
    yield_write.writerow(['Data_hist', Data_hist.Integral() ])
    yield_write.writerow(['ZTT',  ZTT_hist.Integral()  ])
    yield_write.writerow(['ZLL',  ZLL_hist.Integral()  ])
    yield_write.writerow(['Fake', F_bkg.Integral() ])
    yield_write.writerow(['TT',   TT_hist.Integral() ])
    yield_write.writerow(['ggh',  GluGluH_hist.Integral() ])
    yield_write.writerow(['VV' ,  VV_hist.Integral() ])
