from ROOT import TROOT, TCanvas, TH1F, TChain, TLine
	#constants
m_D_meson 		= 1869.65
m_D_meson_down 	= m_D_meson-100.0
m_D_meson_up 	= m_D_meson+100.0
err=30
halferr=err/2.0

	# load files
chain= TChain('D2KSPiPiPi_DD/DecayTree')
chain.Add('DV_Prompt_D2Kshhh_2016_DT_0.root')
	#you can add as many files as needed
print 'Number of entries: ' , chain.GetEntries() 
	#prints in the command line the number of entries


L0_set=('D_L0HadronDecision_TOS==1' or 'D_L0MuonDecision_Tis==1' or 'D_L0DiMuonDecision_Tis==1' or \
'D_L0ElectronDecision_Tis==1' or 'D_L0PhotonDecision_Tis==1')

Hlt1_set=('D_Hlt1TrackMVADecision_TOS==1' or 'D_Hlt1TwoTrackMVADecision_TOS==1')

trigger=L0_set and Hlt1_set

SignalWindow = 'D_M > 1839.65 ' and 'D_M < 1899.65 '
SideBand = ('D_M > 1754.65 ' and 'D_M < 1784.65 ') \
or ('D_M > 1954.65' and 'D_M < 1984.65')


#D_PT
h1_D_PT = TH1F("h1_D_PT", "D_PT", 100, 0, 20000)
	#from 0 to 600000, 100 bins
chain.Project("h1_D_PT", 'D_PT', trigger and SignalWindow) 
	#this fills the histogram h1 with variable D_PT
c1_D_PT = TCanvas( 'c1_D_PT', 'c1_D_PT', 200, 10, 700, 500 )
h1_D_PT.SetFillColor(5)
h1_D_PT.SetFillStyle(3144)
h1_D_PT.Draw()
c1_D_PT.Print("D_PT1.pdf")

h2_D_PT = TH1F("h2_D_PT", "D_PT", 100, 0, 20000)
	#from 0 to 600000, 100 bins
chain.Project("h2_D_PT", 'D_PT', trigger and SideBand) 
	#this fills the histogram h1 with variable D_PT
c2_D_PT = TCanvas( 'c2_D_PT', 'c2_D_PT', 200, 10, 700, 500 )
h2_D_PT.SetFillColor(8)
h2_D_PT.SetFillStyle(3005)
h2_D_PT.Draw()
c2_D_PT.Print("D_PT2.pdf")

h3_D_PT = h1_D_PT
h3_D_PT.Add(h2_D_PT, -1)
c3_D_PT = TCanvas( 'c3_D_PT', 'c3_D_PT', 200, 10, 700, 500 )
h3_D_PT.SetFillColor(4)
h3_D_PT.SetFillStyle(3016)
h3_D_PT.Draw()
c3_D_PT.Print("D_PT3.pdf")

c4_D_PT = TCanvas( 'c4_D_PT', 'c4_D_PT', 200, 10, 700, 500 )
h2_D_PT.SetFillColor(3)
h2_D_PT.SetFillStyle(3335)
h3_D_PT.SetFillColor(5)
h3_D_PT.SetFillStyle(3353)
h2_D_PT.DrawNormalized()
h3_D_PT.DrawNormalized('same')
c4_D_PT.Print("D_PT4.pdf")


'''
#D_ReFit_decayLength/D_ReFit_decayLengthErr
h2 = TH1F("h2", "D_ReFit_decayLength/D_ReFit_decayLengthErr", 100, 0, 200)
	#from 0 to 600000, 100 bins
chain.Project("h2", 'D_ReFit_decayLength/D_ReFit_decayLengthErr', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable D_ReFit_decayLength
c2 = TCanvas( 'c2', 'c2', 200, 10, 700, 500 )
h2.Draw()
c2.Print("D_ReFit_decayLength__Err.pdf")


#D_ReFit_P
h3 = TH1F("h3", "D_ReFit_P", 100, 0, 500000)
	#from 0 to 600000, 100 bins
chain.Project("h3", 'D_ReFit_P', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable D_ReFit_P
c3 = TCanvas( 'c3', 'c3', 200, 10, 700, 500 )
h3.Draw()
c3.Print("D_ReFit_P.pdf")


#D_ReFit_chi2
h4 = TH1F("h4", "D_ReFit_chi2", 100, 0, 100)
	#from 0 to 600000, 100 bins
chain.Project("h4", 'D_ReFit_chi2', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable D_ReFit_chi2
c4 = TCanvas( 'c4', 'c4', 200, 10, 700, 500 )
h4.Draw()
c4.Print("D_ReFit_chi2.pdf")


#D_ReFit_KS0_M
h5 = TH1F("h5", "D_ReFit_KS0_M", 100, 400, 600)
	#from 0 to 600000, 100 bins
chain.Project("h5", 'D_ReFit_KS0_M', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable D_ReFit_KS0_M
c5 = TCanvas( 'c5', 'c5', 200, 10, 700, 500 )
h5.Draw()
c5.Print("D_ReFit_KS0_M.pdf")


#D_ReFit_KS0_P
h6 = TH1F("h6", "D_ReFit_KS0_P", 100, 0, 100000)
	#from 0 to 600000, 100 bins
chain.Project("h6", 'D_ReFit_KS0_P', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable D_ReFit_KS0_P
c6 = TCanvas( 'c6', 'c6', 200, 10, 700, 500 )
h6.Draw()
c6.Print("D_ReFit_KS0_P.pdf")


#P1_ProbNNpi
h7 = TH1F("h7", "P1_ProbNNpi", 100, 0.8, 1.1)
	#from 0 to 600000, 100 bins
chain.Project("h7", 'P1_ProbNNpi', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable P1_ProbNNpi
c7 = TCanvas( 'c7', 'c7', 200, 10, 700, 500 )
h7.Draw()
c7.Print("P1_ProbNNpi.pdf")


#P2_ProbNNpi
h8 = TH1F("h8", "P2_ProbNNpi", 100, 0.8, 1.1)
	#from 0 to 600000, 100 bins
chain.Project("h8", 'P2_ProbNNpi', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable P2_ProbNNpi
c8 = TCanvas( 'c8', 'c8', 200, 10, 700, 500 )
h8.Draw()
c8.Print("P2_ProbNNpi.pdf")


#P3_ProbNNpi
h9 = TH1F("h9", "P3_ProbNNpi", 100, 0, 1.1)
	#from 0 to 600000, 100 bins
chain.Project("h9", 'P3_ProbNNpi', L0_set and Hlt1_set) 
	#this fills the histogram h1 with variable P3_ProbNNpi
c9 = TCanvas( 'c9', 'c9', 200, 10, 700, 500 )
h9.Draw()
c9.Print("P3_ProbNNpi.pdf")
'''
