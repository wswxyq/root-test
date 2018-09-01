from ROOT import TROOT, TCanvas, TH1F, TChain, TLine

#m_D_meson 		= 1869.65

chain= TChain('D2KSPiPiPi_DD/DecayTree')
chain.Add('/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root')

print 'Number of entries: ' , chain.GetEntries() 
	#prints in the command line the number of entries

trigger='((D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
(D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1)) && ((D_Hlt1TrackMVADecision_TOS==1) || \
(D_Hlt1TwoTrackMVADecision_TOS==1))'

SignalWindow = '(D_M > 1839.65) && (D_M < 1899.65)'
SideBand = '((D_M > 1779.65) && (D_M < 1809.65)) || ((D_M > 1929.65) && (D_M < 1959.65))'

h1='(((D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
(D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1)) && ((D_Hlt1TrackMVADecision_TOS==1) || \
(D_Hlt1TwoTrackMVADecision_TOS==1)))&&((D_M > 1839.65) && (D_M < 1899.65))'

h2='(((D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
(D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1)) && ((D_Hlt1TrackMVADecision_TOS==1) || \
(D_Hlt1TwoTrackMVADecision_TOS==1)))&&(((D_M > 1779.65) && (D_M < 1809.65)) \
|| ((D_M > 1929.65) && (D_M < 1959.65)))'

#D_ReFit_P
h1_D_ReFit_P = TH1F("h1_D_ReFit_P", "D_ReFit_P", 100, 50000, 300000)
	#from 0 to 600000, 100 bins
chain.Project("h1_D_ReFit_P", 'D_ReFit_P', h1) 
	#this fills the histogram h1 with variable D_ReFit_P
c1_D_ReFit_P = TCanvas( 'c1_D_ReFit_P', 'c1_D_ReFit_P', 200, 10, 700, 500 )
h1_D_ReFit_P.SetFillColor(5)
h1_D_ReFit_P.SetFillStyle(3144)
h1_D_ReFit_P.Draw()
c1_D_ReFit_P.Print("D_ReFit_P1.pdf")

h2_D_ReFit_P = TH1F("h2_D_ReFit_P", "D_ReFit_P", 100, 50000, 300000)
	#from 0 to 600000, 100 bins
chain.Project("h2_D_ReFit_P", 'D_ReFit_P', h2) 
	#this fills the histogram h1 with variable D_ReFit_P
c2_D_ReFit_P = TCanvas( 'c2_D_ReFit_P', 'c2_D_ReFit_P', 200, 10, 700, 500 )
h2_D_ReFit_P.SetFillColor(8)
h2_D_ReFit_P.SetFillStyle(3005)
h2_D_ReFit_P.Draw()
c2_D_ReFit_P.Print("D_ReFit_P2.pdf")

h3_D_ReFit_P = h1_D_ReFit_P
h3_D_ReFit_P.Add(h2_D_ReFit_P, -1)
c3_D_ReFit_P = TCanvas( 'c3_D_ReFit_P', 'c3_D_ReFit_P', 200, 10, 700, 500 )
h3_D_ReFit_P.SetFillColor(4)
h3_D_ReFit_P.SetFillStyle(3016)
h3_D_ReFit_P.Draw()
c3_D_ReFit_P.Print("D_ReFit_P3.pdf")

c4_D_ReFit_P = TCanvas( 'c4_D_ReFit_P', 'c4_D_ReFit_P', 200, 10, 700, 500 )
h2_D_ReFit_P.SetFillColor(3)
h2_D_ReFit_P.SetFillStyle(3335)
h3_D_ReFit_P.SetFillColor(5)
h3_D_ReFit_P.SetFillStyle(3353)
h2_D_ReFit_P.DrawNormalized()
h3_D_ReFit_P.DrawNormalized('same')
c4_D_ReFit_P.Print("D_ReFit_P4.pdf")







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
