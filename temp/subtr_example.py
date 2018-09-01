from ROOT import TROOT, TCanvas, TH1F, TChain, TLine


	# load files
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




