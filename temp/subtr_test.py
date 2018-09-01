from ROOT import TROOT, TCanvas, TH1F, TChain, TLine, TCut
	# load files
chain= TChain('D2KSPiPiPi_DD/DecayTree')
chain.Add('/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root')

print 'Number of entries: ' , chain.GetEntries() 
	#prints in the command line the number of entries




SW1 = TCut("D_M>1830")
SW2 = TCut("D_M<1900")
SW = TCut(SW1.GetTitle()+'&&'+SW2.GetTitle())
SW.Print()

#D_M
h0_D_M = TH1F("h0_D_M", "D_M", 100, 1600, 2200)
	#from 0 to 600000, 100 bins
#chain.Project('h0_D_M', 'D_M', SW.GetTitle())
chain.Project('h0_D_M', 'D_M', "((D_M>1600)&&(D_M<1700))||((D_M>1700)&&(D_M>1800))||((D_M>1800)&&(D_M>1900))")

	#this fills the histogram h1 with variable D_M
c0_D_M = TCanvas( 'c0_D_M', 'c0_D_M', 200, 10, 700, 500 )
h0_D_M.Draw()
c0_D_M.Print("D_M0.pdf")
