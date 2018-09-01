from ROOT import TROOT, TCanvas, TH1F, TChain

###################################

	#Load a file 

###################################

chain= TChain('DecayTreeTuple/DecayTree')
chain.Add('Tuples/DVNtuple_mumu_tmp_MU.root')
chain.Add('Tuples/DVNtuple_mumu_tmp_MD.root')
	#you can add as many files as needed

print 'Number of entries: ' , chain.GetEntries() 
	#prints in the command line the number of entries

###################################
	#Step 1: Draw a histogram of the D0 mass and save it to a pdf
###################################

h1 = TH1F("h1", "D0 mass", 50, 1700, 2000)
chain.Project("h1", 'D0_M') 
	#this fills the histogram h1 with variable D0_M

c1 = TCanvas( 'c1', 'c1', 200, 10, 700, 500 )
h1.Draw()
c1.Print("test_MD.pdf")

###################################
	#Step 2: Draw a normalised histogram of both the generated and reconstructed x-momentum of D0
###################################

# 	#first we need to load also the generated tree
# MCchain= TChain('MCDecayTreeTuple/MCDecayTree')
# MCchain.Add('Tuples/DVNtuple_mumu_tmp_MU.root')
# MCchain.Add('Tuples/DVNtuple_mumu_tmp_MD.root')

# h3 = TH1F("h3", " ", 50, -15000, 20000)
# chain.Project("h3", 'D0_PX')

# h2 = TH1F("h2", " ", 50, -15000, 20000)
# MCchain.Project("h2", 'D0_TRUEP_X')

# c2 = TCanvas( 'c2', 'c2', 200, 10, 700, 500 )
# h2.DrawNormalized() 		
# 	#we are drawing normalised histograms (meaning we scale them with the number of events, so the integral is 1),
# 	#because h2 has so many more events, that otherwise we couldn't see h3 at all
# h3.SetLineColor(2)
# h3.DrawNormalized("same")
# c2.Print("test_P.pdf")
