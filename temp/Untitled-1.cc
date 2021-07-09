----fit parameters for signal----
RooRealVar N_tot("N_tot", "N_tot",  500, 0.,1500);
RooRealVar ACP("ACP", "ACP", 0., -0.2, 0.2);
RooFormulaVar Nsig_plus("Nsig_plus", "0.5*@0*(1+@1)",RooArgList(N_tot, ACP));
RooFormulaVar Nsig_minus("Nsig_minus", "0.5*@0*(1-@1)",RooArgList(N_tot, ACP));

----fit parameters for bkg----
RooRealVar Nbkg_komb("Nbkg_komb","comb. bkg",2100.,0.,5000.) ;
RooRealVar A_komb("A_komb", "A_komb", 0., -0.2, 0.2);
RooFormulaVar Nkomb_plus("Nkomb_plus", "0.5*@0*(1+@1)",RooArgList(Nbkg_komb, A_komb));
RooFormulaVar Nkomb_minus("Nkomb_minus", "0.5*@0*(1-@1)",RooArgList(Nbkg_komb, A_komb));

---create RooAddPdf for the two datasets----
RooAddPdf  model_plus("model_plus","model_plus",RooArgList(sigPDF, bkgPDF), RooArgList(Nsig_plus, Nkomb_plus)) ;
RooAddPdf  model_minus("model_minus","model_minus",RooArgList(sigPDF, bkgPDF), RooArgList(Nsig_minus,Nkomb_minus)) ;

---create combined datasample for simultaneous fit----
RooCategory sample("sample","sample") ;
sample.defineType("plus") ;
sample.defineType("minus") ;

RooDataSet combData("combData","combined data",RooArgList(*mass),Index(sample),Import("plus",*dh1),Import("minus",*dh1m)) ;

here, dh1 and dh1m are the two datasets and mass is the variable in which you fit

---define simultaneous PDF----
RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
simPdf.addPdf(model_plus,"plus") ;
simPdf.addPdf(model_minus,"minus") ;

---do the sim. fit---
    simPdf.fitTo(combData) ;


---draw a component----
RooPlot* frame1 = mass->frame(Title(" "));
combData.plotOn(frame1,Cut("sample==sample::plus")) ;
simPdf.plotOn(frame1,Slice(sample,"plus"),ProjWData(sample,combData)) ;
simPdf.paramOn(frame1);
combData.statOn(frame1);