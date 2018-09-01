def Sweights(self, ws, obs='x', data= 'data', model='model', pdfs = [''], nevs = ['']):
        '''
        Calcolate s-weights.
        ws = rooworkspace
        data = data name
        model = model name
        pdfs = pdfs to fix
        nevs = list of variables with nevs
        '''
        for pdf in pdfs: self.FixVariables(ws, pdf, True)
        ws.var(obs).setConstant(False)
        # Re fitting
        RooMsgService.instance().setSilentMode(True)
        ws.pdf('model').fitTo(ws.data(data),RooFit.Extended(True))
        # Make splot!
        sData = RooStats.SPlot("sData","An SPlot", ws.data(data), ws.pdf('model'), self.GetVarsList(ws, nevs) )
        RooMsgService.instance().reset()
        # Declare a tree
        tVars = {}
        for nev in nevs: tVars['w_'+nev] = 'F'
        t = DefineTree(tVars, 'ntp_sweights', 'tree with s-weights')
        # Fill the tree
        for entry in xrange(0,ws.data(data).numEntries()):
            Vars = ws.data(data).get(entry)
            for nev in nevs:
                t['vars']['w_'+nev][0] = Vars.find(nev+'_sw').getVal()
            t['tree'].Fill()
        return t

