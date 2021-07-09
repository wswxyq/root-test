#include <iostream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "subtr_exec.h"
using namespace std;

int main()
{
    /* code */
    int i1;
    int i2;
    cout << "Plot subtr_D2KSPiPiPi_DD? ('1' for yes, others for no)"<<"\n";
    cin >> i1;
    cout << "Plot subtr_D2KSKpPiPi_DD? ('1' for yes, others for no)"<<"\n";
    cin >> i2;

    if (i1==1) {
        subtr_D2KSPiPiPi_DD();
    }

    if (i2==1) {
        subtr_D2KSKpPiPi_DD();
    }

    return 0;
}
