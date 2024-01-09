#include "LAMPSTPC.h"
#include "LTPadPlane.h"

ClassImp(LAMPSTPC);

LAMPSTPC::LAMPSTPC()
{
    fName = "LAMPS-TPC";
    if (fDetectorPlaneArray==nullptr)
        fDetectorPlaneArray = new TObjArray();

    //for (auto type=0; type<fNumPSAType; ++type)
    //    fChannelAnalyzer[type] = nullptr;
    fChannelAnalyzer[0] = nullptr;
}

LAMPSTPC::LAMPSTPC(TString parName) : LAMPSTPC()
{
    AddPar(parName);
    Init();
}

LAMPSTPC::LAMPSTPC(LKParameterContainer* par) : LAMPSTPC()
{
    AddPar(par);
    Init();
}

bool LAMPSTPC::Init()
{
    LKDetector::Init();

    // Put intialization todos here which are not iterative job though event
    e_info << "Initializing LAMPSTPC" << std::endl;

    InitChannelAnalyzer();

    return true;
}

bool LAMPSTPC::InitChannelAnalyzer()
{
    /*
    if (fChannelAnalyzer[0]==nullptr)
        for (auto type=0; type<fNumPSAType; ++type)
        {
            TString parName1 = Form("LAMPSTPC/pulseFile/%s",fTypeNames[type].Data());
            TString parName2 = Form("LAMPSTPC/analysis/%s",fTypeNames[type].Data());

            TString pulseFileName = fPar -> GetParString(parName1);
            int dynamicRange = 4096;
            int threshold = 50;
            int tbStart = 1;
            int tbMax = 350;
            int iterMax = 15;
            double tbStepCut = 0.01;
            int tbStartCut = 330;
            double scaleTbStep = 0.2;
            int thresholdOneStep = 2;
            int numTbAcendingCut = 5;
            if (fPar -> CheckPar(parName2)) {
                dynamicRange = fPar -> GetParInt(parName2,0);
                threshold = fPar -> GetParInt(parName2,1);
                tbStart = fPar -> GetParInt(parName2,2);
                tbMax = fPar -> GetParInt(parName2,3);
                iterMax = fPar -> GetParInt(parName2,4);
                tbStepCut = fPar -> GetParDouble(parName2,5);
                tbStartCut = fPar -> GetParInt(parName2,6);
                scaleTbStep = fPar -> GetParDouble(parName2,7);
                thresholdOneStep = fPar -> GetParInt(parName2,8);
                numTbAcendingCut = fPar -> GetParInt(parName2,9);
            }
            fChannelAnalyzer[type] = new LKChannelAnalyzer();
            fChannelAnalyzer[type] -> SetPulse(pulseFileName);
            fChannelAnalyzer[type] -> SetDynamicRange(dynamicRange);
            fChannelAnalyzer[type] -> SetThreshold(threshold);
            fChannelAnalyzer[type] -> SetTbStart(tbStart);
            fChannelAnalyzer[type] -> SetTbMax(tbMax);
            fChannelAnalyzer[type] -> SetIterMax(iterMax);
            fChannelAnalyzer[type] -> SetTbStepCut(tbStepCut);
            fChannelAnalyzer[type] -> SetTbStartCut(tbStartCut);
            fChannelAnalyzer[type] -> SetScaleTbStep(scaleTbStep);
            fChannelAnalyzer[type] -> SetThresholdOneStep(thresholdOneStep);
            fChannelAnalyzer[type] -> SetNumTbAcendingCut(numTbAcendingCut);
        }
    */
    fChannelAnalyzer[0] = new LKChannelAnalyzer();
    return true;
}

void LAMPSTPC::Print(Option_t *option) const
{
    e_info << "LAMPSTPC" << std::endl;
}

bool LAMPSTPC::BuildGeometry()
{
    return true;
}

bool LAMPSTPC::BuildDetectorPlane()
{
    auto mm = new LTPadPlane;
    AddPlane(mm);
    return true;
}

bool LAMPSTPC::IsInBoundary(Double_t x, Double_t y, Double_t z)
{
    //if (x>-10 and x<10)
    //    return true;
    //return false;
    return true;
}

bool LAMPSTPC::GetEffectiveDimension(Double_t &x1, Double_t &y1, Double_t &z1, Double_t &x2, Double_t &y2, Double_t &z2)
{
    x1 = -750;
    x2 = +750;
    y1 = -750;
    y2 = +750;
    z1 = -250;
    z2 = 1350;
    return true;
}
