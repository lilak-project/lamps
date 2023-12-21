#include "ejungwooA.h"
#include "LKPulse.cpp"
#include "LKPulse.h"
#include "LKPulseAnalyzer.cpp"
#include "LKPulseAnalyzer.h"
#include "setChannels.h"

int fCvsDX = 1200;
int fCvsDY = 720;
//int fCvsDX = 3000;
//int fCvsDY = 2000;
int fCvsGroupDX = 3;
int fCvsGroupDY = 2;
int fNumCvsGroup = 0;
bool fDrawChannel = false;
bool fDrawOnlyGoodChannels = true;
bool fDrawAverage = true;
bool fDrawAccumulate = true;
bool fDrawProjY = false;
bool fDrawHeightWidth = false;
bool fDrawWidth = false;
bool fDrawHeight = false;
bool fDrawPulseTb = false;
bool fDrawPedestal = false;
bool fDrawResidual = false;
bool fDrawReference = false;
bool fDrawPulseParameter = false;
bool fWriteReferencePulse = false;
bool fWriteTree = false;
int  fCAAC = 0;
int  fMaxInputFiles = 1;
bool fSaveAll = false;
const char* fDataPath = "data";

void SetRunAll();
void SetAllFalse();
void SetDrawCheckAllResults();
void SetDrawCheckGoodChannels();
void SetExtractBuffer(int ecaac = 2011044);
void SetDrawResidual();
void SetDebugPulseData(int i);

void anaExtractPulse(int mode=1)
{
    if (mode==0) SetRunAll();
    if (mode==1) SetDrawCheckAllResults();
    if (mode==2) SetDrawCheckGoodChannels();
    if (mode==3) SetExtractBuffer();
    if (mode==4) SetDrawResidual();
    if (mode==5) SetDebugPulseData(0);
    if (mode==6) SetDebugPulseData(efCsI);

    auto run = new LKRun();
    run -> AddPar("config.mac");
    ifstream file_list_conv("/home/ejungwoo/data/texat/conv/list_conv");
    int countFiles = 0;
    TString fileName;
    while (file_list_conv>>fileName) {
        if (countFiles>=fMaxInputFiles) break;
        run -> AddInputFile(fileName);
        countFiles++;
    }

    run -> SetTag("read");
    run -> Init();
    auto channelArray = run -> GetBranchA("RawData");

    LKPulse *pulse[fNumTypes];
    LKPulseAnalyzer *ana[fNumTypes];
    for (auto type : fSelTypes)
    {
        //pulse[type] = new LKPulse(Form("data/pulseReference_%s.root",fTypeNames[type]));

        ana[type] = new LKPulseAnalyzer(fTypeNames[type],fDataPath);
        ana[type] -> SetTbRange(fPulseCuts[type][7],350);
        //ana[type] -> SetTbStart(fPulseCuts[type][7]);
        //ana[type] -> SetTbMax(350);
        ana[type] -> SetThreshold(fPulseCuts[type][0]);
        ana[type] -> SetPulseHeightCuts(fPulseCuts[type][1],fPulseCuts[type][2]);
        ana[type] -> SetPulseWidthCuts(fPulseCuts[type][3],fPulseCuts[type][4]);
        ana[type] -> SetPulseTbCuts(fPulseCuts[type][5],fPulseCuts[type][6]);
        ana[type] -> SetCvsGroup(fCvsDX,fCvsDY,fCvsGroupDX,fCvsGroupDY);
        ana[type] -> SetInvertChannel(fPulseCuts[type][8]);
        ana[type] -> SetFixPedestal(fPulseCuts[type][9]);
    }

    auto numEvents = run -> GetNumEvents();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        run -> GetEntry(iEvent);

        int numChannels = channelArray -> GetEntries();
        for (int iChannel=0; iChannel<numChannels; iChannel++)
        {
            auto channel = (MMChannel *) channelArray -> At(iChannel);
            auto cobo = channel -> GetCobo();
            auto asad = channel -> GetAsad();
            auto aget = channel -> GetAget();
            auto chan = channel -> GetChan();
            auto data = channel -> GetWaveformY();
            auto type = GetType(cobo,asad,aget,chan);
            auto caac = iEvent*1000000 + cobo*100000 + asad*10000 + aget*1000 + chan;

            bool selected = false;
            for (auto type0 : fSelTypes)
                if (type==type0)
                    selected = true;
            if (selected==false)
                continue;

            ana[type] -> AddChannel(data, iEvent, cobo, asad, aget, chan);

            if (fDrawChannel)
            {
                bool addCvsAll = (ana[type]->GetNumHistChannels()<fCvsGroupDX*fCvsGroupDY*fNumCvsGroup);
                bool addCvsGood = (ana[type]->IsCollected() && ana[type]->GetNumGoodChannels()<=fCvsGroupDX*fCvsGroupDY*fNumCvsGroup);
                bool addCvs = (fDrawOnlyGoodChannels?addCvsGood:addCvsAll);
                if (addCvs)
                {
                    bool cvsIsNew = ana[type] -> DrawChannel();
                    auto cvs = ana[type] -> GetGroupCanvas();
                    if (cvsIsNew) {
                        e_add(cvs,fTypeNames[type]);
                    }
                }
            }

            if (fCAAC>0 && caac==fCAAC) {
                ana[type] -> DumpChannel();
            }
        }
    }

    if (fCAAC>0)
        return;

    int divX = 1;
    int divY = 1;
    int divAX = 1;
    int divAY = 2;
    int fNumSel = 0; for (auto type : fSelTypes) fNumSel++;
    if (fNumSel>1)  { divX = 1; divY = 2; }
    if (fNumSel>2)  { divX = 2; divY = 2; }
    if (fNumSel>4)  { divX = 3; divY = 2; }
    if (fNumSel>6)  { divX = 3; divY = 3; }
    if (fNumSel>9)  { divX = 4; divY = 3; }
    if (fNumSel>12) { divX = 5; divY = 3; }
    if (fNumSel>=1)  { divAX = 1; divAY = 2; }
    if (fNumSel>=2)  { divAX = 2; divAY = 2; }
    if (fNumSel>=4)  { divAX = 3; divAY = 2; }
    if (fNumSel>=6)  { divAX = 3; divAY = 3; }
    if (fNumSel>=9)  { divAX = 4; divAY = 3; }
    if (fNumSel>=12) { divAX = 5; divAY = 3; }

    if (fDrawAverage) {
        auto cvsAverageAll = new TCanvas("cvsAverageAll","cvsAverageAll",fCvsDX,fCvsDY);
        cvsAverageAll -> Divide(divAX,divAY);
        e_add(cvsAverageAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawAverage(cvsAverageAll->cd(iType+1));
        }

        cvsAverageAll -> cd(divAX*divAY);
        auto frame = new TH2D("frameWidthAll",";y;width",100,0,100,100,0,60);
        frame -> SetStats(0);
        frame -> Draw();
        auto legend = new TLegend(0.5,0.5,0.85,0.95);
        for (auto iType=0; iType<fNumSel; ++iType)
        {
            auto type = fSelTypes[iType];
            auto graph = new TGraph();
            graph -> SetMarkerStyle(20+type);
            graph -> SetMarkerSize(0.6);
            graph -> SetMarkerColor(type+1);
            graph -> SetLineColor  (type+1);
            if (type==4) {
                graph -> SetMarkerColor(1);
                graph -> SetLineColor  (1);
                graph -> SetLineStyle(2);
            }
            if (type>8) {
                graph -> SetMarkerColor(type+1-8);
                graph -> SetLineColor  (type+1-8);
                graph -> SetLineStyle(2);
            }
            for (auto ratio : {0.05, 0.25, 0.50, 0.75})
            {
                double x0,x1,error;
                auto width = ana[type]->FullWidthRatioMaximum(ana[type]->GetHistAverage(),ratio,4,x0,x1,error);
                graph -> SetPoint(graph->GetN(),ratio*100,width);
            }
            graph -> Draw("samelp");
            legend -> AddEntry(graph,fTypeNames[type],"lp");
        }
        legend -> Draw();
    }

    if (fDrawAccumulate) {
        auto cvsAccumulateAll = new TCanvas("cvsAccumulateAll","cvsAccumulateAll",fCvsDX,fCvsDY);
        cvsAccumulateAll -> Divide(divX,divY);
        e_add(cvsAccumulateAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawAccumulate(cvsAccumulateAll->cd(iType+1));
        }
    }

    if (fDrawWidth) {
        auto cvsWidthAll = new TCanvas("cvsWidthAll","cvsWidthAll",fCvsDX,fCvsDY);
        cvsWidthAll -> Divide(divX,divY);
        e_add(cvsWidthAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawWidth(cvsWidthAll->cd(iType+1));
        }
    }

    if (fDrawHeight) {
        auto cvsHeightAll = new TCanvas("cvsHeightAll","cvsHeightAll",fCvsDX,fCvsDY);
        cvsHeightAll -> Divide(divX,divY);
        e_add(cvsHeightAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawHeight(cvsHeightAll->cd(iType+1));
        }
    }

    if (fDrawPulseTb) {
        auto cvsPulseTbAll = new TCanvas("cvsPulseTbAll","cvsPulseTbAll",fCvsDX,fCvsDY);
        cvsPulseTbAll -> Divide(divX,divY);
        e_add(cvsPulseTbAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawPulseTb(cvsPulseTbAll->cd(iType+1));
        }
    }

    if (fDrawPedestal) {
        auto cvsPedestalAll = new TCanvas("cvsPedestalAll","cvsPedestalAll",fCvsDX,fCvsDY);
        cvsPedestalAll -> Divide(divX,divY);
        e_add(cvsPedestalAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawPedestal(cvsPedestalAll->cd(iType+1));
        }
    }

    if (fDrawResidual) {
        auto cvsResidualAll = new TCanvas("cvsResidualAll","cvsResidualAll",fCvsDX,fCvsDY);
        cvsResidualAll -> Divide(divX,divY);
        e_add(cvsResidualAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawResidual(cvsResidualAll->cd(iType+1));
        }
    }

    if (fDrawReference) {
        auto cvsReferenceAll = new TCanvas("cvsReferenceAll","cvsReferenceAll",fCvsDX,fCvsDY);
        cvsReferenceAll -> Divide(divX,divY);
        e_add(cvsReferenceAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawReference(cvsReferenceAll->cd(iType+1));
        }
    }

    if (fDrawHeightWidth) {
        auto cvsHeightWidthAll = new TCanvas("cvsHeightWidthAll","cvsHeightWidthAll",fCvsDX,fCvsDY);
        cvsHeightWidthAll -> Divide(divX,divY);
        e_add(cvsHeightWidthAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType) {
            auto type = fSelTypes[iType];
            ana[type] -> DrawHeightWidth(cvsHeightWidthAll->cd(iType+1));
        }
    }

    if (fDrawProjY)
    {
        auto cvsAverageEAll = new TCanvas("cvsAverageEAll","cvsAverageEAll",fCvsDX,fCvsDY);
        cvsAverageEAll -> Divide(divX,divY);
        e_add(cvsAverageEAll,"Summary");
        for (auto iType=0; iType<fNumSel; ++iType)
        {
            auto type = fSelTypes[iType];
            cvsAverageEAll -> cd(iType+1);
            //auto histReference2 = new TH2D(Form("histReference2_%s",fTypeNames[iType]),";tb;y",100,0,350,140,-0.2,1.2);
            //histReference2 -> SetStats(0);
            //histReference2 -> Draw();
            ana[type] -> GetReferencePulse() -> Draw("aplz");
        }
    }

    if (fDrawPulseParameter)
    {
        //auto histArray = ana[type] -> GetHistArray();
        for (auto type : fSelTypes)
        {
            ana[type] -> WriteReferencePulse(20,40);
            //ana[type] -> WriteReferencePulse(0,0);
            //ana[type] -> GetTree();
            auto cvsPulseParameter = new TCanvas(Form("cvsPulseParameter_%s",fTypeNames[type]),"cvsAverageEAll",fCvsDX,fCvsDY);
            cvsPulseParameter -> Divide(2,2);
            e_add(cvsPulseParameter,"Debug");
            cvsPulseParameter -> cd(1); ana[type] -> GetHistPedestal() -> Draw();
            cvsPulseParameter -> cd(2); ana[type] -> GetHistPedestalResidual() -> Draw();
            cvsPulseParameter -> cd(3); ana[type] -> GetHistYFluctuation() -> Draw();
            e_info << fTypeNames[type]
                << " bgl " << ana[type] -> GetBackGroundLevel()
                << " bge " << ana[type] -> GetBackGroundError()
                << " fll " << ana[type] -> GetFluctuationLevel()  << endl;
        }
    }

    if (!fDrawPulseParameter && fWriteReferencePulse) {
        for (auto type : fSelTypes)
            ana[type] -> WriteReferencePulse(20,40);
            //ana[type] -> WriteReferencePulse(0,0);
    }

    if (fWriteTree) {
        for (auto type : fSelTypes)
            ana[type] -> WriteTree();
    }

    if (fSaveAll) e_save_all();
}

void SetAllFalse()
{
    SetTypeAll();
    fDataPath = "data";
    fMaxInputFiles = 1;
    fNumCvsGroup = 1;
    fCAAC = 0;
    fDrawChannel = false;
    fDrawOnlyGoodChannels = false;
    fDrawAverage = false;
    fDrawAccumulate = false;
    fDrawProjY = false;
    fDrawHeightWidth = false;
    fDrawWidth = false;
    fDrawHeight = false;
    fDrawPulseTb = false;
    fDrawPedestal = false;
    fDrawResidual = false;
    fDrawReference = false;
    fWriteReferencePulse = false;
    fWriteTree = false;
    fSaveAll = false;
}

void SetRunAll() 
{
    SetAllFalse();
    fDataPath = "data100";
    fMaxInputFiles = 100;
    fWriteReferencePulse = true;
    fWriteTree = true;
}

void SetDrawCheckAllResults()
{
    //SetAllFalse();
    //SetMMType();
    e_batch();
    fMaxInputFiles = 5;
    fNumCvsGroup = 4;
    fCAAC = 0;
    fDrawChannel = true;
    fDrawAverage = true;
    fDrawAccumulate = true;
    fDrawProjY = true;
    fDrawHeightWidth = true;
    fDrawWidth = true;
    fDrawHeight = true;
    fDrawPulseTb = true;
    fDrawPedestal = true;
    fDrawResidual = true;
    fDrawReference = true;
    fSaveAll = true;
}

void SetDrawCheckGoodChannels()
{
    SetAllFalse();
    e_batch();
    SetType(0);
    fNumCvsGroup = 5;
    fDrawChannel = true;
    fDrawOnlyGoodChannels = true;
    fDrawAverage = true;
    fDrawAccumulate = true;
    fDrawProjY = true;
    fDrawHeightWidth = true;
    fDrawWidth = true;
    fDrawHeight = true;
    fDrawPulseTb = true;
    fDrawPedestal = true;
    fDrawResidual = true;
    fDrawReference = true;
}

void SetExtractBuffer(int ecaac)
{
    SetAllFalse();
    e_batch();
    SetType(0);
    fMaxInputFiles = 1;
    fNumCvsGroup = 10;
    fCAAC = ecaac;
    fDrawChannel = true;
    fDrawOnlyGoodChannels = true;
}

void SetDrawResidual()
{
    SetAllFalse();
    e_batch();
    SetType();
    fMaxInputFiles = 1;
    fNumCvsGroup = 1;
    fDrawAverage = true;
    fDrawAccumulate = true;
    fDrawProjY = true;
    fDrawResidual = true;
}

void SetDebugPulseData(int i)
{
    SetAllFalse();
    e_batch();
    SetType(i);
    fMaxInputFiles = 5;
    fNumCvsGroup = 10;
    fDrawChannel = true;
    //fDrawAverage = true;
    //fDrawAccumulate = true;
    fDrawPulseParameter = true;
    fWriteReferencePulse = true;
}
