#include "LTPadPlane.h"
#include "TPaveText.h"
#include "LAMPSTPC.h"
#include "LKEventHeader.h"
#include "LKWindowManager.h"

ClassImp(LTPadPlane);

LTPadPlane* LTPadPlane::fInstance = nullptr;
LTPadPlane* LTPadPlane::GetPlane() { return fInstance; }

LTPadPlane::LTPadPlane()
{
    fInstance = this;
    fName = "LTPadPlane";
    fAxis1 = LKVector3::kX;
    fAxis2 = LKVector3::kZ;
}

bool LTPadPlane::Init()
{
    e_info << "Initializing LTPadPlane" << std::endl;

    fHistChannel1 = new TH1D("hist_channel_1","channel buffer 1;time-bucket;charge",512,0,512);
    fHistChannel2 = new TH1D("hist_channel_2","channel buffer 2;time-bucket;charge",512,0,512);
    for (auto histChannel : {fHistChannel1,fHistChannel2}) {
        histChannel -> SetStats(0);
        histChannel -> GetXaxis() -> SetLabelSize(0.065);
        histChannel -> GetYaxis() -> SetLabelSize(0.065);
        histChannel -> GetXaxis() -> SetTitleSize(0.065);
        histChannel -> GetYaxis() -> SetTitleSize(0.065);
        histChannel -> GetXaxis() -> SetTitleOffset(1.20);
        histChannel -> GetYaxis() -> SetTitleOffset(0.80);
    }

    TString mapFileName = fPar -> GetParString("LTPadPlane/position_map");
    if (mapFileName.IsNull()) {
        lk_error << mapFileName << " is null!" << endl;
        return false;
    }

    fHistPadPlaneFrame = new TH2Poly("frameLTPP","LAMPS TPC Pad Plane;x (mm);y (mm)",-750,+750,-750,+750);
    fHistSideViewFrame = new TH2Poly("frameLTSV","LAMPS TPC Side View;y (mm);z (mm)",-350,1350,-750,+750);
    fHistPadPlaneFrame -> SetStats(0);
    fHistSideViewFrame -> SetStats(0);
    fHistPadPlane = new TH2Poly("histLTPP","LAMPS TPC Pad Plane;x (mm);y (mm)",-750,+750,-750,+750);
    fHistSideView = new TH2Poly("histLTSV","LAMPS TPC Side View;z (mm);y (mm)",-350,1350,-750,+750);
    fHistPadPlane -> SetStats(0);
    fHistSideView -> SetStats(0);
    ifstream filePositionMap(mapFileName);

    int numChannels = 2698*8;
    int cobo, asad, aget, ch, section, npad, nlayer, nasad;
    double x, y;
    for(int i=0; i<numChannels; i++)
    {
        filePositionMap >> cobo >> asad >> aget >> ch >> x >> y >> section >> npad >> nlayer >> nasad;

        double x1 = x-2;
        double x2 = x+2;
        double y1 = y-2;
        double y2 = y+2;
        auto bin = fHistPadPlaneFrame -> AddBin(x1,y1,x2,y2);
        fHistPadPlane -> AddBin(x1,y1,x2,y2);
        int caac = cobo*10000 + asad*1000 + aget*100 + ch;
        fMapCAACToBin.insert(std::pair<int, int>(caac,bin));
        fMapBinToCAAC.insert(std::pair<int, int>(bin,caac));
        fMapBinToX1.insert(std::pair<int, double>(bin,x1));
        fMapBinToX2.insert(std::pair<int, double>(bin,x2));
        fMapBinToY1.insert(std::pair<int, double>(bin,y1));
        fMapBinToY2.insert(std::pair<int, double>(bin,y2));
    }
    filePositionMap.close();

    SetDataFromBranch();

    fGraphBoundaryPadPlane = new TGraph();
    fGraphBoundaryPadPlane -> SetName("fGraphBoundaryPadPlane");
    fGraphBoundaryPadPlane -> SetLineColor(kRed);
    fGraphBoundaryPadPlane -> SetLineWidth(2);

    fGraphBoundarySideView = new TGraph();
    fGraphBoundarySideView -> SetName("fGraphBoundarySideView");
    fGraphBoundarySideView -> SetLineColor(kRed);
    fGraphBoundarySideView -> SetLineWidth(2);

    return true;
}

void LTPadPlane::Clear(Option_t *option)
{
    LKPadPlane::Clear(option);
}

void LTPadPlane::Print(Option_t *option) const
{
    e_info << "LTPadPlane" << std::endl;
}

bool LTPadPlane::IsInBoundary(Double_t x, Double_t y)
{
    return true;
}

Int_t LTPadPlane::FindChannelID(Double_t x, Double_t y)
{
    //return fHistPadPlane -> FindBin(x,y);
    return -1;
}

Int_t LTPadPlane::FindChannelID(Int_t section, Int_t row, Int_t layer)
{
    // int id = 10000*section + 100*row + layer;
    // return id;
    return -1;
}

TCanvas* LTPadPlane::GetCanvas(Option_t *option)
{
    if (fCanvas==nullptr) {
        fCanvas = LKWindowManager::GetWindowManager() -> CanvasResize("LTPadPlane",1100,700,0.9);
        auto pad1 = new TPad("pad1","",0,230./700,0.5,1);
        pad1 -> SetMargin(0.12,0.15,0.1,0.1);
        pad1 -> SetNumber(1);
        pad1 -> AddExec("ex", "LTPadPlane::MouseClickEvent1()");
        pad1 -> Draw();
        auto pad2 = new TPad("pad2","",0,0,0.5,230./700);
        pad2 -> SetMargin(0.12,0.05,0.20,0.12);
        pad2 -> SetNumber(2);
        pad2 -> Draw();
        auto pad3 = new TPad("pad1","",0.5,230./700,1,1);
        pad3 -> SetMargin(0.12,0.15,0.1,0.1);
        pad3 -> SetNumber(3);
        pad3 -> AddExec("ex", "LTPadPlane::MouseClickEvent2()");
        pad3 -> Draw();
        auto pad4 = new TPad("pad2","",0.5,0,1,230./700);
        pad4 -> SetMargin(0.12,0.05,0.20,0.12);
        pad4 -> SetNumber(4);
        pad4 -> Draw();
        fCanvas -> Modified();
        fCanvas -> Update();
    }
    return fCanvas;
}

TH2* LTPadPlane::GetHist(Option_t *option)
{
    return (TH2Poly *) fHistPadPlane;
    //return (TH2Poly *) fHistSideView;
}

bool LTPadPlane::SetDataFromBranch()
{
    if (fRun==nullptr)
        return false;

    fBufferArray = fRun -> GetBranchA("RawData");
    fHitArray = fRun -> GetBranchA("Hit");
    fEventHeaderHolder = fRun -> GetBranchA("EventHeader");
    return true;
}

void LTPadPlane::FillDataToHist()
{
    if (fHitArray==nullptr)
        return;

    fHistPadPlane -> Reset("ICES");
    fHistSideView -> Reset("ICES");

    Long64_t eventIDFromRun = -1;
    int evenIDFromHeader = -1;
    if (fRun!=nullptr) eventIDFromRun = fRun -> GetCurrentEventID();
    if (fEventHeaderHolder!=nullptr) {
        auto eventHeader = (LKEventHeader*) fEventHeaderHolder -> At(0);
        if (eventHeader!=nullptr) {
            evenIDFromHeader = eventHeader -> GetEventNumber();
        }
    }
    lk_info << "LKRun eventID: " << eventIDFromRun << endl;
    lk_info << "Experiment eventID: " << evenIDFromHeader << endl;
    if (eventIDFromRun>=0 && evenIDFromHeader>=0) {
        fHistPadPlaneFrame -> SetTitle(Form("LAMPS TPC Pad Plane, Event %lld (%d)",eventIDFromRun,evenIDFromHeader));
        fHistSideViewFrame -> SetTitle(Form("LAMPS TPC Pad Plane, Event %lld (%d)",eventIDFromRun,evenIDFromHeader));
    }
    else if (eventIDFromRun>=0 && evenIDFromHeader<0) {
        fHistPadPlaneFrame -> SetTitle(Form("LAMPS TPC Pad Plane, Event %lld",eventIDFromRun));
        fHistSideViewFrame -> SetTitle(Form("LAMPS TPC Pad Plane, Event %lld",eventIDFromRun));
    }

    lk_info << "# of hits : " << fHitArray -> GetEntriesFast() << endl;

    for (auto hitArray : {fHitArray})
    {
        auto numHits = hitArray -> GetEntries();
        for (auto iHit=0; iHit<numHits; ++iHit)
        {
            auto hit = (LKHit *) hitArray -> At(iHit);
            auto position = hit -> GetPosition();
            fHistPadPlane -> Fill(position.X(), position.Y(), hit->GetCharge());
            fHistSideView -> Fill(position.Z(), position.Y(), hit->GetCharge());
        }
    }

    SelectAndDrawChannel1();
    SelectAndDrawChannel2();
}


void LTPadPlane::SelectAndDrawChannel(bool isChannel1, Int_t bin)
{
    if (isChannel1==false)
        return;

    auto cvsPlane = fCanvas -> cd(3);
    auto cvsChannel = fCanvas -> cd(4);
    if (isChannel1) {
        cvsPlane = fCanvas -> cd(1);
        cvsChannel = fCanvas -> cd(2);
    }

    bool existHitArray = false;
    bool existBufferArray = (fBufferArray!=nullptr);
    if (isChannel1) existHitArray = (fHitArray!=nullptr);
    else         existHitArray = (fHitArray!=nullptr);

    fHitArrayList.clear();
    fHitArrayList.push_back(fHitArray);
    if (bin<0) {
        if (existHitArray) {
            LKHit *hit = nullptr;
            for (auto hitArray : fHitArrayList) {
                auto numHits = hitArray -> GetEntries();
                for (auto iHit=0; iHit<numHits; ++iHit)
                    hit = (LKHit *) hitArray -> At(iHit);
            }
            if (hit==nullptr)
                return;
            auto caac = hit -> GetChannelID();
            if (isChannel1) bin = fMapCAACToBin[caac];
            if (isChannel1) lk_info << "SelectAndDrawChannel CAAC=" << Form("%05d",bin) << endl;
        }
    }

    if (bin<0) {
        lk_warning << "bin<0" << bin << endl;
        return;
    }

    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
    TGraph* graphBoundary = fGraphBoundarySideView;
    if (isChannel1) graphBoundary = fGraphBoundaryPadPlane;

    x1 = fMapBinToX1[bin];
    x2 = fMapBinToX2[bin];
    y1 = fMapBinToY1[bin];
    y2 = fMapBinToY2[bin];
    z1 = -350;
    z2 = 1350;
    x0 = (x1 + x2)/2.;
    y0 = (y1 + y2)/2.;
    z0 = (z1 + z2)/2.;
    if (isChannel1) {
        graphBoundary -> Set(0);
        graphBoundary -> SetPoint(0,x1,y1);
        graphBoundary -> SetPoint(1,x2,y1);
        graphBoundary -> SetPoint(2,x2,y2);
        graphBoundary -> SetPoint(3,x1,y2);
        graphBoundary -> SetPoint(4,x1,y1);
        graphBoundary -> SetLineColor(kBlue);
    }
    else {
        graphBoundary -> Set(0);
        graphBoundary -> SetPoint(0,z1,y1);
        graphBoundary -> SetPoint(1,z2,y1);
        graphBoundary -> SetPoint(2,z2,y2);
        graphBoundary -> SetPoint(3,z1,y2);
        graphBoundary -> SetPoint(4,z1,y1);
        graphBoundary -> SetLineColor(kBlue);
    }
    cvsPlane -> cd();
    graphBoundary -> Draw("samel");

    ////////////////////////////////////////////////////////////////////

    fCAAC = fMapBinToCAAC[bin];

    auto texat = (LAMPSTPC*) fDetector;
    
    TH1D* histChannel = nullptr;
    if (isChannel1) histChannel = fHistChannel1;
    else         histChannel = fHistChannel2;
    histChannel -> Reset();
    fCurrElectronicsID = texat -> GetElectronicsID(fCAAC);
    histChannel -> SetTitle(Form("CAAC=%d, eID=%d, position=(%.2f, %.2f), #Hits=%d",fCAAC,fCurrElectronicsID,x0,z0,0));
    
    cvsChannel -> Modified();
    cvsChannel -> Update();

    if (!existBufferArray)
        return;

    fCurrSelChannel = nullptr;
    auto numChannels = fBufferArray -> GetEntries();
    for (auto iChannel=0; iChannel<numChannels; ++iChannel)
    {
        auto channel0 = (GETChannel* ) fBufferArray -> At(iChannel);
        if (fCAAC==channel0->GetCAAC()) {
            fCurrSelChannel = channel0;
            break;
        }
    }

    if (fCurrSelChannel==nullptr) {
        cvsChannel -> cd();
        histChannel -> Draw();
        return;
    }

    auto buffer = fCurrSelChannel -> GetWaveformY();
    for (auto tb=0; tb<360; ++tb)
        histChannel -> SetBinContent(tb+1,buffer[tb]);
    cvsChannel -> cd();
    histChannel -> Draw();
    //if (texat==nullptr) ...;
    if (existHitArray)
    {
        auto chAna = texat -> GetChannelAnalyzer(fCurrElectronicsID);
        int countHitsInChannel = 0;
        for (auto hitArray : fHitArrayList) {
            auto numHits = hitArray -> GetEntries();
            for (auto iHit=0; iHit<numHits; ++iHit)
            {
                auto hit = (LKHit *) hitArray -> At(iHit);
                if (fCAAC==hit -> GetChannelID())
                {
                    ++countHitsInChannel;
                    if (texat!=nullptr)
                    {
                        auto pulse = chAna -> GetPulse();
                        auto graph = pulse -> GetPulseGraph(hit->GetY(), hit->GetCharge(), hit->GetPedestal());
                        cvsChannel -> cd();
                        graph -> Draw("samelx");
                    }
                    graphBoundary -> SetLineColor(kRed);
                    //graphBoundary -> Draw("samel");
                }
            }
            //histChannel -> Draw();
            histChannel -> SetTitle(Form("CAAC=%d, eID=%d, position=(%.2f, %.2f), #Hits=%d",fCAAC,fCurrElectronicsID,x0,z0,countHitsInChannel));
        }
    }

    fCanvas -> Modified();
    fCanvas -> Update();
}

void LTPadPlane::WriteCurrentChannel(TString name)
{
    if (name.IsNull()) name = "selected_channel.root";
    lk_info << "Creating " << name << endl;
    auto fileTestChannel = new TFile(name,"recreate");
    fCurrSelChannel -> Write("channel",TObject::kSingleKey);
    for (auto hitArray : fHitArrayList) {
        int countHitsInChannel = 0;
        auto numHits = hitArray -> GetEntries();
        for (auto iHit=0; iHit<numHits; ++iHit) {
            auto hit = (LKHit *) hitArray -> At(iHit);
            if (fCAAC==hit -> GetChannelID()) {
                fileTestChannel -> cd();
                hit -> Write(Form("hit%d",countHitsInChannel++),TObject::kSingleKey);
            }
        }
    }
    if (fPar!=nullptr) {
        auto run = LKRun::GetRun();
        if (run!=nullptr)
            fPar -> AddPar("LTPadPlane/WriteCurrentChannel/eventID",LKRun::GetRun()->GetCurrentEventID());
        fPar -> AddPar("LTPadPlane/WriteCurrentChannel/electronicsID",fCurrElectronicsID);
        fPar -> Write(fPar->GetName(),TObject::kSingleKey);
    }
    if (fEventHeaderHolder!=nullptr) {
        auto eventHeader = (LKEventHeader*) fEventHeaderHolder -> At(0);
        if (eventHeader!=nullptr) {
            eventHeader -> Write("EventHeader",TObject::kSingleKey);
        }
    }
}

void LTPadPlane::DrawFrame(Option_t *option)
{
    ;
}

void LTPadPlane::Draw(Option_t *option)
{
    //SetDataFromBranch();
    FillDataToHist();

    auto cvs = GetCanvas();

    cvs -> cd(1);
    if (fHistPadPlane->GetEntries()==0) {
        fHistPadPlaneFrame -> Draw();
    }
    else {
        fHistPadPlaneFrame -> Draw();
        fHistPadPlane -> SetMinimum(0.1);
        fHistPadPlane -> Draw("same colz");
    }

    cvs -> cd(2);
    fHistChannel1 -> Draw();
    cvs -> cd(2) -> Modified();
    cvs -> cd(2) -> Update();
    auto ttt1 = (TPaveText*) (cvs->cd(2)->GetListOfPrimitives()) -> FindObject("title");
    ttt1 -> SetTextSize(0.065);
    ttt1 -> SetTextAlign(12);
    cvs -> cd(2) -> Modified();
    cvs -> cd(2) -> Update();

    cvs -> cd(3);
    if (fHistSideView->GetEntries()==0)
        fHistSideViewFrame -> Draw();
    else {
        fHistSideViewFrame -> Draw();
        fHistSideView -> SetMinimum(0.1);
        fHistSideView -> Draw("same colz");
    }

    cvs -> cd(4);
    fHistChannel2 -> Draw();
    cvs -> cd(4) -> Modified();
    cvs -> cd(4) -> Update();
    auto ttt2 = (TPaveText*) (cvs->cd(4)->GetListOfPrimitives()) -> FindObject("title");
    ttt2 -> SetTextSize(0.065);
    ttt2 -> SetTextAlign(32);
    cvs -> cd(4) -> Modified();
    cvs -> cd(4) -> Update();
}

void LTPadPlane::MouseClickEvent1()
{
    if (gPad==nullptr)
        return;

    TObject* select = gPad -> GetCanvas() -> GetClickSelected();
    if (select == nullptr)
        return;

    bool isNotH2 = !(select -> InheritsFrom(TH1::Class()));
    //bool isNotGraph = !(select -> InheritsFrom(TGraph::Class()));
    //if (isNotH2 && isNotGraph)
    if (isNotH2)
        return;

    int xEvent = gPad -> GetEventX();
    int yEvent = gPad -> GetEventY();
    int xAbs = gPad -> AbsPixeltoX(xEvent);
    int yAbs = gPad -> AbsPixeltoY(yEvent);
    double xOnClick = gPad -> PadtoX(xAbs);
    double yOnClick = gPad -> PadtoY(yAbs);

    TH2* hist = dynamic_cast<TH2*> (select);
    int binCurr = hist -> FindBin(xOnClick, yOnClick);
    int binLast = gPad -> GetUniqueID();
    if (binCurr==binLast)
        return;

    gPad -> SetUniqueID(binCurr);
    gPad -> GetCanvas() -> SetClickSelected(nullptr);

    if (binCurr<=0)
        return;

    LTPadPlane::GetPlane() -> SelectAndDrawChannel1(binCurr);
}

void LTPadPlane::MouseClickEvent2()
{
    if (gPad==nullptr)
        return;

    TObject* select = gPad -> GetCanvas() -> GetClickSelected();
    if (select == nullptr) {
        gPad -> GetCanvas() -> SetClickSelected(nullptr);
        return;
    }

    bool isNotH2 = !(select -> InheritsFrom(TH1::Class()));
    //bool isNotGraph = !(select -> InheritsFrom(TGraph::Class()));
    //if (isNotH2 && isNotGraph)
    if (isNotH2) {
        gPad -> GetCanvas() -> SetClickSelected(nullptr);
        return;
    }

    int xEvent = gPad -> GetEventX();
    int yEvent = gPad -> GetEventY();
    int xAbs = gPad -> AbsPixeltoX(xEvent);
    int yAbs = gPad -> AbsPixeltoY(yEvent);
    double xOnClick = gPad -> PadtoX(xAbs);
    double yOnClick = gPad -> PadtoY(yAbs);

    TH2* hist = dynamic_cast<TH2*> (select);
    int binCurr = hist -> FindBin(xOnClick, yOnClick);
    int binLast = gPad -> GetUniqueID();
    if (binCurr==binLast || binCurr<=0) {
        gPad -> GetCanvas() -> SetClickSelected(nullptr);
        return;
    }

    gPad -> SetUniqueID(binCurr);
    gPad -> GetCanvas() -> SetClickSelected(nullptr);
    LTPadPlane::GetPlane() -> SelectAndDrawChannel2(binCurr);
}
