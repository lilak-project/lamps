#include "LTPadPlane.h"

#include "TGraph.h"
#include "TMath.h"
#include "TVector2.h"
#include "TH2Poly.h"
#include "TCollection.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPaveText.h"

#include "LKWindowManager.h"
#include "GETChannel.h"

#include <iostream>
using namespace std;

ClassImp(LTPadPlane)

LTPadPlane* LTPadPlane::fInstance = nullptr;
LTPadPlane* LTPadPlane::GetPlane() { return fInstance; }

LTPadPlane::LTPadPlane()
    :LKPadPlane("LTPadPlane", "pad plane with rectangular pads following the circle line for LAMPS-TPC")
{
    fInstance = this;
    fName = "LTPadPlane";
    fAxis1 = LKVector3::kX;
    fAxis2 = LKVector3::kY;
}

bool LTPadPlane::Init()
{
    if (fPar -> CheckPar("LTPadPlane/rMinTPC"))      fRMin         = fPar -> GetParDouble("LTPadPlane/rMinTPC");
    if (fPar -> CheckPar("LTPadPlane/rMaxTPC"))      fRMax         = fPar -> GetParDouble("LTPadPlane/rMaxTPC");
    if (fPar -> CheckPar("LTPadPlane/padGap"))       fPadGap       = fPar -> GetParDouble("LTPadPlane/padGap");
    if (fPar -> CheckPar("LTPadPlane/rTopCut"))      fRTopCut      = fPar -> GetParDouble("LTPadPlane/rTopCut");
    if (fPar -> CheckPar("LTPadPlane/padWidth"))     fPadWidth     = fPar -> GetParDouble("LTPadPlane/padWidth");
    if (fPar -> CheckPar("LTPadPlane/padHeight"))    fPadHeight    = fPar -> GetParDouble("LTPadPlane/padHeight");
    if (fPar -> CheckPar("LTPadPlane/radiusLayer0")) fRadiusLayer0 = fPar -> GetParDouble("LTPadPlane/radiusLayer0");
    if (fPar -> CheckPar("LTPadPlane/numLayers"))    fNumLayers    = fPar -> GetParInt   ("LTPadPlane/numLayers");

    fXSpacing = fPadGap + fPadWidth;
    fRSpacing = fPadGap + fPadHeight;

    fTanPi1o8 = TMath::Tan(TMath::Pi()*1./8.);
    fTanPi3o8 = TMath::Tan(TMath::Pi()*3./8.);
    fTanPi5o8 = TMath::Tan(TMath::Pi()*5./8.);
    fTanPi7o8 = TMath::Tan(TMath::Pi()*7./8.);
    for (Int_t i = 0; i < 8; i++) {
        fCosPiNo4[i] = TMath::Cos(TMath::Pi()*i/4.);
        fSinPiNo4[i] = TMath::Sin(TMath::Pi()*i/4.);
    }

    for (auto section=0; section<8; ++section) {

        Double_t x1 = -208;
        Double_t x2 = 208;
        Double_t x3 = -38;
        Double_t x4 = 38;
        Double_t y1 = 100;
        Double_t y2 = 510;
        TVector2 sectionCorner1(x1, y2);
        TVector2 sectionCorner2(x3, y1);
        TVector2 sectionCorner3(x4, y1);
        TVector2 sectionCorner4(x2, y2);

        Double_t phiSection = section * TMath::Pi()/4.;
        fPosSectionCorner[section][0] = sectionCorner1.Rotate(phiSection);
        fPosSectionCorner[section][1] = sectionCorner2.Rotate(phiSection);
        fPosSectionCorner[section][2] = sectionCorner3.Rotate(phiSection);
        fPosSectionCorner[section][3] = sectionCorner4.Rotate(phiSection);
    }

    Double_t xCorner[5] = {0};
    Double_t yCorner[5] = {0};

    for (int layer=0; layer<fNumLayers; ++layer)
    {
        double radius = fRadiusLayer0 + layer*fRSpacing;
        double radius1 = fRadiusLayer0 + layer*fRSpacing - 0.5*fPadHeight;
        double radius2 = fRadiusLayer0 + layer*fRSpacing + 0.5*fPadHeight;
        int isInnerPad = true;

        LKPad *padAtCurrentRow[8][2]; // section, R/L
        LKPad *padAtPreviousRow[8][2]; // section, R/L

        int firstRow = 0;
        int row = 0;
        for (int iRow=1; isInnerPad; ++iRow)
        {
            double xPad = (iRow -.5) * fXSpacing;
            double xPad1 = (iRow -.5) * fXSpacing - 0.5*fPadWidth;
            double xPad2 = (iRow -.5) * fXSpacing + 0.5*fPadWidth;
            double yPad = sqrt(radius*radius - xPad*xPad);

            bool continueRow = false;
            if (layer==41&&iRow<55) continueRow = true;
            if (layer==40&&iRow<46) continueRow = true;
            if (layer==39&&iRow<35) continueRow = true;
            if (layer==38&&iRow<17) continueRow = true;
            //if (yPad-fPadHeight*0.3 > fRTopCut) continueRow = true;

            if (continueRow) {
                fNumSkippedHalfRows[layer]++;
                continue;
            }
            else {
                if (firstRow==0) {
                    firstRow = iRow;
                }
                row++;
            }

            if (yPad < GetPadCenterYBoundAtX(xPad))
            {
                fNumHalfRowsInLayer[layer] = row;
                isInnerPad = false;
            }

            ///////////////////////////////////////////////////////////////////

            int numCorners = 4;
            Double_t yCornerTemp;
            for (auto i : {0,1,2,3})
            {
                int pm = +1;
                double rr = radius;
                if (i==0) { pm = +1; rr = radius2; }
                if (i==1) { pm = -1; rr = radius2; }
                if (i==2) { pm = -1; rr = radius1; }
                if (i==3) { pm = +1; rr = radius1; }
                xCorner[i] = xPad + pm*.5*fPadWidth;
                yCornerTemp = sqrt(rr*rr-xCorner[i]*xCorner[i]);
                if (yCornerTemp>fRTopCut) yCornerTemp=fRTopCut;
                yCorner[i] = yCornerTemp;
            }

            bool sideBoundaryCutWasMade = false;

            if (fDoCutSideBoundary && (yCorner[3] < GetPadCutBoundaryYAtX(xCorner[3]))) {
                sideBoundaryCutWasMade = true;
                if (xCorner[0] > GetPadCutBoundaryXAtR(radius2)) {
                    xCorner[0] = GetPadCutBoundaryXAtR(radius2);
                    yCorner[2] = GetPadCutBoundaryYAtX(xCorner[2]);
                    numCorners = 3;
                }
                else if (yCorner[2] < GetPadCutBoundaryYAtX(xCorner[2])) {
                    yCorner[2] = GetPadCutBoundaryYAtX(xCorner[2]);
                    yCorner[3] = GetPadCutBoundaryYAtX(xCorner[3]);
                }
                else if (yCorner[3] < GetPadCutBoundaryYAtX(xCorner[3])) {
                    xCorner[3] = GetPadCutBoundaryXAtR(radius1);
                    yCorner[3] = fTanPi3o8*xCorner[3];
                    xCorner[4] = xCorner[0];
                    yCorner[4] = GetPadCutBoundaryYAtX(xCorner[4]);
                    numCorners = 5;
                }
            }
            if (fDoCutTopBoundary && yCorner[0] > fRTopCut) {
                yCorner[0] = fRTopCut;
                yCorner[1] = fRTopCut;
            }

            ///////////////////////////////////////////////////////////////////

            for (auto section=0; section<8; ++section)
            {
                Double_t phiSection = section * TMath::Pi()/4.;

                //TVector2 point(xPad, yPad);
                //point = point.Rotate(phiSection);

                TVector2 point;
                for (int iRL : {0,1})
                {
                    int signRL = iRL==0?1:-1;
                    auto pad = NewPad(section, signRL*row, layer);
                    auto padID = pad -> GetPadID();
                    auto slr = iRL + 10*row + 1000*layer + 100000*section;
                    //fMapSLRToPadID.insert(std::pair<int,int>(slr,padID));
                    double xSum = 0;
                    double ySum = 0;
                    for (auto iC=0; iC<numCorners; ++iC) {
                        point = TVector2(signRL*xCorner[iC],yCorner[iC]);
                        point = point.Rotate(phiSection);
                        pad -> AddPadCorner(point.X(), point.Y());
                        xSum += point.X();
                        ySum += point.Y();
                    }
                    pad -> SetPosition(xSum/numCorners, ySum/numCorners);
                    padAtCurrentRow[section][iRL] = pad;
                }
            }

            if (sideBoundaryCutWasMade) {
                SetNeighborPads(padAtCurrentRow[0][1], padAtCurrentRow[1][0]);
                SetNeighborPads(padAtCurrentRow[1][1], padAtCurrentRow[2][0]);
                SetNeighborPads(padAtCurrentRow[2][1], padAtCurrentRow[3][0]);
                SetNeighborPads(padAtCurrentRow[3][1], padAtCurrentRow[4][0]);
                SetNeighborPads(padAtCurrentRow[4][1], padAtCurrentRow[5][0]);
                SetNeighborPads(padAtCurrentRow[5][1], padAtCurrentRow[6][0]);
                SetNeighborPads(padAtCurrentRow[6][1], padAtCurrentRow[7][0]);
                SetNeighborPads(padAtCurrentRow[7][1], padAtCurrentRow[0][0]);
            }

            if (row==1) {
                // row=-1 and row=1 are neighbor
                if (fNumSkippedHalfRows[layer]==0)
                    for (auto section=0; section<8; ++section)
                        SetNeighborPads(padAtCurrentRow[section][0], padAtCurrentRow[section][1]);
            }
            else {
                // [row] and [row-1] are neighbor
                for (auto section=0; section<8; ++section)
                    for (int iRL : {0,1})
                        SetNeighborPads(padAtPreviousRow[section][iRL],padAtCurrentRow[section][iRL]);
            }

            for (auto section=0; section<8; ++section)
                for (int iRL : {0,1})
                    padAtPreviousRow[section][iRL] = padAtCurrentRow[section][iRL];
        }
    }

    for (auto iLayer=fNumLayers-1; iLayer>=0; --iLayer) {
        fNumPadsDownToLayer[iLayer] = 2*fNumHalfRowsInLayer[iLayer] + fNumPadsDownToLayer[iLayer+1];
    }

    fChannelArray -> Sort(); // see LKPad::Compare();

    for (auto i=0; i<fChannelArray->GetEntries(); ++i)
    {
        auto pad = (LKPad *) fChannelArray -> At(i);
        pad -> SetPadID(i);
    }

    for (auto section=0; section<8; ++section) {
        for (auto layer=0; layer<fNumLayers; ++layer) {
            for (auto row=1; row<=fNumHalfRowsInLayer[layer]; ++row)
            {
                auto padL = LKPadPlane::GetPad(section,layer,-row);
                auto padR = LKPadPlane::GetPad(section,layer,+row);

                // pad below
                if (layer>0) {
                    auto rowNb = row + fNumSkippedHalfRows[layer] - fNumSkippedHalfRows[layer-1];
                    if (rowNb<=0||rowNb>fNumHalfRowsInLayer[layer-1]) continue;
                    auto padBelowL = LKPadPlane::GetPad(section,layer-1,-rowNb);
                    auto padBelowR = LKPadPlane::GetPad(section,layer-1,+rowNb);
                    SetNeighborPads(padL,padBelowL);
                    SetNeighborPads(padR,padBelowR);
                }
            }
        }
    }

    fMapCAACToPadID = new int***[22];
    for (int i=0; i<22; ++i) {
        fMapCAACToPadID[i] = new int**[4];
        for (int j=0; j<4; ++j) {
            fMapCAACToPadID[i][j] = new int*[12];
            for (int k=0; k<12; ++k) {
                fMapCAACToPadID[i][j][k] = new int[68];
                for (int l=0; l<68; ++l) {
                    fMapCAACToPadID[i][j][k][l] = -1;
                }
            }
        }
    }


    if (fPar -> CheckPar("LTPadPlane/position_map"))
    {
        double x;
        double y;
        int cobo; // 0-21
        int slot; // 0-3
        int aget; // 0-3
        int channelID; // 0-67
        int section; // 1-8
        int padIDLocal; // 0-2697
        int padID; // global padID create from this class LTPadPlane. This is the order of pad that has been created and stored separately.
        int layer; // 0-41
        int asad; // 1-11
        auto mapFileName = fPar -> GetParString("LTPadPlane/position_map");
        lk_info << "pad electronics mapping: " << mapFileName << endl;
        ifstream fileMap(mapFileName);
        //for (auto i=0; i<21584; ++i) fileMap >> cobo >> slot >> aget >> channelID >> x >> y >> section >> padIDLocal >> layer >> asad >> padID;
        while (fileMap >> cobo >> slot >> aget >> channelID >> x >> y >> section >> padIDLocal >> layer >> asad >> padID)
        {
            if (padID>=0) {
                auto pad = GetPad(padID);
                pad -> SetCoboID(cobo);
                pad -> SetAGETID(aget);
                pad -> SetAsAdID(asad);
                pad -> SetChannelID(channelID);
                auto eleID = GetElectronicsID(cobo,aget,asad,channelID);
                fMapCAACToPadID[cobo][aget][asad][channelID] = padID;
            }
        }
    }

    SetDataFromBranch();

    return true;
}

Int_t LTPadPlane::FindPadID(Int_t section, Int_t layer, Int_t row)
{
    int idLayer = section*fNumPadsDownToLayer[layer];
    if (layer<fNumLayers) idLayer = idLayer + (8-section)*fNumPadsDownToLayer[layer+1];
    int idRow = ((row>0) ? (fNumHalfRowsInLayer[layer] - row) : (fNumHalfRowsInLayer[layer] - row - 1));
    int id = idLayer + idRow;

    return id;
}

Int_t LTPadPlane::FindPadID(Int_t cobo, Int_t aget, Int_t asad, Int_t chan)
{
    //auto eleID = GetElectronicsID(cobo,aget,asad,channelID);
    auto padID = fMapCAACToPadID[cobo][aget][asad][chan];
    return padID;
}

Int_t LTPadPlane::FindPadID(Double_t i, Double_t j)
{
    Int_t section = FindSection(i,j);

    Double_t xRotatedToSec0 =  i*fCosPiNo4[section] + j*fSinPiNo4[section];
    Double_t yRotatedToSec0 = -i*fSinPiNo4[section] + j*fCosPiNo4[section];
    Double_t rRotatedToSec0 =  sqrt(xRotatedToSec0*xRotatedToSec0 + yRotatedToSec0*yRotatedToSec0);

    Double_t rFromSectionBottom = rRotatedToSec0 - fRadiusLayer0 + .5*fPadHeight;
    if (rFromSectionBottom < 0)
        return -1;

    Int_t layer = (Int_t)(rFromSectionBottom/fRSpacing);
    if (layer > fNumLayers)
        return -1;

    Int_t pm = 1;
    if (xRotatedToSec0 < 0) {
        xRotatedToSec0 = -xRotatedToSec0;
        pm = -1;
    }

    if (xRotatedToSec0 < .5*fPadGap)
        return -1;

    Double_t xFromRow0LeftEdge = xRotatedToSec0 - .5*fPadGap;
    Int_t row = (Int_t) (xFromRow0LeftEdge / fXSpacing) + 1;
    if (xFromRow0LeftEdge - (row-1)*fXSpacing > fPadWidth)
        return -1;

    row = row - fNumSkippedHalfRows[layer];

    if (row > fNumHalfRowsInLayer[layer])
        return -1;

    return FindPadID(section,layer,pm*row);
}

LKPad *LTPadPlane::GetPadFromEleID(Int_t cobo, Int_t aget, Int_t asad, Int_t chan)
{
    LKPad *pad = nullptr;
    auto padID = FindPadID(cobo,aget,asad,chan);
    if (padID>=0)
        pad = LKPadPlane::GetPad(padID);
    return pad;
}

Double_t LTPadPlane::PadDisplacement() const
{
    return sqrt(fXSpacing*fXSpacing + fRSpacing*fRSpacing);
}

bool LTPadPlane::IsInBoundary(Double_t i, Double_t j)
{
    Double_t r = TMath::Sqrt(i*i+j*j);
    if (r < fRMin || r > fRMax)
        return false;

    return true;
}

void LTPadPlane::GetSectionParameters(Int_t section, Double_t &t1, Double_t &t2, Double_t &xmin, Double_t &xmax, Double_t &ymin, Double_t &ymax)
{
    xmin=DBL_MAX, xmax=-DBL_MAX, ymin=DBL_MAX, ymax=-DBL_MAX;
         if (section==0) { xmin =-250; xmax = 250; ymin =    0; ymax = 550; t1 = +fTanPi3o8; t2 = +fTanPi5o8; }
    else if (section==1) { xmin =-550; xmax =   0; ymin =    0; ymax = 550; t1 = +fTanPi5o8; t2 = +fTanPi7o8; }
    else if (section==2) { xmin =-550; xmax =   0; ymin = -250; ymax = 250; t1 = +fTanPi7o8; t2 = +fTanPi1o8; }
    else if (section==3) { xmin =-550; xmax =   0; ymin = -550; ymax =   0; t1 = +fTanPi1o8; t2 = -fTanPi3o8; }
    else if (section==4) { xmin =-250; xmax = 250; ymin = -550; ymax =   0; t1 = -fTanPi3o8; t2 = -fTanPi5o8; }
    else if (section==5) { xmin =   0; xmax = 550; ymin = -550; ymax =   0; t1 = -fTanPi5o8; t2 = -fTanPi7o8; }
    else if (section==6) { xmin =   0; xmax = 550; ymin = -250; ymax = 250; t1 = -fTanPi7o8; t2 = +fTanPi1o8; }
    else                 { xmin =   0; xmax = 550; ymin =    0; ymax = 550; t1 = +fTanPi1o8; t2 = +fTanPi3o8; }
    //double r1 = 550;
    //double r2 = 550;
    //auto x1 = r1/sqrt(1+t1*t1); auto y1 = x1*t1;
    //auto x2 = r1/sqrt(1+t2*t2); auto y2 = x2*t2;
    //auto x3 = r2/sqrt(1+t1*t1); auto y3 = x3*t1;
    //auto x4 = r2/sqrt(1+t2*t2); auto y4 = x4*t2;
    //for (auto x : {x1,x2,x3,x4}) {
    //    if (xmin>x) xmin = x;
    //    if (xmax<x) xmax = x;
    //}
    //for (auto y : {y1,y2,y3,y4}) {
    //    if (ymin>y) ymin = y;
    //    if (ymax<y) ymax = y;
    //}
}

void LTPadPlane::CreateHistograms()
{
    if (fFramePadPlane!=nullptr)
        return;

    fFramePadPlane = new TH2Poly("frameLTPP","LAMPS TPC Pad Plane;x (mm);y (mm)",fXMin,fXMax,fYMin,fYMax);
    fHistPadPlane  = new TH2Poly("histLTPP", "LAMPS TPC Pad Plane;x (mm);y (mm)",fXMin,fXMax,fYMin,fYMax);
    fHistSideView  = new TH2D("histLTSV", "LAMPS TPC Side View;z (mm);y (mm)",fZBins,fZMin,fZMax,3*fNumLayers,fYMin,fYMax);
    fFramePadPlane -> SetStats(0);
    fHistPadPlane  -> SetStats(0);
    fHistSideView  -> SetStats(0);
    fFramePadPlane -> GetXaxis() -> SetTickLength(0.01);
    fFramePadPlane -> GetYaxis() -> SetTickLength(0.01);
    fHistPadPlane  -> GetXaxis() -> SetTickLength(0.01);
    fHistPadPlane  -> GetYaxis() -> SetTickLength(0.01);
    fHistSideView  -> GetXaxis() -> SetTickLength(0.01);
    fHistSideView  -> GetYaxis() -> SetTickLength(0.01);
    for (auto section=0; section<8; ++section) {
        Double_t t1, t2, xmin, xmax, ymin, ymax;
        GetSectionParameters(section,t1,t2,xmin,xmax,ymin,ymax);
        fHistPadPlaneSection[section] = new TH2Poly(Form("frameLTPP%d",section),Form("LAMPS TPC Pad Plane section-%d;x (mm);y (mm)",section),xmin,xmax,ymin,ymax);
        fHistPadPlaneSection[section]  -> SetStats(0);
        fHistPadPlaneSection[section]  -> SetTitle(";x (mm); y (mm)");
        //fHistPadPlaneSection[section]  -> GetXaxis() -> SetRangeUser(xmin,xmax);
        //fHistPadPlaneSection[section]  -> GetYaxis() -> SetRangeUser(xmin,xmax);
        fHistPadPlaneSection[section]  -> GetXaxis() -> SetTickLength(0.01);
        fHistPadPlaneSection[section]  -> GetYaxis() -> SetTickLength(0.01);
    }

    LKPad *pad;
    Double_t xPoints[6] = {0};
    Double_t yPoints[6] = {0};
    TIter iterPads(fChannelArray);
    while ((pad = (LKPad *) iterPads.Next())) 
    {
        auto corners = pad -> GetPadCorners();
        Int_t numCorners = corners->size();
        for (auto iCorner = 0; iCorner < numCorners; ++iCorner)
        {
            TVector2 corner = corners->at(iCorner);
            xPoints[iCorner] = corner.X();
            yPoints[iCorner] = corner.Y();
        }
        TVector2 corner = corners->at(0);
        xPoints[numCorners] = corner.X();
        yPoints[numCorners] = corner.Y();
        fHistPadPlane  -> AddBin(numCorners+1, xPoints, yPoints);
        auto section = pad -> GetSection();
        fHistPadPlaneSection[section] -> AddBin(numCorners+1, xPoints, yPoints);
    }

    for (auto section=0; section<8; ++section) {
        xPoints[0] = fPosSectionCorner[section][0].X(); yPoints[0] = fPosSectionCorner[section][0].Y();
        xPoints[1] = fPosSectionCorner[section][1].X(); yPoints[1] = fPosSectionCorner[section][1].Y();
        xPoints[2] = fPosSectionCorner[section][2].X(); yPoints[2] = fPosSectionCorner[section][2].Y();
        xPoints[3] = fPosSectionCorner[section][3].X(); yPoints[3] = fPosSectionCorner[section][3].Y();
        xPoints[4] = fPosSectionCorner[section][0].X(); yPoints[4] = fPosSectionCorner[section][0].Y();
        fFramePadPlane -> AddBin(5, xPoints, yPoints);
    }

    xPoints[0] = fXMin+50;  yPoints[0] = fYMax-50;
    xPoints[1] = fXMin+50;  yPoints[1] = fYMax-250;
    xPoints[2] = fXMin+250; yPoints[2] = fYMax-50;
    xPoints[3] = fXMin+50;  yPoints[3] = fYMax-50;
    fBinNumberSideView = fFramePadPlane -> AddBin(4, xPoints, yPoints); 

    TIter nextBin(fFramePadPlane -> GetBins());
    while (auto bin = (TH2PolyBin*) nextBin()) {
        auto graphBin = (TGraph *) bin -> GetPolygon();
        graphBin -> SetLineColor(kGray+1);
        //if (removeBinBoundary) graphBin -> SetLineWidth(0);
    }

    fHistChannel1 = new TH1D("hist_channel_1","channel buffer;time-bucket;charge",512,0,512);
    fHistChannel2 = new TH1D("hist_channel_2","channel buffer;time-bucket;charge",512,0,512);
    for (auto histChannel : {fHistChannel1,fHistChannel2}) {
        histChannel -> SetStats(0);
        histChannel -> GetXaxis() -> SetLabelSize(0.065);
        histChannel -> GetYaxis() -> SetLabelSize(0.065);
        histChannel -> GetXaxis() -> SetTitleSize(0.065);
        histChannel -> GetYaxis() -> SetTitleSize(0.065);
        histChannel -> GetXaxis() -> SetTitleOffset(1.20);
        //histChannel -> GetYaxis() -> SetTitleOffset(0.68);
        histChannel -> GetYaxis() -> SetTitleOffset(0.80);
    }
}

TH2* LTPadPlane::GetHist(Option_t *option)
{
    CreateHistograms();

    TString optionString(option);
    Int_t selectSection = -1;
    if (optionString.IsNull())
        selectSection = optionString.Atoi();
    return GetHist(selectSection);
}

TH2* LTPadPlane::GetHist(Int_t selectSection)
{
    CreateHistograms();

    if (selectSection>=0)
        return fHistPadPlaneSection[selectSection];
    return fHistPadPlane;
}

void LTPadPlane::DrawFrame(Option_t *)
{
}

TCanvas *LTPadPlane::GetCanvas(Option_t *)
{
    if (fCanvas==nullptr) {
        fCanvas = LKWindowManager::GetWindowManager() -> CanvasResize("LAMPS_TPC_PadPlane",1100,800,0.9);
        auto pad1 = new TPad("pad1","",0,230./700,0.5,1);
        pad1 -> SetMargin(0.12,0.15,0.1,0.1);
        pad1 -> SetNumber(1);
        pad1 -> AddExec("ex", "LTPadPlane::MouseClickEvent1()");
        pad1 -> Draw();
        auto pad2 = new TPad("pad1","",0.5,230./700,1,1);
        pad2 -> SetMargin(0.12,0.15,0.1,0.1);
        pad2 -> SetNumber(2);
        pad2 -> AddExec("ex", "LTPadPlane::MouseClickEvent2()");
        pad2 -> Draw();
        auto pad3 = new TPad("pad2","",0,0,0.5,230./700);
        pad3 -> SetMargin(0.12,0.05,0.20,0.12);
        pad3 -> SetNumber(3);
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

LKPad *LTPadPlane::NewPad(Int_t s, Int_t r, Int_t l)
{
    auto pad = new LKPad();
    pad -> SetSectionRowLayer(s, r, l);
    auto padID = fChannelArray -> GetEntriesFast();
    pad -> SetPadID(padID);
    fChannelArray -> Add(pad);
    return pad;
}

void LTPadPlane::SetNeighborPads(LKPad *pad0, LKPad *pad1)
{
    pad0 -> AddNeighborPad(pad1);
    pad1 -> AddNeighborPad(pad0);
}

Int_t LTPadPlane::FindSection(Double_t i, Double_t j)
{
    if (j > fTanPi3o8*i) {
        if (j > fTanPi1o8*i) {
            if (j > fTanPi7o8*i) {
                if (j > fTanPi5o8*i) {
                    return 0;
                } else return 1;
            } else return 2;
        } else return 3;
    }
    else
    {
        if (j < fTanPi1o8*i) {
            if (j < fTanPi7o8*i) {
                if (j < fTanPi5o8*i) {
                    return 4;
                } else return 5;
            } else return 6;
        } else return 7;
    }
}

bool LTPadPlane::SetDataFromBranch()
{
    if (fRun==nullptr)
        return false;
    fBufferArray = fRun -> GetBranchA("RawData");
    fHitArray = fRun -> GetBranchA("Hit");
    return true;
}

void LTPadPlane::FillDataToHist()
{
    lk_info << endl;
    if (fPar -> CheckPar("eve/planeFillType")) {
        fFillType = fPar -> GetParString("eve/planeFillType");
        if (fFillType!="Hit" || fFillType!="Buffer") {
            lk_warning << "eve/planeFillType must be Buffer or Hit" << endl;
            if (fBufferArray!=nullptr) fFillType = "Buffer";
            else if (fHitArray!=nullptr) fFillType = "Hit";
        }
    }

    if (fFillType=="Buffer")
    {
        if (fBufferArray==nullptr) {
            lk_error << "Fill type is " << fFillType << " but buffer array is null!" << endl;
            return;
        }

        fHistPadPlane -> Reset("ICES");
        fHistSideView -> Reset("ICES");
        for (auto section=0; section<8; ++section)
            fHistPadPlaneSection[section] -> Reset("ICES");

        auto numChannels = fBufferArray -> GetEntries();
        lk_info << "# of channels: " << numChannels << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (GETChannel*) fBufferArray -> At(iChannel);
            auto padID = channel -> GetChan2();
            if (padID<0) {
                auto cobo = channel -> GetCobo();
                auto asad = channel -> GetAsad();
                auto aget = channel -> GetAget();
                auto chan = channel -> GetChan();
                padID = FindPadID(cobo,aget,asad,chan);
            }
            if (padID<0) {
                lk_debug << padID << endl;
                continue;
            }
            auto charge = channel -> GetEnergy();
            int max[3] = {0};
            int bin[3] = {0};
            if (charge<=0) {
                charge = 0;
                auto buffer = channel -> GetWaveformY();
                for (auto i=0; i<512; ++i) {
                    auto value = buffer[i];
                    charge += value;
                         if (max[0]<value) { max[0] = value; bin[0] = i; }
                    else if (max[1]<value) { max[1] = value; bin[1] = i; }
                    else if (max[2]<value) { max[2] = value; bin[2] = i; }
                }
            }
            lk_debug << padID << " " << charge << " " << max[0] << " " << max[1] << " " << max[2] << endl;
            if (charge<0)
                continue;
            auto pad = LKPadPlane::GetPad(padID);
            auto x = pad -> GetX();
            auto y = pad -> GetY();
            auto section = pad -> GetSection();
            auto hbin = padID+5;
            //fHistPadPlane -> //SetBinContent(hbin,charge);
            fHistPadPlane -> Fill(x,y,charge);
            fHistPadPlaneSection[section] -> Fill(x,y,charge);
            //fHistSideView -> Fill(z,y,charge);
        }
    }
    else if (fFillType=="Hit")
    {
        if (fHitArray==nullptr) {
            lk_error << "Fill type is " << fFillType << " but hit array is null!" << endl;
            return;
        }
    }
}

void LTPadPlane::Draw(Option_t *option)
{
    CreateHistograms();
    FillDataToHist();

    auto cvs = GetCanvas();

    cvs -> cd(1);
    if (fHistPadPlane->GetEntries()==0)
        fFramePadPlane -> Draw();
    else {
        fFramePadPlane -> Draw();
        fHistPadPlane -> SetMinimum(0.1);
        fHistPadPlane -> Draw("same colz");
    }

    cvs -> cd(2);
    fHistPadPlaneSection[0] -> SetMinimum(0.1);
    fHistPadPlaneSection[0] -> Draw("same colz");

    cvs -> cd(3);
    fHistChannel1 -> Draw();
    cvs -> cd(3) -> Modified();
    cvs -> cd(3) -> Update();
    auto ttt1 = (TPaveText*) (cvs->cd(3)->GetListOfPrimitives()) -> FindObject("title");
    ttt1 -> SetTextSize(0.065);
    ttt1 -> SetTextAlign(12);
    cvs -> cd(3) -> Modified();
    cvs -> cd(3) -> Update();

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

    if (binCurr<=0)
        return;

    gPad -> SetUniqueID(binCurr);
    gPad -> GetCanvas() -> SetClickSelected(nullptr);
    LTPadPlane::GetPlane() -> ZoomInWindow(binCurr, xOnClick, yOnClick);
}

void LTPadPlane::MouseClickEvent2()
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

    LTPadPlane::GetPlane() -> SelectAndDrawChannel(binCurr);
}

void LTPadPlane::ZoomInWindow(Int_t bin, Double_t x, Double_t y)
{
    GetCanvas() -> cd(2);
    if (bin==fBinNumberSideView) {
        fHistSideView -> Draw("colz");

        fGraphSectionBoundary -> Set(0);
        fGraphSectionBoundary -> SetPoint(0, fXMin+50,  fYMax-50);
        fGraphSectionBoundary -> SetPoint(1, fXMin+50,  fYMax-250);
        fGraphSectionBoundary -> SetPoint(2, fXMin+250, fYMax-50);
        fGraphSectionBoundary -> SetPoint(3, fXMin+50,  fYMax-50);
        fGraphSectionBoundary -> SetPoint(4, fXMin+50,  fYMax-50);
    }
    else {
        auto section = FindSection(x,y);
        if (fHistPadPlaneSection[section]->GetEntries()==0)
            fHistPadPlaneSection[section] -> Draw();
        else
            fHistPadPlaneSection[section] -> Draw("colz");
        Double_t t1, t2, xmin, xmax, ymin, ymax;
        GetSectionParameters(section,t1,t2,xmin,xmax,ymin,ymax);
        if (fGraphSectionBoundary==nullptr) {
            fGraphSectionBoundary = new TGraph();
            fGraphSectionBoundary -> SetLineColor(kBlue);
            fGraphSectionBoundary -> SetLineWidth(3);
        }
        fGraphSectionBoundary -> Set(0);
        fGraphSectionBoundary -> SetPoint(0,fPosSectionCorner[section][0].X(), fPosSectionCorner[section][0].Y());
        fGraphSectionBoundary -> SetPoint(1,fPosSectionCorner[section][1].X(), fPosSectionCorner[section][1].Y());
        fGraphSectionBoundary -> SetPoint(2,fPosSectionCorner[section][2].X(), fPosSectionCorner[section][2].Y());
        fGraphSectionBoundary -> SetPoint(3,fPosSectionCorner[section][3].X(), fPosSectionCorner[section][3].Y());
        fGraphSectionBoundary -> SetPoint(4,fPosSectionCorner[section][0].X(), fPosSectionCorner[section][0].Y());
    }
    GetCanvas() -> cd(1);
    fGraphSectionBoundary -> Draw("samel");
}

void LTPadPlane::SelectAndDrawChannel(int bin)
{
}

/*
{
    //if (isChain) lk_info << "SelectAndDrawChannel (Chain) bin = " << bin << endl;
    //else         lk_info << "SelectAndDrawChannel (Strip) bin = " << bin << endl;

    auto cvsPlane = fCanvas -> cd(3);
    auto cvsChannel = fCanvas -> cd(4);
    if (isChain) {
        cvsPlane = fCanvas -> cd(1);
        cvsChannel = fCanvas -> cd(2);
    }

    bool existHitArray = false;
    bool existBufferArray = (fBufferArray!=nullptr);
    if (isChain) existHitArray = (fHitCenterArray!=nullptr&&fHitLChainArray!=nullptr&&fHitRChainArray!=nullptr);
    else         existHitArray = (fHitCenterArray!=nullptr&&fHitLStripArray!=nullptr&&fHitRStripArray!=nullptr);

    fHitArrayList.clear();
    fHitArrayList.push_back(fHitCenterArray);
    if (isChain) {
        fHitArrayList.push_back(fHitLChainArray);
        fHitArrayList.push_back(fHitRChainArray);
    }
    else {
        fHitArrayList.push_back(fHitLStripArray);
        fHitArrayList.push_back(fHitRStripArray);
    }
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
            if (isChain) bin = fMapCAACToBinChain[caac];
            else         bin = fMapCAACToBinStrip[caac];
            if (isChain) lk_info << "SelectAndDrawChannel (Chain) CAAC=" << Form("%05d",bin) << endl;
            else         lk_info << "SelectAndDrawChannel (Strip) CAAC=" << Form("%05d",bin) << endl;
        }
    }

    if (bin<0) {
        lk_warning << "bin<0" << bin << endl;
        return;
    }

    double x0,x1,x2,z0,z1,z2;
    TGraph* graphBoundary;
    if (isChain) {
        x1 = fMapBinToX1Chain[bin];
        x2 = fMapBinToX2Chain[bin];
        z1 = fMapBinToZ1Chain[bin];
        z2 = fMapBinToZ2Chain[bin];
        graphBoundary = fGraphChannelBoundaryChain;
    }
    else {
        x1 = fMapBinToX1Strip[bin];
        x2 = fMapBinToX2Strip[bin];
        z1 = fMapBinToZ1Strip[bin];
        z2 = fMapBinToZ2Strip[bin];
        graphBoundary = fGraphChannelBoundaryStrip;
    }
    x0 = (x1 + x2)/2.;
    z0 = (z1 + z2)/2.;
    graphBoundary -> Set(0);
    graphBoundary -> SetPoint(0,x1,z1);
    graphBoundary -> SetPoint(1,x2,z1);
    graphBoundary -> SetPoint(2,x2,z2);
    graphBoundary -> SetPoint(3,x1,z2);
    graphBoundary -> SetPoint(4,x1,z1);
    graphBoundary -> SetLineColor(kBlue);
    cvsPlane -> cd();
    graphBoundary -> Draw("samel");

    fCAAC = -1;
    if (isChain) fCAAC = fMapBinToCAACChain[bin];
    else         fCAAC = fMapBinToCAACStrip[bin];

    auto texat = (TexAT2*) fDetector;
    
    TH1D* histChannel = nullptr;
    if (isChain) histChannel = fHistChannel1;
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
                    //cvsPlane -> cd();
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
*/
