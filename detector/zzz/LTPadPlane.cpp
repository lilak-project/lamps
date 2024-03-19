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
#include "LKTracklet.h"

#include "LKWindowManager.h"
#include "GETChannel.h"
#include "LKHit.h"
#include "TStyle.h"

#include <iostream>
using namespace std;

ClassImp(LTPadPlane)

LTPadPlane* LTPadPlane::fInstance = nullptr;
LTPadPlane* LTPadPlane::GetPlane() { return fInstance; }

LTPadPlane::LTPadPlane()
    :LKDetectorPlane("LTPadPlane", "pad plane with rectangular pads following the circle line for LAMPS-TPC")
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

    fMapSLRPToPadID = new int***[8];
    for (int section=0; section<8; ++section) {
        fMapSLRPToPadID[section] = new int**[42];
        for (int layer=0; layer<42; ++layer) {
            fMapSLRPToPadID[section][layer] = new int*[56];
            for (int row=0; row<56; ++row) {
                fMapSLRPToPadID[section][layer][row] = new int[2];
                for (int pm=0; pm<2; ++pm) {
                    fMapSLRPToPadID[section][layer][row][pm] = -1;
                }
            }
        }
    }

    for (auto section=0; section<8; ++section)
    {
        Double_t r1 = fRMin-10;
        Double_t r2 = 560;
        Double_t x1 = -r2/sqrt(1+fTanPi3o8*fTanPi3o8);
        Double_t x2 = -r1/sqrt(1+fTanPi3o8*fTanPi3o8);
        Double_t x3 = -x2;
        Double_t x4 = -x1;
        Double_t y1 = -fTanPi3o8*x1;
        Double_t y2 = -fTanPi3o8*x2;
        Double_t y3 = +fTanPi3o8*x3;
        Double_t y4 = +fTanPi3o8*x4;
        TVector2 sectionCorner1(x1, y1);
        TVector2 sectionCorner2(x2, y2);
        TVector2 sectionCorner3(x3, y3);
        TVector2 sectionCorner4(x4, y4);

        Double_t phiSection = section * TMath::Pi()/4.;
        fPosSectionCorner[section][0] = sectionCorner1.Rotate(phiSection);
        fPosSectionCorner[section][1] = sectionCorner2.Rotate(phiSection);
        fPosSectionCorner[section][2] = sectionCorner3.Rotate(phiSection);
        fPosSectionCorner[section][3] = sectionCorner4.Rotate(phiSection);
    }

    Double_t xCorner[5] = {0};
    Double_t yCorner[5] = {0};

    for (auto layer=0; layer<fNumLayers; ++layer)
        fStartRowFrom[layer] = 1;
    fStartRowFrom[41] = 55;
    fStartRowFrom[40] = 46;
    fStartRowFrom[39] = 35;
    fStartRowFrom[38] = 17;

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
            if (layer==41&&iRow<fStartRowFrom[41]) continueRow = true;
            if (layer==40&&iRow<fStartRowFrom[40]) continueRow = true;
            if (layer==39&&iRow<fStartRowFrom[39]) continueRow = true;
            if (layer==38&&iRow<fStartRowFrom[38]) continueRow = true;
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

            if (row>=fNumHalfRowsInLayerInput[layer])
            //if (yPad < GetPadCenterYBoundAtX(xPad))
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

            if (fDoCutSideBoundary && (yCorner[3] < GetPadCutBoundaryYAtX(xCorner[3])))
            {
                sideBoundaryCutWasMade = true;
                if (xCorner[0] > GetPadCutBoundaryXAtR(radius2)) {
                    xCorner[0] = GetPadCutBoundaryXAtR(radius2);
                    yCorner[0] = GetPadCutBoundaryYAtX(xCorner[0]);
                    yCorner[2] = GetPadCutBoundaryYAtX(xCorner[2]);
                    numCorners = 3;
                }
                else if (yCorner[2] < GetPadCutBoundaryYAtX(xCorner[2])) {
                    yCorner[2] = GetPadCutBoundaryYAtX(xCorner[2]);
                    yCorner[3] = GetPadCutBoundaryYAtX(xCorner[3]);
                }
                else if (yCorner[3] < GetPadCutBoundaryYAtX(xCorner[3])) {
                    xCorner[3] = GetPadCutBoundaryXAtR(radius1);
                    yCorner[3] = GetPadCutBoundaryYAtX(xCorner[3]);
                    if (layer>37&&row>fNumHalfRowsInLayerInput[layer]-3) {
                        //+170  /home/ejungwoo/lilak/lamps_2023/detector/LTPadPlane.cpp # 38 39 40
                    }
                    if (row==fNumHalfRowsInLayerInput[layer]) {
                        xCorner[0] = GetPadCutBoundaryXAtR(radius2);
                        yCorner[0] = GetPadCutBoundaryYAtX(xCorner[0]);
                    }
                    else {
                        xCorner[4] = xCorner[0];
                        yCorner[4] = GetPadCutBoundaryYAtX(xCorner[4]);
                        numCorners = 5;
                    }
                }
            }
            if (fDoCutTopBoundary && yCorner[0] > fRTopCut)
            {
                yCorner[0] = fRTopCut;
                yCorner[1] = fRTopCut;
            }

            if (layer<fNumLayers-1 && fStartRowFrom[layer+1]>iRow)
            {
                yCorner[0] = fRTopCut;
                yCorner[1] = fRTopCut;
            }

            ///////////////////////////////////////////////////////////////////

            for (auto section=0; section<8; ++section)
            {
                Double_t phiSection = section * TMath::Pi()/4.;

                TVector2 point;
                for (int iRL : {0,1})
                {
                    int signRL = (iRL==0?1:-1);
                    int pm = (iRL==0?1:0);
                    auto pad = NewPad(section, layer, signRL*row);
                    auto padID = pad -> GetPadID();
                    fMapSLRPToPadID[section][layer][row][pm] = padID;
                    //auto slr = iRL + 10*row + 1000*layer + 100000*section;
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
                auto padL = LKDetectorPlane::GetPad(section,layer,-row);
                auto padR = LKDetectorPlane::GetPad(section,layer,+row);

                // pad below
                if (layer>0)
                {
                    auto rowNb = row + fNumSkippedHalfRows[layer] - fNumSkippedHalfRows[layer-1];

                    auto rowL1 = -rowNb - 1;
                    if (abs(rowL1)<=fNumHalfRowsInLayer[layer-1]) {
                        auto padBelowL1 = LKDetectorPlane::GetPad(section,layer-1,rowL1);
                        if (padBelowL1!=nullptr) SetNeighborPads(padL,padBelowL1);
                    }

                    auto rowL2 = -rowNb + 1; if (rowL2==0) rowL2 += 1;
                    auto padBelowL2 = LKDetectorPlane::GetPad(section,layer-1,rowL2);
                    if (padBelowL2!=nullptr) SetNeighborPads(padL,padBelowL2);

                    auto rowR1 = rowNb - 1; if (rowR1==0) rowR1 -= 1;
                    auto padBelowR1 = LKDetectorPlane::GetPad(section,layer-1,rowR1);
                    if (padBelowR1!=nullptr) SetNeighborPads(padR,padBelowR1);

                    auto rowR2 = rowNb + 1;
                    if (abs(rowR2)<=fNumHalfRowsInLayer[layer-1]) {
                        auto padBelowR2 = LKDetectorPlane::GetPad(section,layer-1,rowR2);
                        if (padBelowR2!=nullptr) SetNeighborPads(padR,padBelowR2);
                    }

                    if (rowNb<=0||rowNb>fNumHalfRowsInLayer[layer-1])
                        continue;

                    auto padBelowL = LKDetectorPlane::GetPad(section,layer-1,-rowNb);
                    SetNeighborPads(padL,padBelowL);

                    auto padBelowR = LKDetectorPlane::GetPad(section,layer-1,+rowNb);
                    SetNeighborPads(padR,padBelowR);

                }
            }
        }
    }

    /*
    auto numPads = GetNumPads();
    for (auto iPad=0; iPad<numPads; ++iPad) {
        auto pad = GetPad(iPad);
        auto nbPadArray = pad -> GetNeighborPadArray();
        auto numNbPads = nbPadArray -> size();
        for (auto iPad=0; iPad<numNbPads; ++iPad)
        {
            auto pad1 = nbPadArray -> at(iPad);
            for (auto jPad=0; jPad<numNbPads; ++jPad)
            {
                if (iPad==jPad) continue;
                auto pad2 = nbPadArray -> at(jPad);
            }
        }
    }
    */

    fMapCAACToPadID = new int***[22];
    for (int cobo=0; cobo<22; ++cobo) {
        fMapCAACToPadID[cobo] = new int**[4];
        for (int asad=0; asad<4; ++asad) {
            fMapCAACToPadID[cobo][asad] = new int*[4];
            for (int aget=0; aget<4; ++aget) {
                fMapCAACToPadID[cobo][asad][aget] = new int[68];
                for (int chan=0; chan<68; ++chan) {
                    fMapCAACToPadID[cobo][asad][aget][chan] = -1;
                }
            }
        }
    }

    if (fPar -> CheckPar("LTPadPlane/position_map"))
    {
        double x;
        double y;
        int cobo; // 0-21
        int asad; // 0-3
        int aget; // 0-3
        int channelID; // 0-67
        int section; // 1-8
        int padIDLocal; // 0-2697
        int padID; // global padID create from this class LTPadPlane. This is the order of pad that has been created and stored separately.
        int layer; // 0-41
        int asad2; // 1-11
        auto mapFileName = fPar -> GetParString("LTPadPlane/position_map");
        lk_info << "pad electronics mapping: " << mapFileName << endl;
        ifstream fileMap(mapFileName);
        //for (auto i=0; i<21584; ++i) fileMap >> cobo >> asad >> aget >> channelID >> x >> y >> section >> padIDLocal >> layer >> asad2 >> padID;
        while (fileMap >> cobo >> asad >> aget >> channelID >> x >> y >> section >> padIDLocal >> layer >> asad2 >> padID)
        {
            if (padID>=0) {
                auto pad = GetPad(padID);
                pad -> SetCoboID(cobo);
                pad -> SetAsadID(asad);
                pad -> SetAgetID(aget);
                pad -> SetChannelID(channelID);
                padID = FindPadID(x,y); // XXX
                fMapCAACToPadID[cobo][asad][aget][channelID] = padID;
                //e_cout << "input " << cobo << " " << aget << " " << asad << " " << channelID << " " << fMapCAACToPadID[cobo][aget][asad][channelID]  << endl;
            }
        }
    }

    //for (int cobo=0; cobo<22; ++cobo)
    //    for (int asad=0; asad<4; ++asad)
    //        for (int aget=0; aget<4; ++aget)
    //            for (int chan=0; chan<68; ++chan)
    //                e_cout << cobo << " " << asad << " " << aget << " " << chan << " " << fMapCAACToPadID[cobo][asad][aget][chan] << endl;

    if (fRun!=nullptr)
    {
        fRawDataArray = fRun -> GetBranchA("RawData");
        fHitArray = fRun -> GetBranchA("Hit");
        fTrackArray = fRun -> GetBranchA("Track");
    }

    fAxis1 = LKVector3::kX;
    fAxis2 = LKVector3::kY;
    fAxis3 = LKVector3::kZ;
    fAxisDrift = LKVector3::kZ;
    fTbToLength = 2.34375;
    fPosition = 900;

    fChannelAnalyzer = GetChannelAnalyzer();

    return true;
}

Int_t LTPadPlane::FindPadID(Int_t section, Int_t layer, Int_t row)
{
    int idLayer = section*fNumPadsDownToLayer[layer];
    if (layer<fNumLayers) idLayer = idLayer + (8-section)*fNumPadsDownToLayer[layer+1];
    int idRow = ((row>0) ? (fNumHalfRowsInLayer[layer] - row) : (fNumHalfRowsInLayer[layer] - row - 1));
    int padID = idLayer + idRow;

    //int padID2 = -1;
    //int pm = int (row>0?1:0);
    //row = abs(row);
    //if (section<0||section>=8) padID2 = -1;
    //if (layer<0||layer>=fNumLayers) padID2 = -2;
    //if (row==0||abs(row)>fNumHalfRowsInLayer[layer]) padID2 = -3;
    //padID2 = fMapSLRPToPadID[section][layer][row][pm];
    //if (padID2!=padID)
    //    lk_debug << section << " " << layer << " " << row << " " << pm << " " << padID2 << " " << padID << endl;

    return padID;
}

Int_t LTPadPlane::FindPadID(Int_t cobo, Int_t asad, Int_t aget, Int_t chan)
{
    auto padID = fMapCAACToPadID[cobo][asad][aget][chan];
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
        return -2;

    Int_t pm = 1;
    if (xRotatedToSec0 < 0) {
        xRotatedToSec0 = -xRotatedToSec0;
        pm = -1;
    }

    if (xRotatedToSec0 < .5*fPadGap)
        return -3;

    Double_t xFromRow0LeftEdge = xRotatedToSec0 - .5*fPadGap;
    Int_t row = (Int_t) (xFromRow0LeftEdge / fXSpacing) + 1;
    if (xFromRow0LeftEdge - (row-1)*fXSpacing > fPadWidth)
        return -4;

    row = row - fNumSkippedHalfRows[layer];

    if (row > fNumHalfRowsInLayer[layer])
        return -5;

    return FindPadID(section,layer,pm*row);
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
    auto x1 = fPosSectionCorner[section][0].X();
    auto x2 = fPosSectionCorner[section][1].X();
    auto x3 = fPosSectionCorner[section][2].X();
    auto x4 = fPosSectionCorner[section][3].X();
    auto y1 = fPosSectionCorner[section][0].Y();
    auto y2 = fPosSectionCorner[section][1].Y();
    auto y3 = fPosSectionCorner[section][2].Y();
    auto y4 = fPosSectionCorner[section][3].Y();
    xmin=DBL_MAX, xmax=-DBL_MAX, ymin=DBL_MAX, ymax=-DBL_MAX;
    for (auto x : {x1,x2,x3,x4}) {
        if (xmin>x) xmin = x;
        if (xmax<x) xmax = x;
    }
    for (auto y : {y1,y2,y3,y4}) {
        if (ymin>y) ymin = y;
        if (ymax<y) ymax = y;
    }
}

void LTPadPlane::CreateHistograms()
{
    if (fFramePadPlane!=nullptr)
        return;

    fPar -> UpdatePar(fHistZMin,"LTPadPlane/histZMin");

    fFramePadPlane = new TH2Poly("frameLTPP","LAMPS TPC Pad Plane;x (mm);y (mm)",fXMin,fXMax,fYMin,fYMax);
    fFramePadPlane -> SetStats(0);
    fFramePadPlane -> GetXaxis() -> SetTickLength(0.01);
    fFramePadPlane -> GetYaxis() -> SetTickLength(0.01);

    fHistPadPlane = new TH2Poly("histLTPP", "LAMPS TPC Pad Plane;x (mm);y (mm)",fXMin,fXMax,fYMin,fYMax);
    fHistPadPlane -> SetStats(0);
    fHistPadPlane -> GetXaxis() -> SetTickLength(0.01);
    fHistPadPlane -> GetYaxis() -> SetTickLength(0.01);
    fHistPadPlane -> SetMinimum(fHistZMin);

    fHistTopView = new TH2D("histLTTV", "LAMPS TPC Top View;z (mm);x (mm)",fZBins,fZMin,fZMax,fNumLayers,fXMin,fXMax);
    fHistTopView -> SetStats(0);
    fHistTopView -> GetXaxis() -> SetTickLength(0.01);
    fHistTopView -> GetYaxis() -> SetTickLength(0.01);
    fHistTopView -> SetMinimum(fHistZMin);

    fHistSideView = new TH2D("histLTSV", "LAMPS TPC Side View;z (mm);y (mm)",fZBins,fZMin,fZMax,fNumLayers,fYMin,fYMax);
    fHistSideView -> SetStats(0);
    fHistSideView -> GetXaxis() -> SetTickLength(0.01);
    fHistSideView -> GetYaxis() -> SetTickLength(0.01);
    fHistSideView -> SetMinimum(fHistZMin);

    fHistWaveform = new TH2D("histWaveform", "All channels;tb;charge", fNumTbs, 0, fNumTbs, fMaxWaveformY, 0, fMaxWaveformY);
    fHistWaveform -> SetStats(0);
    fHistWaveform -> GetXaxis() -> SetTickLength(0.01);
    fHistWaveform -> GetYaxis() -> SetTickLength(0.01);

    for (auto section=0; section<8; ++section) {
        Double_t t1, t2, xmin, xmax, ymin, ymax;
        GetSectionParameters(section,t1,t2,xmin,xmax,ymin,ymax);
        fHistPadPlaneSection[section] = new TH2Poly(Form("histLTPP_S%d",section),Form("S%d;x (mm);y (mm)",section),xmin,xmax,ymin,ymax);
        fHistPadPlaneSection[section] -> SetStats(0);
        fHistPadPlaneSection[section] -> GetXaxis() -> SetTickLength(0.01);
        fHistPadPlaneSection[section] -> GetYaxis() -> SetTickLength(0.01);
        fHistPadPlaneSection[section] -> SetMinimum(fHistZMin);

        fFramePadPlaneSection[section] = new TH2Poly(Form("frameLTPP_S%d",section),Form("S%d;x (mm);y (mm)",section),xmin,xmax,ymin,ymax);
        fFramePadPlaneSection[section] -> SetStats(0);
        fFramePadPlaneSection[section] -> GetXaxis() -> SetTickLength(0.01);
        fFramePadPlaneSection[section] -> GetYaxis() -> SetTickLength(0.01);
        fFramePadPlaneSection[section] -> SetMinimum(fHistZMin);
    }

    fGraphSectionBoundary1 = new TGraph();
    fGraphSectionBoundary1 -> SetLineColor(kBlue);
    fGraphSectionBoundary1 -> SetLineWidth(3);

    fGraphSectionBoundary2 = new TGraph();
    fGraphSectionBoundary2 -> SetLineColor(kRed);
    fGraphSectionBoundary2 -> SetLineWidth(3);

    fGraphPadBoundary = new TGraph();
    fGraphPadBoundary -> SetLineColor(kRed);
    fGraphPadBoundary -> SetLineWidth(3);
    for (auto i=0; i<10; ++i) {
        fGraphPadBoundaryNb[i] = new TGraph();
        fGraphPadBoundaryNb[i] -> SetLineColor(kRed);
        fGraphPadBoundaryNb[i] -> SetLineWidth(2);
        fGraphPadBoundaryNb[i] -> SetLineStyle(2);
    }

    LKPad *pad;
    Double_t xPoints[10] = {0};
    Double_t yPoints[10] = {0};
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
        fFramePadPlaneSection[section] -> AddBin(numCorners+1, xPoints, yPoints);
    }

    for (auto section=0; section<8; ++section) {
        xPoints[0] = fPosSectionCorner[section][0].X(); yPoints[0] = fPosSectionCorner[section][0].Y();
        xPoints[1] = fPosSectionCorner[section][1].X(); yPoints[1] = fPosSectionCorner[section][1].Y();
        xPoints[2] = fPosSectionCorner[section][2].X(); yPoints[2] = fPosSectionCorner[section][2].Y();
        xPoints[3] = fPosSectionCorner[section][3].X(); yPoints[3] = fPosSectionCorner[section][3].Y();
        xPoints[4] = fPosSectionCorner[section][0].X(); yPoints[4] = fPosSectionCorner[section][0].Y();
        fFramePadPlane -> AddBin(5, xPoints, yPoints);
    }

    //for (auto section=0; section<8; ++section) {
    //    xPoints[section] = fPosSectionCorner[section][2].X();
    //    yPoints[section] = fPosSectionCorner[section][2].Y();
    //}
    //xPoints[8] = fPosSectionCorner[0][2].X();
    //yPoints[8] = fPosSectionCorner[0][2].Y();
    //fBinZoomZoomButton = fFramePadPlane -> AddBin(9, xPoints, yPoints);

    xPoints[0] = fXMax-50;  yPoints[0] = fYMax-50;
    xPoints[1] = fXMax-50;  yPoints[1] = fYMax-250;
    xPoints[2] = fXMax-250; yPoints[2] = fYMax-50;
    xPoints[3] = fXMax-50;  yPoints[3] = fYMax-50;
    fBinNumberTopView = fFramePadPlane -> AddBin(4, xPoints, yPoints); 

    xPoints[0] = fXMin+50;  yPoints[0] = fYMin+50;
    xPoints[1] = fXMin+50;  yPoints[1] = fYMin+250;
    xPoints[2] = fXMin+250; yPoints[2] = fYMin+50;
    xPoints[3] = fXMin+50;  yPoints[3] = fYMin+50;
    fBinNumberSideView = fFramePadPlane -> AddBin(4, xPoints, yPoints); 

    TIter nextBin(fFramePadPlane -> GetBins());
    while (auto bin = (TH2PolyBin*) nextBin()) {
        auto graphBin = (TGraph *) bin -> GetPolygon();
        graphBin -> SetLineColor(kGray+1);
    }

    for (auto section=0; section<8; ++section) {
        TIter nextBin(fFramePadPlaneSection[section] -> GetBins());
        while (auto bin = (TH2PolyBin*) nextBin()) {
            auto graphBin = (TGraph *) bin -> GetPolygon();
            graphBin -> SetLineColor(kGray+1);
        }
    }

    fHistChannel1 = new TH1D("hist_channel_1","channel buffer;time-bucket;charge",fNumTbs,0,fNumTbs);
    fHistChannel2 = new TH1D("hist_channel_2","channel buffer;time-bucket;charge",fNumTbs,0,fNumTbs);
    for (TH1* histChannel : {(TH1*)fHistChannel1,(TH1*)fHistChannel2,(TH1*)fHistWaveform,(TH1*)fHistTopView,(TH1*)fHistSideView}) {
        histChannel -> SetStats(0);
        histChannel -> GetXaxis() -> SetLabelSize(0.062);
        histChannel -> GetYaxis() -> SetLabelSize(0.062);
        histChannel -> GetXaxis() -> SetTitleSize(0.062);
        histChannel -> GetYaxis() -> SetTitleSize(0.062);
        histChannel -> GetXaxis() -> SetTitleOffset(1.20);
        histChannel -> GetYaxis() -> SetTitleOffset(1.00);
    }

    double binTextSize = 6.0;
    double ctrlLabelSize = 0.15;

    fHistCtrlEv1 = new TH2D("LTPadPlane_CtrlEv1","",7,0,7,1,0,1);
    fHistCtrlEv1 -> SetStats(0);
    fBinCtrlFrst = fHistCtrlEv1 -> GetBin(1,1);
    fBinCtrlPr50 = fHistCtrlEv1 -> GetBin(2,1);
    fBinCtrlPrev = fHistCtrlEv1 -> GetBin(3,1);
    fBinCtrlCurr = fHistCtrlEv1 -> GetBin(4,1);
    fBinCtrlNext = fHistCtrlEv1 -> GetBin(5,1);
    fBinCtrlNe50 = fHistCtrlEv1 -> GetBin(6,1);
    fBinCtrlLast = fHistCtrlEv1 -> GetBin(7,1);
    fHistCtrlEv1 -> GetXaxis() -> SetTickSize(0);
    fHistCtrlEv1 -> GetYaxis() -> SetTickSize(0);
    fHistCtrlEv1 -> GetYaxis() -> SetBinLabel(1,"");
    fHistCtrlEv1 -> GetXaxis() -> SetLabelSize(ctrlLabelSize);
    fHistCtrlEv1 -> GetXaxis() -> SetBinLabel(1,"First");
    fHistCtrlEv1 -> GetXaxis() -> SetBinLabel(2,"-50");
    fHistCtrlEv1 -> GetXaxis() -> SetBinLabel(3,"Prev.");
    fHistCtrlEv1 -> GetXaxis() -> SetBinLabel(4,"Current");
    fHistCtrlEv1 -> GetXaxis() -> SetBinLabel(5,"Next");
    fHistCtrlEv1 -> GetXaxis() -> SetBinLabel(6,"+50");
    fHistCtrlEv1 -> GetXaxis() -> SetBinLabel(7,"Last");
    fHistCtrlEv1 -> SetBinContent(fBinCtrlFrst,0);
    if (fRun!=nullptr)
        fHistCtrlEv1 -> SetBinContent(fBinCtrlLast,fRun->GetNumEvents()-1);
    else {
        fHistCtrlEv1 -> SetBinContent(fBinCtrlLast,0);
    }
    fHistCtrlEv1 -> SetMarkerSize(binTextSize);
    fHistCtrlEv1 -> SetMinimum(0);

    fHistCtrlEv2 = new TH2D("LTPadPlane_CtrlEv2","",7,0,7,1,0,1);
    fHistCtrlEv2 -> SetStats(0);
    fBinCtrlZZZZZZZ = fHistCtrlEv2 -> GetBin(1,1);
    fBinCtrlEngyMax = fHistCtrlEv2 -> GetBin(2,1);
    fBinCtrlAcmltEv = fHistCtrlEv2 -> GetBin(3,1);
    fBinCtrlAcmltCh = fHistCtrlEv2 -> GetBin(4,1);
    fBinCtrlFitChan = fHistCtrlEv2 -> GetBin(5,1);
    fBinCtrlNEEL500 = fHistCtrlEv2 -> GetBin(6,1);
    fBinCtrlNEEL203 = fHistCtrlEv2 -> GetBin(7,1);
    fHistCtrlEv2 -> GetXaxis() -> SetTickSize(0);
    fHistCtrlEv2 -> GetYaxis() -> SetTickSize(0);
    fHistCtrlEv2 -> GetYaxis() -> SetBinLabel(1,"");
    fHistCtrlEv2 -> GetXaxis() -> SetLabelSize(ctrlLabelSize);
    fHistCtrlEv2 -> GetXaxis() -> SetBinLabel(0,"");
    fHistCtrlEv2 -> GetXaxis() -> SetBinLabel(2,"E_{max}");  // fBinCtrlEngyMax
    fHistCtrlEv2 -> GetXaxis() -> SetBinLabel(3,"++Event");  // fBinCtrlAcmltEv
    fHistCtrlEv2 -> GetXaxis() -> SetBinLabel(4,"++Channel");// fBinCtrlAcmltCh
    fHistCtrlEv2 -> GetXaxis() -> SetBinLabel(5,"Fit Ch.");  // fBinCtrlFitChan
    fHistCtrlEv2 -> GetXaxis() -> SetBinLabel(6,"@E>=500");  // fBinCtrlNEEL500
    fHistCtrlEv2 -> GetXaxis() -> SetBinLabel(7,"@E>=2000"); // fBinCtrlNEEL203
    fHistCtrlEv2 -> SetBinContent(fBinCtrlZZZZZZZ, 0);
    fHistCtrlEv2 -> SetBinContent(fBinCtrlEngyMax, 0);
    fHistCtrlEv2 -> SetBinContent(fBinCtrlAcmltEv, 0);
    fHistCtrlEv2 -> SetBinContent(fBinCtrlAcmltCh, 0);
    fHistCtrlEv2 -> SetBinContent(fBinCtrlFitChan, 0);
    fHistCtrlEv2 -> SetBinContent(fBinCtrlNEEL500, 500);
    fHistCtrlEv2 -> SetBinContent(fBinCtrlNEEL203, 2000);
    fHistCtrlEv2 -> SetMarkerSize(binTextSize);
    fHistCtrlEv2 -> SetMinimum(0);
    fHistCtrlEv2 -> SetMaximum(3000);
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
        fPadSection = new TPad("pad_2devent_zoom","",0,233./700,0.5,1);
        fPadSection -> SetMargin(0.12,0.15,0.1,0.1);
        fPadSection -> SetNumber(1);
        //fPadSection -> AddExec("ex", "LTPadPlane::MouseClickEvent2()");
        fPadSection -> Draw();
        fPadChannel = new TPad("pad_channel","",0,0,0.5,233./700);
        fPadChannel -> SetMargin(0.12,0.15,0.20,0.12);
        fPadChannel -> SetNumber(2);
        fPadChannel -> Draw();
        fPad2DEvent = new TPad("pad_2devent_full","",0.5,233./700,1,1);
        fPad2DEvent -> SetMargin(0.12,0.15,0.1,0.1);
        //fPad2DEvent -> AddExec("ex", "LTPadPlane::MouseClickEvent1()");
        fPad2DEvent -> SetNumber(3);
        fPad2DEvent -> Draw();
        double yCtrl1 = 0;
        double yCtrl2 = 230./700;
        double y1 = 0;
        double y2 = y1 + 0.5*(yCtrl2-yCtrl1);
        fPadCtrlEv1 = new TPad("pad_control1","",0.5,y1,1,y2);
        fPadCtrlEv1 -> SetMargin(0.02,0.02,0.30,0.02);
        fPadCtrlEv1 -> SetNumber(4);
        fPadCtrlEv1 -> Draw();
        y1 = y2;
        y2 = y1 + 0.5*(yCtrl2-yCtrl1);
        fPadCtrlEv2 = new TPad("pad_control2","",0.5,y1,1,y2);
        fPadCtrlEv2 -> SetMargin(0.02,0.02,0.30,0.02);
        fPadCtrlEv2 -> SetNumber(4);
        fPadCtrlEv2 -> Draw();

        AddInteractivePad(fPad2DEvent);
        AddInteractivePad(fPadSection);
        AddInteractivePad(fPadCtrlEv1);
        AddInteractivePad(fPadCtrlEv2);

        fCanvas -> Modified();
        fCanvas -> Update();
    }

    return fCanvas;
}
/*
{
    if (fCanvas==nullptr)
    {
        fCanvas = LKWindowManager::GetWindowManager() -> CanvasResize("TTMicromegas",1100,700,0.9);
        fPad2DEvent = new TPad("pad_2DEvent","",0,230./700,0.5,1);
        fPad2DEvent -> SetMargin(0.12,0.15,0.1,0.1);
        fPad2DEvent -> SetNumber(1);
        fPad2DEvent -> Draw();
        fPadChannel = new TPad("pad_channel","",0,0,0.5,230./700);
        fPadChannel -> SetMargin(0.12,0.05,0.20,0.12);
        fPadChannel -> SetNumber(2);
        fPadChannel -> Draw();
        fPad3DEvent = new TPad("pad_3DEvent","",0.5,230./700,1,1);
        fPad3DEvent -> SetMargin(0.12,0.15,0.1,0.1);
        fPad3DEvent -> SetNumber(3);
        fPad3DEvent -> Draw();

        double yCtrl1 = 0;
        double yCtrl2 = 230./700;
        double y1 = 0;
        double y2 = y1 + 0.5*(yCtrl2-yCtrl1);
        fPadCtrlEv1 = new TPad("pad_control","",0.5,y1,1,y2);
        //fPadCtrlEv1 -> SetMargin(0.12,0.05,0.20,0.12);
        fPadCtrlEv1 -> SetMargin(0.02,0.02,0.30,0.02);
        fPadCtrlEv1 -> SetNumber(4);
        fPadCtrlEv1 -> Draw();
        y1 = y2;
        y2 = y1 + 0.5*(yCtrl2-yCtrl1);
        fPadCtrlEv2 = new TPad("pad_control","",0.5,y1,1,y2);
        //fPadCtrlEv2 -> SetMargin(0.12,0.05,0.20,0.12);
        fPadCtrlEv2 -> SetMargin(0.02,0.02,0.30,0.02);
        fPadCtrlEv2 -> SetNumber(4);
        fPadCtrlEv2 -> Draw();

        fCanvas -> Modified();
        fCanvas -> Update();

        AddInteractivePad(fPad2DEvent);
        AddInteractivePad(fPadChannel);
        AddInteractivePad(fPadCtrlEv1);
        AddInteractivePad(fPadCtrlEv2);
    }

    return fCanvas;
}
*/

LKPad *LTPadPlane::NewPad(Int_t section, Int_t layer, Int_t row)
{
    auto pad = new LKPad();
    pad -> SetSectionLayerRow(section, layer, row);
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

void LTPadPlane::Clear(Option_t *)
{
    LKDetectorPlane::Clear();
    if (fAccumulateEvents==0) {
        fHistPadPlane -> Reset("ICES");
        fHistTopView -> Reset("ICES");
        fHistSideView -> Reset("ICES");
        fHistPadPlaneSection[0] -> Reset("ICES");
        fHistPadPlaneSection[1] -> Reset("ICES");
        fHistPadPlaneSection[2] -> Reset("ICES");
        fHistPadPlaneSection[3] -> Reset("ICES");
        fHistPadPlaneSection[4] -> Reset("ICES");
        fHistPadPlaneSection[5] -> Reset("ICES");
        fHistPadPlaneSection[6] -> Reset("ICES");
        fHistPadPlaneSection[7] -> Reset("ICES");
    }
    fSelectedPad = nullptr;
}

void LTPadPlane::FillDataToHist()
{
    Clear();

    if (fPar -> CheckPar("eve/planeFillType")) {
        fFillType = fPar -> GetParString("eve/planeFillType");
        if (fFillType!="Hit" && fFillType!="Raw") {
            lk_warning << "eve/planeFillType must be Raw or Hit" << endl;
            if (fRawDataArray!=nullptr) fFillType = "Raw";
            else if (fHitArray!=nullptr) fFillType = "Hit";
        }
    }

    if (fRawDataArray!=nullptr)
    {
        fHistWaveform -> Reset("ICES");

        auto numChannels = fRawDataArray -> GetEntries();
        lk_info << "# of channels: " << numChannels << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (GETChannel*) fRawDataArray -> At(iChannel);
            auto padID = channel -> GetChan2();
            if (padID<0) {
                auto cobo = channel -> GetCobo();
                auto asad = channel -> GetAsad();
                auto aget = channel -> GetAget();
                auto chan = channel -> GetChan();
                padID = FindPadID(cobo,asad,aget,chan);
            }
            if (padID<0) {
                auto cobo = channel -> GetCobo();
                auto asad = channel -> GetAsad();
                auto aget = channel -> GetAget();
                auto chan = channel -> GetChan();
                lk_warning << "This pad(CAAC) is not mapped!: " << cobo << " " << asad << " " << aget << " " << chan << endl;
                continue;
            }
            auto buffer = channel -> GetWaveformY();
            for (auto tb=0; tb<fNumTbs; ++tb)
                fHistWaveform -> Fill(tb,buffer[tb]);
            auto pad = LKDetectorPlane::GetPad(padID);
            pad -> SetBuffer(buffer);
        }
    }

    if (fHitArray!=nullptr)
    {
        //fHitArray -> Print();
        //fHistWaveform -> Reset("ICES");

        auto numHits = fHitArray -> GetEntries();
        lk_info << "# of hits: " << numHits << endl;
        for (auto iHit=0; iHit<numHits; ++iHit)
        {
            auto hit = (LKHit*) fHitArray -> At(iHit);
            auto padID = hit -> GetChannelID();
            if (padID<0) {
                lk_warning << "padID is <0!" << endl;
                continue;
            }
            auto pad = LKDetectorPlane::GetPad(padID);
            pad -> AddHit(hit);
        }
    }

    if (fFillType=="Raw")
    {
        if (fRawDataArray==nullptr) {
            lk_error << "Fill type is " << fFillType << " but buffer array is null!" << endl;
            return;
        }

        fHistPadPlane -> Reset("ICES");
        fHistSideView -> Reset("ICES");
        fHistTopView  -> Reset("ICES");
        for (auto section=0; section<8; ++section)
            fHistPadPlaneSection[section] -> Reset("ICES");

        auto numChannels = fRawDataArray -> GetEntries();
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (GETChannel*) fRawDataArray -> At(iChannel);
            auto padID = channel -> GetChan2();
            if (padID<0) {
                auto cobo = channel -> GetCobo();
                auto asad = channel -> GetAsad();
                auto aget = channel -> GetAget();
                auto chan = channel -> GetChan();
                padID = FindPadID(cobo,asad,aget,chan);
            }
            if (padID<0) {
                //auto cobo = channel -> GetCobo();
                //auto asad = channel -> GetAsad();
                //auto aget = channel -> GetAget();
                //auto chan = channel -> GetChan();
                //lk_warning << "This pad(CAAC) is not mapped!: " << cobo << " " << asad << " " << aget << " " << chan << endl;
                continue;
            }
            auto charge = channel -> GetEnergy();
            const int numCollectBin = 10;
            int max[numCollectBin] = {0};
            int bin[numCollectBin] = {0};
            if (charge<=0) {
                charge = 0;
                auto buffer = channel -> GetWaveformY();
                for (auto tb=0; tb<fNumTbs; ++tb) {
                    auto value = buffer[tb];
                    charge += value;
                         if (max[0]<value) { max[0] = value; bin[0] = tb; }
                    else if (max[1]<value) { max[1] = value; bin[1] = tb; }
                    else if (max[2]<value) { max[2] = value; bin[2] = tb; }
                    else if (max[3]<value) { max[3] = value; bin[3] = tb; }
                    else if (max[4]<value) { max[4] = value; bin[4] = tb; }
                    else if (max[5]<value) { max[5] = value; bin[5] = tb; }
                    else if (max[6]<value) { max[6] = value; bin[6] = tb; }
                    else if (max[7]<value) { max[7] = value; bin[7] = tb; }
                    else if (max[8]<value) { max[8] = value; bin[8] = tb; }
                    else if (max[9]<value) { max[9] = value; bin[9] = tb; }
                }
            }
            if (charge<0)
                continue;
            auto pad = LKDetectorPlane::GetPad(padID);
            auto x = pad -> GetX();
            auto y = pad -> GetY();
            auto section = pad -> GetSection();
            auto hbin = padID+5;
            //fHistPadPlane -> //SetBinContent(hbin,charge);
            fHistPadPlane -> Fill(x,y,charge);
            fHistPadPlaneSection[section] -> Fill(x,y,charge);
            for (auto icb=0; icb<numCollectBin; ++icb) {
                fHistTopView ->  Fill(bin[icb],x,max[icb]);
                fHistSideView -> Fill(bin[icb],y,max[icb]);
            }
        }
    }

    else if (fFillType=="Hit")
    {
        if (fHitArray==nullptr) {
            lk_error << "Fill type is " << fFillType << " but hit array is null!" << endl;
            return;
        }
        auto numHits = fHitArray -> GetEntries();
        for (auto iHit=0; iHit<numHits; ++iHit)
        {
            auto hit = (LKHit*) fHitArray -> At(iHit);
            auto padID = hit -> GetChannelID();
            if (padID<0) {
                lk_warning << "padID is <0!" << endl;
                continue;
            }
            auto pad = LKDetectorPlane::GetPad(padID);
            auto charge = hit -> GetCharge();
            auto x = hit -> GetX();
            auto y = hit -> GetY();
            auto z = hit -> GetZ();
            auto section = pad -> GetSection();
            auto hbin = padID+5;
            //fHistPadPlane -> //SetBinContent(hbin,charge);
            fHistPadPlane -> Fill(x,y,charge);
            fHistPadPlaneSection[section] -> Fill(x,y,charge);
            fHistTopView -> Fill(z,x,charge);
            fHistSideView -> Fill(z,y,charge);
        }
    }

    if (fTrackArray!=nullptr)
    {
        auto numTracks = fTrackArray -> GetEntries();
        lk_info << "# of tracks: " << numTracks << endl;
    }
}

void LTPadPlane::Draw(Option_t *option)
{
    if (TString(option)!="r") {
        CreateHistograms();
        FillDataToHist();
    }

    GetCanvas();

    UpdateChannel();
    UpdateSection();
    Update2DEvent();
    UpdateCtrlEv1();
    UpdateCtrlEv2();
}

void LTPadPlane::ZoomInWindow(Int_t bin, Double_t x, Double_t y)
{
    auto cvs = GetCanvas();
    if (bin==fBinNumberTopView)
        fZoomTopSide = 1;//true;
    else if (bin==fBinNumberSideView)
        fZoomTopSide = 2;//false;
    //else if (bin==fBinZoomZoomButton)
        //fZoomZoomPressed = true;
    else {
        fZoomTopSide = 0;//false;
        fSelectedSection = FindSection(x,y);
    }
    LTPadPlane::Draw("r");
}

void LTPadPlane::SelectAndDrawChannel(Double_t x, Double_t y)
{
    auto padID = fHistPadPlane -> FindBin(x,y) - 1;
    fSelectedPad = LKDetectorPlane::GetPad(padID);
    fSelectedPad -> Print();
    if (fSelectedPad==nullptr) {
        lk_error << endl;
        return;
    }

    double x0,x1,x2,z0,z1,z2;
    fGraphPadBoundary -> Set(0);
    for (auto i=0; i<10; ++i)
        fGraphPadBoundaryNb[i] -> Set(0);

    { // draw pad boundary
        if (0) { // draw neighbor pads
            auto nbPadArray = fSelectedPad -> GetNeighborPadArray();
            auto numNbPads = nbPadArray -> size();
            for (auto iPad=0; iPad<numNbPads; ++iPad)
            {
                auto nbPad = nbPadArray -> at(iPad);
                auto graphNb = fGraphPadBoundaryNb[iPad];
                auto corners = nbPad -> GetPadCorners();
                Int_t numCorners = corners->size();
                for (auto iCorner=0; iCorner<numCorners+1; ++iCorner) {
                    auto iAt = iCorner;
                    if (iCorner==numCorners) iAt = 0;
                    TVector2 corner = corners->at(iAt);
                    graphNb -> SetPoint(graphNb->GetN(),corner.X(),corner.Y());
                }
                GetCanvas() -> cd(1);
                graphNb -> Draw("samel");
            }
        }

        auto corners = fSelectedPad -> GetPadCorners();
        Int_t numCorners = corners->size();
        for (auto iCorner=0; iCorner<numCorners+1; ++iCorner)
        {
            auto iAt = iCorner;
            if (iCorner==numCorners) iAt = 0;
            TVector2 corner = corners->at(iAt);
            fGraphPadBoundary -> SetPoint(fGraphPadBoundary->GetN(),corner.X(),corner.Y());
        }
        GetCanvas() -> cd(1);
        fGraphPadBoundary -> Draw("samel");
    }

    GetCanvas() -> cd(2);
    auto buffer = fSelectedPad -> GetBuffer();
    TH1D* histChannel = fHistChannel1;
    TH1D* histChannelPre = fHistChannel2;
    if (!fUseChannel1) {
        histChannel = fHistChannel2;
        histChannelPre = fHistChannel1;
    }
    histChannel -> SetLineColor(602);
    histChannel -> SetLineStyle(1);
    histChannelPre -> SetLineColor(kGray+1);
    histChannelPre -> SetLineStyle(2);
    histChannel -> Reset();
    for (auto i=0; i<fNumTbs; ++i)
        histChannel -> SetBinContent(i+1,buffer[i]);
    histChannel -> Draw();
    histChannelPre -> Draw("same");
    fUseChannel1 = !fUseChannel1;

    auto numHits = fSelectedPad -> GetNumHits();
    for (auto iHit=0; iHit<numHits; ++iHit)
    {
        auto hit = fSelectedPad -> GetHit(iHit);
        hit -> Print();
        auto pulse = fChannelAnalyzer -> GetPulse();
        auto graph = pulse -> GetPulseGraph(hit->GetTb(), hit->GetCharge(), hit->GetPedestal());
        graph -> Draw("samelx");
    }
    histChannel -> SetTitle(Form("PadID=%d, position=(%.2f, %.2f), #Hits=%d",padID,fSelectedPad->GetI(),fSelectedPad->GetJ(),numHits));
}

LKChannelAnalyzer* LTPadPlane::GetChannelAnalyzer(int)
{
    if (fChannelAnalyzer==nullptr)
    {
        if (fPar->CheckPar("LTPadPlane/pulseFile")==false)
            fPar -> AddLine("LTPadPlane/pulseFile {lilak_common}/pulseReference.root");
        TString pulseFileName = fPar -> GetParString("LTPadPlane/pulseFile");
        fChannelAnalyzer = new LKChannelAnalyzer();
        fChannelAnalyzer -> SetPulse(pulseFileName);
        fChannelAnalyzer -> Print();
    }
    return fChannelAnalyzer;
}

void LTPadPlane::ExecMouseClickEventOnPad(TVirtualPad *pad, double xOnClick, double yOnClick)
{
    if (pad==fPad2DEvent) ClickedPad2DEvent(xOnClick, yOnClick);
    if (pad==fPadSection) ClickedPadSection(xOnClick, yOnClick);
    if (pad==fPadCtrlEv1) ClickedPadCtrlEv1(xOnClick, yOnClick);
    if (pad==fPadCtrlEv2) ClickedPadCtrlEv2(xOnClick, yOnClick);
}

void LTPadPlane::ClickedPad2DEvent(double xOnClick, double yOnClick)
{
    int selectedBin = fHistPadPlane -> FindBin(xOnClick, yOnClick);
    ZoomInWindow(selectedBin, xOnClick, yOnClick);
}

void LTPadPlane::ClickedPadSection(double xOnClick, double yOnClick)
{
    SelectAndDrawChannel(xOnClick, yOnClick);
}

void LTPadPlane::ClickedPadCtrlEv1(double xOnClick, double yOnClick)
{
    if (fHistCtrlEv1==nullptr)
        return;

    if (fRun==nullptr)
        return;

    int selectedBin = fHistCtrlEv1 -> FindBin(xOnClick, yOnClick);

    auto currentEventID = fRun -> GetCurrentEventID();
    auto lastEventID = fRun -> GetNumEvents() - 1;

    if (fAccumulateEvents>0)
    {
        if (selectedBin==fBinCtrlFrst) { lk_info << "(First event) option is not available in acuumulate-event-mode" << endl; return; }
        if (selectedBin==fBinCtrlPr50) { lk_info << "(Event -50)   option is not available in acuumulate-event-mode" << endl; return; }
        if (selectedBin==fBinCtrlPrev) { lk_info << "(Prev. event  option is not available in acuumulate-event-mode" << endl; return; }
        if (selectedBin==fBinCtrlCurr) { return; }
        if (selectedBin==fBinCtrlNext) {
            lk_info << "Next event"  << endl;
            fRun -> ExecuteNextEvent();
            fAccumulateEvent2 = fRun -> GetCurrentEventID();
            ++fAccumulateEvents;
        }
        if (selectedBin==fBinCtrlNe50 || selectedBin==fBinCtrlLast)
        {
            Long64_t testEventTo = currentEventID + 50;
            if ((selectedBin==fBinCtrlNe50 && (currentEventID+50>lastEventID)) || (selectedBin==fBinCtrlLast))
                testEventTo = lastEventID;
            lk_info << "Accumulating events: " << currentEventID+1 << " - " << testEventTo << " (" << testEventTo-currentEventID << ")" << endl;
            for (Long64_t eventID=currentEventID+1; eventID<=testEventTo; ++eventID) {
                fRun -> GetEvent(eventID);
                SetDataFromBranch();
                FillDataToHist();
                fAccumulateEvent2 = fRun -> GetCurrentEventID();
                ++fAccumulateEvents;
            }
            fHistCtrlEv2 -> SetBinContent(fBinCtrlAcmltEv, fAccumulateEvents);
        }
    }
    else
    {
        if (selectedBin==fBinCtrlFrst) { lk_info << "First event" << endl; fRun -> ExecuteFirstEvent(); }
        if (selectedBin==fBinCtrlPr50) { lk_info << "Event -50"   << endl; fRun -> ExecuteEvent((currentEventID-50<0?0:currentEventID-50)); }
        if (selectedBin==fBinCtrlPrev) { lk_info << "Prev. event" << endl; fRun -> ExecutePreviousEvent(); }
        if (selectedBin==fBinCtrlCurr) { return; }
        if (selectedBin==fBinCtrlNext) { lk_info << "Next event"  << endl; fRun -> ExecuteNextEvent(); }
        if (selectedBin==fBinCtrlNe50) { lk_info << "Event +50"   << endl; fRun -> ExecuteEvent((currentEventID+50>lastEventID?lastEventID:currentEventID+50)); }
        if (selectedBin==fBinCtrlLast) { lk_info << "Last event"  << endl; fRun -> ExecuteLastEvent(); }
    }

    Draw();
}

void LTPadPlane::ClickedPadCtrlEv2(double xOnClick, double yOnClick)
{
    if (fHistCtrlEv2==nullptr)
        return;

    int selectedBin = fHistCtrlEv2 -> FindBin(xOnClick, yOnClick);

    Long64_t currentEventID;
    Long64_t lastEventID;

    if (selectedBin==fBinCtrlZZZZZZZ) { return; }
    if (selectedBin==fBinCtrlNEEL500 || selectedBin==fBinCtrlNEEL203)
    {
        if (fRun==nullptr)
            return;

        currentEventID = fRun -> GetCurrentEventID();
        lastEventID = fRun -> GetNumEvents() - 1;

        double energyCut = 500;
        if (selectedBin==fBinCtrlNEEL500) energyCut = 500;
        else if (selectedBin==fBinCtrlNEEL203) energyCut = 2000;

        if (fRawDataArray==nullptr)
            return;

        auto testEventID = currentEventID;
        while (currentEventID<=lastEventID+1)
        {
            testEventID++;
            lk_info << "Testing " << testEventID << endl;

            //fRun -> ExecuteNextEvent();
            fRun -> GetEvent(testEventID);

            double maxEnergy = 0;
            auto numChannels = fRawDataArray -> GetEntries();
            for (auto iRawData=0; iRawData<numChannels; ++iRawData)
            {
                auto channel = (GETChannel*) fRawDataArray -> At(iRawData);
                if (channel->GetEnergy()>maxEnergy) {
                    maxEnergy = channel->GetEnergy();
                    break;
                }
            }
            if (maxEnergy>energyCut)
                break;

            if (testEventID>=lastEventID)
                break;
        }
        if (testEventID==lastEventID) {
            lk_error << "No event with energy " << energyCut << endl;
            return;
        }

        fRun -> ExecuteEvent(testEventID);
        lk_info << "Event with energy " << energyCut << " : " << currentEventID << endl;
    }
    if (selectedBin==fBinCtrlEngyMax) {
        if (fEnergyMaxMode==0) {
            fEnergyMaxMode = 1;
            fHistCtrlEv2 -> SetBinContent(fBinCtrlEngyMax, 2500);
            lk_info << "Set energy range automatic" << endl;
        }
        else if (fEnergyMaxMode==1) {
            fEnergyMaxMode = 2;
            fHistCtrlEv2 -> SetBinContent(fBinCtrlEngyMax, 4200);
            lk_info << "Set energy range to 2500" << endl;
        }
        else //if (fEnergyMaxMode==2)
        {
            fEnergyMaxMode = 0;
            fHistCtrlEv2 -> SetBinContent(fBinCtrlEngyMax, 0);
            lk_info << "Set energy range to 4200" << endl;
        }
    }
    if (selectedBin==fBinCtrlAcmltEv)
    {
        if (fRun==nullptr)
            return;
        currentEventID = fRun -> GetCurrentEventID();
        lastEventID = fRun -> GetNumEvents() - 1;

        if (fAccumulateEvents>0)
            fAccumulateEvents = 0;
        else {
            fAccumulateEvents = 1;
            fAccumulateEvent1 = currentEventID;
            fAccumulateEvent2 = currentEventID;
        }
        fHistCtrlEv2 -> SetBinContent(fBinCtrlAcmltEv, fAccumulateEvents);
        return;
    }
    if (selectedBin==fBinCtrlAcmltCh)
    {
        if (fAccumulateChannel) {
            fCountChannelGraph = 0;
            fAccumulateChannel = false;
            fHistCtrlEv2 -> SetBinContent(fBinCtrlAcmltCh, 0);
        }
        else {
            fCountChannelGraph = 0;
            fAccumulateChannel = true;
            fHistCtrlEv2 -> SetBinContent(fBinCtrlAcmltCh, 1);
        }
        UpdateChannel();
        return;
    }
    if (selectedBin==fBinCtrlFitChan)
    {
        lk_info << "Fit channel" << endl;
        fFitChannel = true;
        UpdateChannel();
        return;
    }

    Draw();
}

void LTPadPlane::UpdateChannel()
{
    fPadChannel -> cd();
    fHistWaveform -> Draw("colz");
}

void LTPadPlane::UpdateSection()
{
    fPadSection -> cd();

    if (fZoomTopSide==1) {
        fHistTopView -> SetMinimum(fHistZMin);
        fHistTopView -> Draw("colz");
        fGraphSectionBoundary2 -> Set(0);
        fGraphSectionBoundary2 -> SetPoint(0, fXMax-50,  fYMax-50);
        fGraphSectionBoundary2 -> SetPoint(1, fXMax-50,  fYMax-250);
        fGraphSectionBoundary2 -> SetPoint(2, fXMax-250, fYMax-50);
        fGraphSectionBoundary2 -> SetPoint(3, fXMax-50,  fYMax-50);
        if (fTrackArray!=nullptr) {
            LKVector3::Axis axis1 = LKVector3::kZ;
            LKVector3::Axis axis2 = LKVector3::kX;
            auto numTracks = fTrackArray -> GetEntries();
            for (auto iTrack = 0; iTrack < numTracks; ++iTrack) {
                auto track = (LKTracklet *) fTrackArray -> At(iTrack);
                auto graph = track -> TrajectoryOnPlane(axis1,axis2);
                auto graphNew = (TGraphErrors*) graph -> Clone(); //TODO
                graphNew -> Draw("samel");
            }
        }
    }
    else if (fZoomTopSide==2)
    {
        fHistSideView -> SetMinimum(fHistZMin);
        fHistSideView -> Draw("colz");
        fGraphSectionBoundary2 -> Set(0);
        fGraphSectionBoundary2 -> SetPoint(0, fXMin+50,  fYMin+50);
        fGraphSectionBoundary2 -> SetPoint(1, fXMin+50,  fYMin+250);
        fGraphSectionBoundary2 -> SetPoint(2, fXMin+250, fYMin+50);
        fGraphSectionBoundary2 -> SetPoint(3, fXMin+50,  fYMin+50);
        if (fTrackArray!=nullptr) {
            LKVector3::Axis axis1 = LKVector3::kZ;
            LKVector3::Axis axis2 = LKVector3::kY;
            auto numTracks = fTrackArray -> GetEntries();
            for (auto iTrack = 0; iTrack < numTracks; ++iTrack) {
                auto track = (LKTracklet *) fTrackArray -> At(iTrack);
                auto graph = track -> TrajectoryOnPlane(axis1,axis2);
                auto graphNew = (TGraphErrors*) graph -> Clone(); //TODO
                graphNew -> Draw("samel");
            }
        }
    }
    else {
        if (fHistPadPlaneSection[fSelectedSection]->GetEntries()==0)
            fHistPadPlaneSection[fSelectedSection] -> Draw();
        else {
            fFramePadPlaneSection[fSelectedSection] -> Draw();
            fHistPadPlaneSection[fSelectedSection] -> SetMinimum(fHistZMin);
            fHistPadPlaneSection[fSelectedSection] -> Draw("same colz");
        }
        if (fTrackArray!=nullptr) {
            LKVector3::Axis axis1 = LKVector3::kX;
            LKVector3::Axis axis2 = LKVector3::kY;
            auto numTracks = fTrackArray -> GetEntries();
            for (auto iTrack = 0; iTrack < numTracks; ++iTrack) {
                auto track = (LKTracklet *) fTrackArray -> At(iTrack);
                auto graph = track -> TrajectoryOnPlane(axis1,axis2);
                auto graphNew = (TGraphErrors*) graph -> Clone(); //TODO
                graphNew -> Draw("samel");
            }
        }
    }
}

void LTPadPlane::Update2DEvent()
{
    fPad2DEvent -> cd();

    for (TH1* hist : {(TH1*) fHistPadPlane, (TH1*) fHistSideView, (TH1*) fHistTopView,
            (TH1*) fHistPadPlaneSection[0], (TH1*) fHistPadPlaneSection[1], (TH1*) fHistPadPlaneSection[2], (TH1*) fHistPadPlaneSection[3],
            (TH1*) fHistPadPlaneSection[4], (TH1*) fHistPadPlaneSection[5], (TH1*) fHistPadPlaneSection[6], (TH1*) fHistPadPlaneSection[7]}
        )
    {
        if      (fEnergyMaxMode==0) { hist -> SetMinimum(-1111);   hist -> SetMaximum(-1111); }
        else if (fEnergyMaxMode==1) { hist -> SetMinimum(0);       hist -> SetMaximum(2500); }
        else if (fEnergyMaxMode==2) { hist -> SetMinimum(0);       hist -> SetMaximum(4200); }
        else
        {
            hist -> SetMinimum(0);
            hist -> SetMaximum(fEnergyMaxMode);
            if (fEnergyMaxMode>100)
                gStyle -> SetNumberContours(100);
            else
                gStyle -> SetNumberContours(fEnergyMaxMode);
        }
    }

    if (fHistPadPlane->GetEntries()==0)
        fFramePadPlane -> Draw();
    else {
        fFramePadPlane -> Draw();
        fHistPadPlane -> SetMinimum(fHistZMin);
        fHistPadPlane -> Draw("same colz");
    }
    fGraphSectionBoundary2 -> Draw("samel");

    Double_t t1, t2, xmin, xmax, ymin, ymax;
    GetSectionParameters(fSelectedSection,t1,t2,xmin,xmax,ymin,ymax);
    fGraphSectionBoundary1 -> Set(0);
    fGraphSectionBoundary1 -> SetPoint(0,fPosSectionCorner[fSelectedSection][0].X(), fPosSectionCorner[fSelectedSection][0].Y());
    fGraphSectionBoundary1 -> SetPoint(1,fPosSectionCorner[fSelectedSection][1].X(), fPosSectionCorner[fSelectedSection][1].Y());
    fGraphSectionBoundary1 -> SetPoint(2,fPosSectionCorner[fSelectedSection][2].X(), fPosSectionCorner[fSelectedSection][2].Y());
    fGraphSectionBoundary1 -> SetPoint(3,fPosSectionCorner[fSelectedSection][3].X(), fPosSectionCorner[fSelectedSection][3].Y());
    fGraphSectionBoundary1 -> SetPoint(4,fPosSectionCorner[fSelectedSection][0].X(), fPosSectionCorner[fSelectedSection][0].Y());
    fGraphSectionBoundary1 -> Draw("samel");
    if (fTrackArray!=nullptr) {
        LKVector3::Axis axis1 = LKVector3::kX;
        LKVector3::Axis axis2 = LKVector3::kY;
        auto numTracks = fTrackArray -> GetEntries();
        for (auto iTrack = 0; iTrack < numTracks; ++iTrack) {
            auto track = (LKTracklet *) fTrackArray -> At(iTrack);
            auto graph = track -> TrajectoryOnPlane(axis1,axis2);
            auto graphNew = (TGraphErrors*) graph -> Clone(); //TODO
            graphNew -> Draw("samel");
        }
    }
}

void LTPadPlane::UpdateCtrlEv1()
{
    if (fHistCtrlEv1==nullptr)
        return;
    fPadCtrlEv1 -> cd();
    fPadCtrlEv1 -> SetGrid();
    if (fRun!=nullptr) {
        auto currentEventID = fRun -> GetCurrentEventID();
        auto lastEventID = fRun -> GetNumEvents() - 1;
        fHistCtrlEv1 -> SetBinContent(fBinCtrlPr50, (currentEventID-50<0?0:currentEventID-50));
        fHistCtrlEv1 -> SetBinContent(fBinCtrlPrev, (currentEventID==0?0:currentEventID-1));
        fHistCtrlEv1 -> SetBinContent(fBinCtrlCurr, currentEventID);
        fHistCtrlEv1 -> SetBinContent(fBinCtrlNext, (currentEventID==lastEventID?lastEventID:currentEventID+1));
        fHistCtrlEv1 -> SetBinContent(fBinCtrlNe50, (currentEventID+50>lastEventID?lastEventID:currentEventID+50));
        fHistCtrlEv1 -> Draw("col text");
    }
    else
        fHistCtrlEv1 -> Draw("text");
}

void LTPadPlane::UpdateCtrlEv2()
{
    if (fHistCtrlEv2==nullptr)
        return;
    fPadCtrlEv2 -> cd();
    fPadCtrlEv2 -> SetGrid();
    if (fRun!=nullptr)
        fHistCtrlEv2 -> Draw("col text");
    else
        fHistCtrlEv2 -> Draw("text");
}
