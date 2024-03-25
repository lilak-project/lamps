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

LTPadPlane::LTPadPlane()
    :LKEvePlane("LTPadPlane", "pad plane with rectangular pads following the circle line for LAMPS-TPC")
{
    fName = "LTPadPlane";
    fAxis1 = LKVector3::kX;
    fAxis2 = LKVector3::kY;
    fAxis3 = LKVector3::kZ;
    fAxisDrift = LKVector3::kZ;
    fChannelAnalyzer = nullptr;

    fDXCanvas = 1200;
    fDYCanvas = 800;
    fYCCanvas = 0.32;

    fPaletteNumber = 1;
    fEnergyMin = 1;
    fEnergyMax = 0;
}

void LTPadPlane::Clear(Option_t *)
{
    LKEvePlane::Clear();
}

void LTPadPlane::Print(Option_t *option) const
{
    lk_info << endl;
}

bool LTPadPlane::Init()
{
    LKEvePlane::Init();

    if (fPar -> CheckPar("lamps/rMinTPC"))      fRMin         = fPar -> GetParDouble("lamps/rMinTPC");
    if (fPar -> CheckPar("lamps/rMaxTPC"))      fRMax         = fPar -> GetParDouble("lamps/rMaxTPC");
    if (fPar -> CheckPar("lamps/padGap"))       fPadGap       = fPar -> GetParDouble("lamps/padGap");
    if (fPar -> CheckPar("lamps/rTopCut"))      fRTopCut      = fPar -> GetParDouble("lamps/rTopCut");
    if (fPar -> CheckPar("lamps/padWidth"))     fPadWidth     = fPar -> GetParDouble("lamps/padWidth");
    if (fPar -> CheckPar("lamps/padHeight"))    fPadHeight    = fPar -> GetParDouble("lamps/padHeight");
    if (fPar -> CheckPar("lamps/radiusLayer0")) fRadiusLayer0 = fPar -> GetParDouble("lamps/radiusLayer0");
    if (fPar -> CheckPar("lamps/numLayers"))    fNumLayers    = fPar -> GetParInt   ("lamps/numLayers");

    fXSpacing = fPadGap + fPadWidth;
    fRSpacing = fPadGap + fPadHeight;

    fTanPi1o8 = TMath::Tan(TMath::Pi()*1./8.);
    fTanPi3o8 = TMath::Tan(TMath::Pi()*3./8.);
    fTanPi5o8 = TMath::Tan(TMath::Pi()*5./8.);
    fTanPi7o8 = TMath::Tan(TMath::Pi()*7./8.);
    for (int i = 0; i < 8; i++) {
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
        double r1 = fRMin-10;
        double r2 = 560;
        double x1 = -r2/sqrt(1+fTanPi3o8*fTanPi3o8);
        double x2 = -r1/sqrt(1+fTanPi3o8*fTanPi3o8);
        double x3 = -x2;
        double x4 = -x1;
        double y1 = -fTanPi3o8*x1;
        double y2 = -fTanPi3o8*x2;
        double y3 = +fTanPi3o8*x3;
        double y4 = +fTanPi3o8*x4;
        TVector2 sectionCorner1(x1, y1);
        TVector2 sectionCorner2(x2, y2);
        TVector2 sectionCorner3(x3, y3);
        TVector2 sectionCorner4(x4, y4);

        double phiSection = section * TMath::Pi()/4.;
        fPosSectionCorner[section][0] = sectionCorner1.Rotate(phiSection);
        fPosSectionCorner[section][1] = sectionCorner2.Rotate(phiSection);
        fPosSectionCorner[section][2] = sectionCorner3.Rotate(phiSection);
        fPosSectionCorner[section][3] = sectionCorner4.Rotate(phiSection);
    }

    double xCorner[5] = {0};
    double yCorner[5] = {0};

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
            double yCornerTemp;
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
                double phiSection = section * TMath::Pi()/4.;

                TVector2 point;
                for (int iRL : {0,1})
                {
                    int signRL = (iRL==0?1:-1);
                    int pm = (iRL==0?1:0);

                    auto pad = new LKPad();
                    auto padID = fChannelArray -> GetEntriesFast();
                    pad -> SetSectionLayerRow(section, layer, signRL*row);
                    pad -> SetPadID(padID);
                    fChannelArray -> Add(pad);

                    //auto padID = pad -> GetPadID();
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
                auto padL = LKEvePlane::GetPad(section,layer,-row);
                auto padR = LKEvePlane::GetPad(section,layer,+row);

                // pad below
                if (layer>0)
                {
                    auto rowNb = row + fNumSkippedHalfRows[layer] - fNumSkippedHalfRows[layer-1];

                    auto rowL1 = -rowNb - 1;
                    if (abs(rowL1)<=fNumHalfRowsInLayer[layer-1]) {
                        auto padBelowL1 = LKEvePlane::GetPad(section,layer-1,rowL1);
                        if (padBelowL1!=nullptr) SetNeighborPads(padL,padBelowL1);
                    }

                    auto rowL2 = -rowNb + 1; if (rowL2==0) rowL2 += 1;
                    auto padBelowL2 = LKEvePlane::GetPad(section,layer-1,rowL2);
                    if (padBelowL2!=nullptr) SetNeighborPads(padL,padBelowL2);

                    auto rowR1 = rowNb - 1; if (rowR1==0) rowR1 -= 1;
                    auto padBelowR1 = LKEvePlane::GetPad(section,layer-1,rowR1);
                    if (padBelowR1!=nullptr) SetNeighborPads(padR,padBelowR1);

                    auto rowR2 = rowNb + 1;
                    if (abs(rowR2)<=fNumHalfRowsInLayer[layer-1]) {
                        auto padBelowR2 = LKEvePlane::GetPad(section,layer-1,rowR2);
                        if (padBelowR2!=nullptr) SetNeighborPads(padR,padBelowR2);
                    }

                    if (rowNb<=0||rowNb>fNumHalfRowsInLayer[layer-1])
                        continue;

                    auto padBelowL = LKEvePlane::GetPad(section,layer-1,-rowNb);
                    SetNeighborPads(padL,padBelowL);

                    auto padBelowR = LKEvePlane::GetPad(section,layer-1,+rowNb);
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

    fMapCAACToPadID = new int***[fMaxCobo];
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

    if (fPar -> CheckPar("lamps/position_map"))
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
        auto mapFileName = fPar -> GetParString("lamps/position_map");
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

int LTPadPlane::FindPadID(int section, int layer, int row)
{
    int idLayer = section*fNumPadsDownToLayer[layer];
    if (layer<fNumLayers) idLayer = idLayer + (8-section)*fNumPadsDownToLayer[layer+1];
    int idRow = ((row>0) ? (fNumHalfRowsInLayer[layer] - row) : (fNumHalfRowsInLayer[layer] - row - 1));
    int padID = idLayer + idRow;

    return padID;
}

int LTPadPlane::FindPadID(int cobo, int asad, int aget, int chan)
{
    auto padID = fMapCAACToPadID[cobo][asad][aget][chan];
    return padID;
}

int LTPadPlane::FindPadID(double i, double j)
{
    int section = FindSection(i,j);

    double xRotatedToSec0 =  i*fCosPiNo4[section] + j*fSinPiNo4[section];
    double yRotatedToSec0 = -i*fSinPiNo4[section] + j*fCosPiNo4[section];
    double rRotatedToSec0 =  sqrt(xRotatedToSec0*xRotatedToSec0 + yRotatedToSec0*yRotatedToSec0);

    double rFromSectionBottom = rRotatedToSec0 - fRadiusLayer0 + .5*fPadHeight;
    if (rFromSectionBottom < 0)
        return -1;

    int layer = (int)(rFromSectionBottom/fRSpacing);
    if (layer > fNumLayers)
        return -2;

    int pm = 1;
    if (xRotatedToSec0 < 0) {
        xRotatedToSec0 = -xRotatedToSec0;
        pm = -1;
    }

    if (xRotatedToSec0 < .5*fPadGap)
        return -3;

    double xFromRow0LeftEdge = xRotatedToSec0 - .5*fPadGap;
    int row = (int) (xFromRow0LeftEdge / fXSpacing) + 1;
    if (xFromRow0LeftEdge - (row-1)*fXSpacing > fPadWidth)
        return -4;

    row = row - fNumSkippedHalfRows[layer];

    if (row > fNumHalfRowsInLayer[layer])
        return -5;

    return FindPadID(section,layer,pm*row);
}


double LTPadPlane::PadDisplacement() const
{
    return sqrt(fXSpacing*fXSpacing + fRSpacing*fRSpacing);
}

bool LTPadPlane::IsInBoundary(double i, double j)
{
    double r = TMath::Sqrt(i*i+j*j);
    if (r < fRMin || r > fRMax)
        return false;

    return true;
}

TH2* LTPadPlane::GetHistEventDisplay1(Option_t* option)
{
    if (fHistHeadView==nullptr)
    {
        if (fDefinePositionByPixelIndex)
        {
            fHistHeadView = new TH2D("histLTTV", "LAMPS-TPC Head View;time-bucket;x (mm)",128,0,512,2*fNumLayers,fXMin,fXMax);
            fHistHeadView -> SetStats(0);
            fHistHeadView -> GetXaxis() -> SetTickLength(0.01);
            fHistHeadView -> GetYaxis() -> SetTickLength(0.01);
            fHistHeadView -> SetMinimum(fHistZMin);

            fHistSideView = new TH2D("histLTSV", "LAMPS-TPC Side View;time-bucket;y (mm)",128,0,512,2*fNumLayers,fYMin,fYMax);
            fHistSideView -> SetStats(0);
            fHistSideView -> GetXaxis() -> SetTickLength(0.01);
            fHistSideView -> GetYaxis() -> SetTickLength(0.01);
            fHistSideView -> SetMinimum(fHistZMin);
        }
        else
        {
            fHistHeadView = new TH2D("histLTTV", "LAMPS-TPC Head View;z (mm);x (mm)",fZBins,fZMin,fZMax,2*fNumLayers,fXMin,fXMax);
            fHistHeadView -> SetStats(0);
            fHistHeadView -> GetXaxis() -> SetTickLength(0.01);
            fHistHeadView -> GetYaxis() -> SetTickLength(0.01);
            fHistHeadView -> SetMinimum(fHistZMin);

            fHistSideView = new TH2D("histLTSV", "LAMPS-TPC Side View;z (mm);y (mm)",fZBins,fZMin,fZMax,2*fNumLayers,fYMin,fYMax);
            fHistSideView -> SetStats(0);
            fHistSideView -> GetXaxis() -> SetTickLength(0.01);
            fHistSideView -> GetYaxis() -> SetTickLength(0.01);
            fHistSideView -> SetMinimum(fHistZMin);
        }

        fGraphPadBoundary = new TGraph();
        fGraphPadBoundary -> SetLineColor(kRed);
        fGraphPadBoundary -> SetLineWidth(3);
        for (auto i=0; i<10; ++i) {
            fGraphPadBoundaryNb[i] = new TGraph();
            fGraphPadBoundaryNb[i] -> SetLineColor(kRed);
            fGraphPadBoundaryNb[i] -> SetLineWidth(2);
            fGraphPadBoundaryNb[i] -> SetLineStyle(2);
        }
    }

    return (TH2*) fHistHeadView;
}

TH2* LTPadPlane::GetHistEventDisplay2(Option_t* option)
{
    if (fGridPadPlane==nullptr)
    {
        fGridPadPlane = new TH2Poly("frameLTPP","LAMPS-TPC Pad Plane;x (mm);y (mm)",fXMin,fXMax,fYMin,fYMax);
        fGridPadPlane -> SetStats(0);
        fGridPadPlane -> GetXaxis() -> SetTickLength(0.01);
        fGridPadPlane -> GetYaxis() -> SetTickLength(0.01);

        fHistPadPlane = new TH2Poly("histLTPP", "LAMPS-TPC Pad Plane;x (mm);y (mm)",fXMin,fXMax,fYMin,fYMax);
        fHistPadPlane -> SetStats(0);
        fHistPadPlane -> GetXaxis() -> SetTickLength(0.01);
        fHistPadPlane -> GetYaxis() -> SetTickLength(0.01);
        fHistPadPlane -> SetMinimum(fHistZMin);

        for (auto section=0; section<8; ++section)
        {
            double t1, t2, xmin, xmax, ymin, ymax;
            GetSectionParameters(section,t1,t2,xmin,xmax,ymin,ymax);

            fGridSectView[section] = new TH2Poly(Form("frameLTPP_S%d",section),Form("LAMPS-TPC Pad Plane Section-%d;x (mm);y (mm)",section),xmin,xmax,ymin,ymax);
            fGridSectView[section] -> SetStats(0);
            fGridSectView[section] -> GetXaxis() -> SetTickLength(0.01);
            fGridSectView[section] -> GetYaxis() -> SetTickLength(0.01);
            fGridSectView[section] -> SetMinimum(fHistZMin);

            fHistSectView[section] = new TH2Poly(Form("histLTPP_S%d",section),Form("LAMPS-TPC Pad Plane Section-%d;x (mm);y (mm)",section),xmin,xmax,ymin,ymax);
            fHistSectView[section] -> SetStats(0);
            fHistSectView[section] -> GetXaxis() -> SetTickLength(0.01);
            fHistSectView[section] -> GetYaxis() -> SetTickLength(0.01);
            fHistSectView[section] -> SetMinimum(fHistZMin);
        }

        fGraphSectionBoundary1 = new TGraph();
        fGraphSectionBoundary1 -> SetLineColor(kBlue);
        fGraphSectionBoundary1 -> SetLineWidth(3);

        fGraphSectionBoundary2 = new TGraph();
        fGraphSectionBoundary2 -> SetLineColor(kRed);
        fGraphSectionBoundary2 -> SetLineWidth(3);

        LKPad *pad;
        double xPoints[10] = {0};
        double yPoints[10] = {0};
        TIter iterPads(fChannelArray);
        while ((pad = (LKPad *) iterPads.Next())) 
        {
            auto corners = pad -> GetPadCorners();
            int numCorners = corners->size();
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
            fHistSectView[section] -> AddBin(numCorners+1, xPoints, yPoints);
            fGridSectView[section] -> AddBin(numCorners+1, xPoints, yPoints);
        }

        for (auto section=0; section<8; ++section) {
            xPoints[0] = fPosSectionCorner[section][0].X(); yPoints[0] = fPosSectionCorner[section][0].Y();
            xPoints[1] = fPosSectionCorner[section][1].X(); yPoints[1] = fPosSectionCorner[section][1].Y();
            xPoints[2] = fPosSectionCorner[section][2].X(); yPoints[2] = fPosSectionCorner[section][2].Y();
            xPoints[3] = fPosSectionCorner[section][3].X(); yPoints[3] = fPosSectionCorner[section][3].Y();
            xPoints[4] = fPosSectionCorner[section][0].X(); yPoints[4] = fPosSectionCorner[section][0].Y();
            fGridPadPlane -> AddBin(5, xPoints, yPoints);
        }

        xPoints[0] = fXMax-50;  yPoints[0] = fYMax-50;
        xPoints[1] = fXMax-50;  yPoints[1] = fYMax-250;
        xPoints[2] = fXMax-250; yPoints[2] = fYMax-50;
        xPoints[3] = fXMax-50;  yPoints[3] = fYMax-50;
        fBinNumberHeadView = fGridPadPlane -> AddBin(4, xPoints, yPoints); 

        xPoints[0] = fXMin+50;  yPoints[0] = fYMin+50;
        xPoints[1] = fXMin+50;  yPoints[1] = fYMin+250;
        xPoints[2] = fXMin+250; yPoints[2] = fYMin+50;
        xPoints[3] = fXMin+50;  yPoints[3] = fYMin+50;
        fBinNumberSideView = fGridPadPlane -> AddBin(4, xPoints, yPoints); 

        TIter nextBin(fGridPadPlane -> GetBins());
        while (auto bin = (TH2PolyBin*) nextBin()) {
            auto graphBin = (TGraph *) bin -> GetPolygon();
            graphBin -> SetLineColor(kGray+1);
        }

        for (auto section=0; section<8; ++section) {
            TIter nextBin(fGridSectView[section] -> GetBins());
            while (auto bin = (TH2PolyBin*) nextBin()) {
                auto graphBin = (TGraph *) bin -> GetPolygon();
                graphBin -> SetLineColor(kGray+1);
            }
        }

        fHistEventDisplay2 = fHistPadPlane;
        fHistEventDisplay1 = fHistSectView[fSelectedSection];
    }

    return (TH2*) fHistPadPlane;
}

TH2Poly* LTPadPlane::GetHistSection(int selectSection)
{
    GetHist();
    if (selectSection>=0)
        return fHistSectView[selectSection];
    return (TH2Poly*) nullptr;
}

void LTPadPlane::FillDataToHist(Option_t* option)
{
    if (fAccumulateEvents==0) {
        fHistPadPlane -> Reset("ICES");
        fHistHeadView -> Reset("ICES");
        fHistSideView -> Reset("ICES");
        fHistSectView[0] -> Reset("ICES");
        fHistSectView[1] -> Reset("ICES");
        fHistSectView[2] -> Reset("ICES");
        fHistSectView[3] -> Reset("ICES");
        fHistSectView[4] -> Reset("ICES");
        fHistSectView[5] -> Reset("ICES");
        fHistSectView[6] -> Reset("ICES");
        fHistSectView[7] -> Reset("ICES");
    }

    Long64_t currentEventID = 0;
    if (fRun!=nullptr)
        currentEventID = fRun -> GetCurrentEventID();

    TString optionString(option);
    if (!fFillOptionSelected.IsNull())
        optionString = fFillOptionSelected;
    if (optionString.IsNull())
        optionString = "preview";
    optionString.ToLower();
    lk_info << "Filling " << optionString << " (" << currentEventID << ")" << endl;

    LKPhysicalPad *pad = nullptr;
    TString title;

    if (optionString.Index("caac")==0) {
        title = "caac";
        int maxCAAC = 0;
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            auto caac = pad -> GetCAAC();
            if (caac>maxCAAC) maxCAAC = caac;
            double x = pad -> GetI();
            double y = pad -> GetJ();
            fHistPadPlane -> Fill(x,y,caac);
            fHistSectView[pad->GetSection()] -> Fill(x,y,caac);
        }
        lk_debug << fEnergyMax << endl;
        fEnergyMax = maxCAAC;
    }
    else if (optionString.Index("cobo")==0) {
        fEnergyMax = fMaxCobo;
        title = "cobo";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetCoboID());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetCoboID());
        }
    }
    else if (optionString.Index("asad")==0) {
        fEnergyMax = 4;
        title = "asad";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetAsadID());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetAsadID());
        }
    }
    else if (optionString.Index("aget")==0) {
        fEnergyMax = 4;
        title = "aget";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetAgetID());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetAgetID());
        }
    }
    else if (optionString.Index("chan")==0) {
        fEnergyMax = 70;
        title = "chan";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetChannelID());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetChannelID());
        }
    }

    else if (optionString.Index("section")==0) {
        fEnergyMax = 100;
        title = "section";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetSection());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetChannelID());
        }
    }
    else if (optionString.Index("layer")==0) {
        fEnergyMax = 100;
        title = "layer";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetLayer());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetLayer());
        }
    }
    else if (optionString.Index("row")==0) {
        fEnergyMax = 100;
        title = "row";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetRow());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetRow());
        }
    }

    else if (optionString.Index("padid")==0) {
        fEnergyMax = 100;
        title = "id";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetPadID());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetPadID());
        }
    }
    else if (optionString.Index("nhit")==0) {
        fEnergyMax = 10;
        title = "nhit";
        TIter nextRawData(fChannelArray);
        while ((pad = (LKPhysicalPad *) nextRawData())) {
            fHistPadPlane -> Fill(pad->GetI(),pad->GetJ(),pad->GetNumHits());
            fHistSectView[pad->GetSection()] -> Fill(pad->GetI(),pad->GetJ(),pad->GetNumHits());
        }
    }
    else if (optionString.Index("hit")==0&&fHitArray!=nullptr)
    {
        title = "Hit";
        TIter nextHit(fHitArray);
        LKHit* hit = nullptr;
        while (hit = (LKHit*) nextHit())
        {
            auto x = hit -> GetX();
            auto y = hit -> GetY();
            auto z = hit -> GetZ();
            auto section = FindSection(x,y);
            auto charge = hit -> GetCharge();
            fHistPadPlane -> Fill(x,y,charge);
            fHistSectView[section] -> Fill(x,y,charge);
            fHistHeadView -> Fill(z,x,charge);
            fHistSideView -> Fill(z,y,charge);
        }
    }
    else if (optionString.Index("preview")==0)
    {
        if (fRawDataArray!=nullptr)
        {
            title = "Preview Data";
            TIter nextPad(fChannelArray);
            while (pad = (LKPhysicalPad*) nextPad())
            {
                auto padID = pad -> GetPadID();
                auto idx = pad -> GetDataIndex();
                if (idx<0)
                    continue;

                auto x = pad -> GetX();
                auto y = pad -> GetY();
                auto section = pad -> GetSection();

                auto channel = (GETChannel*) fRawDataArray -> At(idx);
                auto buffer = channel -> GetWaveformY();
                auto pedestal = channel -> GetPedestal();
                auto energy = channel -> GetEnergy();
                auto time = channel -> GetTime();

                fHistPadPlane -> Fill(x,y,energy);
                fHistSectView[section] -> Fill(x,y,energy);
                if (fDefinePositionByPixelIndex) {
                    for (auto tb=0; tb<512; ++tb) {
                        fHistHeadView -> Fill(tb,x,buffer[tb]);
                        fHistSideView -> Fill(tb,y,buffer[tb]);
                    }
                }
                else {
                    fHistHeadView -> Fill(time,x,energy);
                    fHistSideView -> Fill(time,y,energy);
                }
            }
        }
        else
            lk_error << "Raw-data array is null" << endl;
    }
    else if (optionString.Index("raw")==0)
    {
        if (fRawDataArray!=nullptr)
        {
            title = "Preview Data";
            TIter nextPad(fChannelArray);
            while (pad = (LKPhysicalPad*) nextPad())
            {
                auto padID = pad -> GetPadID();
                auto idx = pad -> GetDataIndex();
                if (idx<0)
                    continue;

                auto x = pad -> GetX();
                auto y = pad -> GetY();
                auto section = pad -> GetSection();

                auto channel = (GETChannel*) fRawDataArray -> At(idx);
                auto buffer = channel -> GetWaveformY();
                double energy = 0.;
                for (auto tb=0; tb<512; ++tb) energy += buffer[tb];
                fHistPadPlane -> Fill(x,y,energy);
                fHistSectView[section] -> Fill(x,y,energy);
                if (fDefinePositionByPixelIndex) {
                    for (auto tb=0; tb<512; ++tb) {
                        fHistHeadView -> Fill(tb,x,buffer[tb]);
                        fHistSideView -> Fill(tb,y,buffer[tb]);
                    }
                }
            }
        }
        else
            lk_error << "Raw-data array is null" << endl;
    }

    TString eventTitle = "";
    if (fRun!=nullptr) {
        auto inputFile = fRun -> GetInputFile();
        if (inputFile!=nullptr)
        {
            if (fAccumulateEvents==0)
                eventTitle = Form("%s (event %lld) [%s]", inputFile->GetName(), currentEventID, optionString.Data());
            else
                eventTitle = Form("%s (event %lld - %lld) [%s]", inputFile->GetName(), fAccumulateEvent1, fAccumulateEvent2, optionString.Data());
        }
        else
            eventTitle = Form("%s (event %lld) [%s]", fRun->GetRunName(), currentEventID, optionString.Data());
    }

    fGridPadPlane -> SetTitle(eventTitle);
}

LKChannelAnalyzer* LTPadPlane::GetChannelAnalyzer(int)
{
    if (fChannelAnalyzer==nullptr)
    {
        if (fPar->CheckPar("lamps/pulseFile")==false)
            fPar -> AddLine("lamps/pulseFile {lilak_common}/pulseReference.root");
        TString pulseFileName = fPar -> GetParString("lamps/pulseFile");
        fChannelAnalyzer = new LKChannelAnalyzer();
        fChannelAnalyzer -> SetPulse(pulseFileName);
        fChannelAnalyzer -> Print();
    }
    return fChannelAnalyzer;
}

int LTPadPlane::FindSection(double i, double j)
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

void LTPadPlane::SetNeighborPads(LKPad *pad0, LKPad *pad1)
{
    pad0 -> AddNeighborPad(pad1);
    pad1 -> AddNeighborPad(pad0);
}

void LTPadPlane::GetSectionParameters(int section, double &t1, double &t2, double &xmin, double &xmax, double &ymin, double &ymax)
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

void LTPadPlane::ClickedEventDisplay1(double xOnClick, double yOnClick)
{
    if (fSelectedSection==8) fHistEventDisplay1 = fHistHeadView;
    else if (fSelectedSection==9) fHistEventDisplay1 = fHistSideView;
    else fHistEventDisplay1 = fHistSectView[fSelectedSection];

    if (fHistEventDisplay1==nullptr)
        return;

    if (fSelectedSection>=8)
        return;

    int padID = FindPadID(xOnClick, yOnClick);

    if (padID<0) {
        lk_error << "Pad index is " << padID << ". (x, y) = (" << xOnClick << ", " << yOnClick << endl;
        return;
    }

    fSelPadID = padID;

    auto pad = (LKPhysicalPad*) fChannelArray -> At(fSelPadID);
    if (pad==nullptr) {
        lk_error << "pad at " << fSelPadID << " is nullptr" << endl;
        return;
    }

    pad -> Print();

    UpdateChannelBuffer();
}

void LTPadPlane::ClickedEventDisplay2(double xOnClick, double yOnClick)
{
    int selectedBin = fGridPadPlane -> FindBin(xOnClick, yOnClick);
    if (selectedBin==fBinNumberSideView) {
        fSelectedSection = 9;
        fHistEventDisplay1 = fHistSideView;
        fPadAxis1[0] = LKVector3::kZ;
        fPadAxis2[0] = LKVector3::kY;
    }
    else if (selectedBin==fBinNumberHeadView) {
        fSelectedSection = 8;
        fHistEventDisplay1 = fHistHeadView;
        fPadAxis1[0] = LKVector3::kZ;
        fPadAxis2[0] = LKVector3::kX;
    }
    else {
        fSelectedSection = FindSection(xOnClick,yOnClick);
        fHistEventDisplay1 = fHistSectView[fSelectedSection];
        fPadAxis1[0] = LKVector3::kX;
        fPadAxis2[0] = LKVector3::kY;
    }
    UpdateEventDisplay1();
}

void LTPadPlane::UpdateEventDisplay1()
{
    if (fHistEventDisplay1==nullptr) 
        return;

    fPadEventDisplay1 -> cd();
    if      (fEnergyMax==0) { fHistEventDisplay1 -> SetMinimum(fEnergyMin); fHistEventDisplay1 -> SetMaximum(-1111); }
    else if (fEnergyMax==1) { fHistEventDisplay1 -> SetMinimum(fEnergyMin); fHistEventDisplay1 -> SetMaximum(4200); }
    else
    {
        fHistEventDisplay1 -> SetMinimum(fEnergyMin);
        fHistEventDisplay1 -> SetMaximum(fEnergyMax);
        if (fEnergyMax>100)
            gStyle -> SetNumberContours(100);
        else
            gStyle -> SetNumberContours(fEnergyMax-fEnergyMin);
    }
    //fGridPadPlane -> Draw(fEventDisplayDrawOption);
    //fHistPadPlane -> Draw(fEventDisplayDrawOption);
    if (fSelectedSection<8)
    {
        if (fHistSectView[fSelectedSection]->GetEntries()==0)
            fGridSectView[fSelectedSection] -> Draw();
        else {
            fGridSectView[fSelectedSection] -> Draw();
            fHistSectView[fSelectedSection] -> Draw(fEventDisplayDrawOption+"same");
        }
        fGSelEventDisplay2 -> Set(0);
        fGSelEventDisplay2 -> SetPoint(0,fPosSectionCorner[fSelectedSection][0].X(), fPosSectionCorner[fSelectedSection][0].Y());
        fGSelEventDisplay2 -> SetPoint(1,fPosSectionCorner[fSelectedSection][1].X(), fPosSectionCorner[fSelectedSection][1].Y());
        fGSelEventDisplay2 -> SetPoint(2,fPosSectionCorner[fSelectedSection][2].X(), fPosSectionCorner[fSelectedSection][2].Y());
        fGSelEventDisplay2 -> SetPoint(3,fPosSectionCorner[fSelectedSection][3].X(), fPosSectionCorner[fSelectedSection][3].Y());
        fGSelEventDisplay2 -> SetPoint(4,fPosSectionCorner[fSelectedSection][0].X(), fPosSectionCorner[fSelectedSection][0].Y());
    }
    else {
        fHistEventDisplay1 -> Draw(fEventDisplayDrawOption);
        if (fSelectedSection==9) {
            fGSelEventDisplay2 -> Set(0);
            fGSelEventDisplay2 -> SetPoint(0, fXMin+50,  fYMin+50);
            fGSelEventDisplay2 -> SetPoint(1, fXMin+50,  fYMin+250);
            fGSelEventDisplay2 -> SetPoint(2, fXMin+250, fYMin+50);
            fGSelEventDisplay2 -> SetPoint(3, fXMin+50,  fYMin+50);
        }
        else if (fSelectedSection==8) {
            fGSelEventDisplay2 -> Set(0);
            fGSelEventDisplay2 -> SetPoint(0, fXMax-50,  fYMax-50);
            fGSelEventDisplay2 -> SetPoint(1, fXMax-50,  fYMax-250);
            fGSelEventDisplay2 -> SetPoint(2, fXMax-250, fYMax-50);
            fGSelEventDisplay2 -> SetPoint(3, fXMax-50,  fYMax-50);
        }
    }
    fPadEventDisplay2 -> cd();
    fGSelEventDisplay2 -> Draw("samel");
    fEventDisplayDrawOption = "colz";
}

void LTPadPlane::UpdateEventDisplay2()
{
    if (fHistEventDisplay2==nullptr) 
        return;
    fPadEventDisplay2 -> cd();
    if      (fEnergyMax==0) { fHistEventDisplay2 -> SetMinimum(fEnergyMin); fHistEventDisplay2 -> SetMaximum(-1111); }
    else if (fEnergyMax==1) { fHistEventDisplay2 -> SetMinimum(fEnergyMin);     fHistEventDisplay2 -> SetMaximum(4200); }
    else
    {
        fHistEventDisplay2 -> SetMinimum(fEnergyMin);
        fHistEventDisplay2 -> SetMaximum(fEnergyMax);
        if (fEnergyMax>100)
            gStyle -> SetNumberContours(100);
        else
            gStyle -> SetNumberContours(fEnergyMax-fEnergyMin);
    }
    //fHistEventDisplay2 -> Draw(fEventDisplayDrawOption);
    fGridPadPlane -> Draw();
    fHistPadPlane -> Draw(fEventDisplayDrawOption+"same");
    fEventDisplayDrawOption = "colz";
}
