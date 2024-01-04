#ifndef LAMPSTPCPADPLANE_HH
#define LAMPSTPCPADPLANE_HH

#include "LKDetectorPlane.h"
#include "TF1.h"
#include "TH2Poly.h"

class LTPadPlane : public LKDetectorPlane
{
    public:
        LTPadPlane();
        virtual ~LTPadPlane() {}

        static LTPadPlane* fInstance;
        static LTPadPlane* GetPlane();

        virtual bool Init();

        virtual Int_t FindPadID(Double_t i, Double_t j);
        virtual Int_t FindPadID(Int_t section, Int_t layer, Int_t row);
        virtual Double_t PadDisplacement() const;
        virtual bool IsInBoundary(Double_t i, Double_t j);
        virtual void DrawFrame(Option_t *option = "");
        virtual TH2* GetHist(Option_t *option = "-1");
        virtual TCanvas *GetCanvas(Option_t *optiont = "");
        virtual bool SetDataFromBranch();
        virtual void FillDataToHist();

        TH2* GetHist(Int_t selectSection);
        void Draw(Option_t *option="");

        Int_t FindPadID(Int_t cobo, Int_t aget, Int_t asad, Int_t chan);
        LKPad *GetPadFromEleID(Int_t cobo, Int_t aget, Int_t asad, Int_t chan);
        int GetElectronicsID(int cobo, int aget, int asad, int channelID) { return 100000*cobo + 10000*aget + 100*asad + channelID; }

        void GetSectionParameters(Int_t section, Double_t &t1, Double_t &t2, Double_t &xmin, Double_t &xmax, Double_t &ymin, Double_t &ymax);
        void CreateHistograms();

        static void MouseClickEvent1();
        static void MouseClickEvent2();

        void ZoomInWindow(Int_t bin, Double_t x, Double_t y);
        void SelectAndDrawChannel(Int_t bin, Double_t x, Double_t y);

        //virtual void DriftElectron(TVector3 posGlobal, TVector3 &posFinal, double &driftLength);
        //virtual void DriftElectronBack(LKPad* pad, double tb, TVector3 &posReco, double &driftLength);

        void WriteCurrentChannel(TString name="");

    private:
        LKPad *NewPad(Int_t s, Int_t r, Int_t l);
        void SetNeighborPads(LKPad *pad0, LKPad *pad1);
        Int_t FindSection(Double_t i, Double_t j);

    private:
        TF1* fFuncXRightBound = nullptr;
        TF1* fFuncXRightBoundInverse = nullptr;

        Double_t fXMin = -550;
        Double_t fXMax = +550;
        Double_t fYMin = -550;
        Double_t fYMax = +550;
        Double_t fZMin = -350;
        Double_t fZMax = 1350;
        Int_t fZBins = 512;

        TVector2 fPosSectionCorner[8][4];

        Double_t fRMin = 105.;
        Double_t fRMax = 545-41.5;
        //Double_t fRMin = 100.;
        //Double_t fRMax = 510;
        Double_t fPadGap = 0.5;
        Double_t fPadWidth = 3.;
        Double_t fPadHeight = 10.;
        Double_t fRadiusLayer0 = fRMin + 0.5*fPadHeight;
        Double_t fRTopCut = fRMax;
        Int_t fNumLayers = 42;

        Int_t fNumTbs = 512;
        Int_t fMaxWaveformY = 4100;

        Int_t fBinNumberTopView;
        Int_t fBinNumberSideView;

        Double_t fTanPi1o8;
        Double_t fTanPi3o8;
        Double_t fTanPi5o8;
        Double_t fTanPi7o8;
        Double_t fCosPiNo4[8];
        Double_t fSinPiNo4[8];

        int fStartRowFrom[42] = {0};
        int fNumHalfRowsInLayerInput[42] = {12,13,14,16,17,18,19,20,21,22,24,25,26,27,28,29,31,32,33,34,35,36,37,39,40,41,42,43,44,45,47,48,49,50,51,52,53,55,40,23,13,5};
        int fNumHalfRowsInLayer[50] = {0};
        int fNumPadsDownToLayer[50] = {0}; ///< in single section
        int fNumSkippedHalfRows[50] = {0};

        Double_t fXSpacing;
        Double_t fRSpacing;
        Double_t fDCX = -1.5; ///< Displacement of pad cut line in x-axis
        //Double_t fDCX = 0; ///< Displacement of pad cut line in x-axis

        const Bool_t fDoCutTopBoundary = true;
        const Bool_t fDoCutSideBoundary = true;

        int ****fMapCAACToPadID;

        int fSelectedSection = 0;

        TH2Poly* fHistPadPlaneSection[8];
        TH2Poly* fFramePadPlane = nullptr;
        TH2Poly* fHistPadPlane = nullptr;
        //TGraph* fFrameSideView = nullptr;
        TH2D* fHistWaveform = nullptr;
        TH2D* fHistTopView = nullptr;
        TH2D* fHistSideView = nullptr;
        TH1D* fHistChannel1 = nullptr;
        TH1D* fHistChannel2 = nullptr;
        TGraph* fGraphSectionBoundary1 = nullptr;
        TGraph* fGraphSectionBoundary2 = nullptr;
        TGraph* fGraphPadBoundary = nullptr;

        bool fUseChannel1 = true;

        TClonesArray* fBufferArray = nullptr;
        TClonesArray* fHitArray = nullptr;

        TString fFillType = "Buffer";

        LKPad *fSelectedPad = nullptr;

        int ConvHP(int bin) { return bin-1; }
        int ConvPH(int bin) { return bin+1; }

    public:
        double GetPadCutBoundaryYAtX(double x) { return fTanPi3o8*(x-fDCX); }
        double GetPadCutBoundaryXAtR(double r) { return (fDCX*fTanPi3o8*fTanPi3o8+sqrt((r*r-fDCX*fDCX)*fTanPi3o8*fTanPi3o8+r*r)) / (fTanPi3o8*fTanPi3o8+1); }
        double GetPadCenterYBoundAtX(double x) { return fTanPi3o8*(x-.5) + 10; }

        Double_t GetRMin()            const { return fRMin; }
        Double_t GetRMax()            const { return fRMax; }
        Double_t GetPadGap()          const { return fPadGap; }
        Double_t GetPadWidth()        const { return fPadWidth; }
        Double_t GetPadHeight()       const { return fPadHeight; }
        Int_t GetNumLayer()        const { return fNumLayers; }
        Double_t GetPadXSpacing()     const { return fXSpacing; }
        Double_t GetPadRSpacing()     const { return fRSpacing; }

        ClassDef(LTPadPlane, 1)
};

#endif
