#ifndef LTPADPLANE_HH
#define LTPADPLANE_HH

#include "LKEvePlane.h"
#include "TH2Poly.h"
#include "TPad.h"
#include "TF1.h"

class LTPadPlane : public LKEvePlane
{
    public:
        LTPadPlane();
        virtual ~LTPadPlane() {};

        virtual void Clear(Option_t *option = "");
        virtual void Print(Option_t *option = "") const;
        virtual bool Init();

        virtual double PadDisplacement() const;
        virtual bool IsInBoundary(double i, double j);

        virtual LKChannelAnalyzer* GetChannelAnalyzer(int i=0);

    public:
        virtual TH2* GetHistEventDisplay1(Option_t *option="-1");
        virtual TH2* GetHistEventDisplay2(Option_t *option="-1");

        virtual int FindPadID(int section, int layer, int row);
        virtual int FindPadID(int cobo, int asad, int aget, int chan);
        virtual int FindPadID(double i, double j);

        virtual void FillDataToHist(Option_t *option=""); ///< Implementation recommanded. Fill data to histograms.
        virtual void ClickedEventDisplay1(double xOnClick, double yOnClick);
        virtual void ClickedEventDisplay2(double xOnClick, double yOnClick);
        virtual void UpdateEventDisplay1();
        virtual void UpdateEventDisplay2();

    private:
        bool fDefinePositionByPixelIndex = true;

        TF1* fFuncXRightBound = nullptr;
        TF1* fFuncXRightBoundInverse = nullptr;

        double fXMin = -550;
        double fXMax = +550;
        double fYMin = -550;
        double fYMax = +550;
        double fZMin = -350;
        double fZMax = 1350;
        int fZBins = 128;

        double fHistZMin = 0.1;

        TVector2 fPosSectionCorner[8][4];

        double fRMin = 105.;
        double fRMax = 545-41.5;
        double fPadGap = 0.0;
        double fPadWidth = 3.5;
        double fPadHeight = 10.5;
        double fRadiusLayer0 = fRMin + 0.5*fPadHeight;
        double fRTopCut = fRMax;
        int fNumLayers = 42;

        int fNumTbs = 512;
        int fMaxWaveformY = 4100;

        bool fZoomZoomPressed = false;
        int fBinZoomZoomButton;
        int fBinNumberHeadView;
        int fBinNumberSideView;

        double fTanPi1o8;
        double fTanPi3o8;
        double fTanPi5o8;
        double fTanPi7o8;
        double fCosPiNo4[8];
        double fSinPiNo4[8];

        int fStartRowFrom[42] = {0};
        int fNumHalfRowsInLayerInput[42] = {12,13,14,16,17,18,19,20,21,22,24,25,26,27,28,29,31,32,33,34,35,36,37,39,40,41,42,43,44,45,47,48,49,50,51,52,53,55,40,23,13,5};
        int fNumHalfRowsInLayer[50] = {0};
        int fNumPadsDownToLayer[50] = {0}; ///< in single section
        int fNumSkippedHalfRows[50] = {0};

        double fXSpacing;
        double fRSpacing;
        double fDCX = -1.5; ///< Displacement of pad cut line in x-axis
        //double fDCX = 0; ///< Displacement of pad cut line in x-axis

        const bool fDoCutTopBoundary = true;
        const bool fDoCutSideBoundary = true;

        int ****fMapCAACToPadID;
        int ****fMapSLRPToPadID;

        int fSelectedSection = 2;

        TH2Poly* fGridSectView[8];
        TH2Poly* fHistSectView[8];
        TH2Poly* fGridPadPlane = nullptr;
        TH2Poly* fHistPadPlane = nullptr;
        TH2D* fHistHeadView = nullptr;
        TH2D* fHistSideView = nullptr;

        TGraph* fGraphSectionBoundary1 = nullptr;
        TGraph* fGraphSectionBoundary2 = nullptr;
        TGraph* fGraphPadBoundary = nullptr;
        TGraph* fGraphPadBoundaryNb[10];

    private:
        int ConvHP(int bin) { return bin-1; }
        int ConvPH(int bin) { return bin+1; }

    public:
        TH2Poly* GetHistSection(int selectSection);

        int FindSection(double i, double j);
        void SetNeighborPads(LKPad *pad0, LKPad *pad1);
        void GetSectionParameters(int section, double &t1, double &t2, double &xmin, double &xmax, double &ymin, double &ymax);

        double GetPadCutBoundaryYAtX(double x) { return fTanPi3o8*(x-fDCX); }
        double GetPadCutBoundaryXAtR(double r) { return (fDCX*fTanPi3o8*fTanPi3o8+sqrt((r*r-fDCX*fDCX)*fTanPi3o8*fTanPi3o8+r*r)) / (fTanPi3o8*fTanPi3o8+1); }
        double GetPadCenterYBoundAtX(double x) { return fTanPi3o8*(x-.5) + 10; }

        double GetRMin()            const { return fRMin; }
        double GetRMax()            const { return fRMax; }
        double GetPadGap()          const { return fPadGap; }
        double GetPadWidth()        const { return fPadWidth; }
        double GetPadHeight()       const { return fPadHeight; }
        int GetNumLayer()        const { return fNumLayers; }
        double GetPadXSpacing()     const { return fXSpacing; }
        double GetPadRSpacing()     const { return fRSpacing; }

    ClassDef(LTPadPlane, 2)
};

#endif
