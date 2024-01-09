#ifndef LAMPSTPC_HH
#define LAMPSTPC_HH

#include "LKDetector.h"
#include "LKLogger.h"
#include "LKChannelAnalyzer.h"
#include "GETChannel.h"

/*
 * Remove this comment block after reading it through
 * or use print_example_comments=False option to omit printing
 *
 * # Example LILAK detector class
 *
 */

class LAMPSTPC : public LKDetector
{
    public:
        LAMPSTPC();
        virtual ~LAMPSTPC() { ; }

        LAMPSTPC(TString parName); ///< Init with parName
        LAMPSTPC(LKParameterContainer* par); ///< Init with par

        void Print(Option_t *option="") const;
        bool Init();
        bool InitChannelAnalyzer();
        bool BuildGeometry();
        bool BuildDetectorPlane();
        bool IsInBoundary(Double_t x, Double_t y, Double_t z);
        bool GetEffectiveDimension(Double_t &x1, Double_t &y1, Double_t &z1, Double_t &x2, Double_t &y2, Double_t &z2);

    private:
        LKChannelAnalyzer* fChannelAnalyzer[18];

    ClassDef(LAMPSTPC,1);
};

#endif
