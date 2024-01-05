#ifndef LTHELIXTRACKFINDINGTASK_HH
#define LTHELIXTRACKFINDINGTASK_HH

#define DEBUG_STEP -1

#include "TClonesArray.h"
#include "TGraphErrors.h"

#include "LKTask.h"
#include "LKHelixTrack.h"
#include "LKHit.h"
#include "LKHitArray.h"

#include "LKDetector.h"
#include "LKDetectorPlane.h"

#include <vector>
using namespace std;

class LTHelixTrackFindingTask : public LKTask
{
    public:
        LTHelixTrackFindingTask() : LKTask("LTHelixTrackFindingTask","LTHelixTrackFindingTask") {}
        virtual ~LTHelixTrackFindingTask() {}

        virtual bool Init();
        virtual void Exec(Option_t*);

        void SetTrackPersistency(bool val) { fPersistency = val; }

        enum StepNo : int {
            kStepInitArray,
            kStepNewTrack,
            kStepRemoveTrack,
            kStepInitTrack,
            kStepInitTrackAddHit,
            kStepContinuum,
            kStepContinuumAddHit,
            kStepExtrapolation,
            kStepExtrapolationAddHit,
            kStepConfirmation,
            kStepFinalizeTrack,
            kStepNextPhase,
            kStepEndEvent,
            kStepEndOfEvent,
        };

        bool ExecStep();

        bool ExecStepUptoTrackNum(Int_t numTracks);

        LKHelixTrack *GetCurrentTrack() const { return fCurrentTrack; }

        void SetHitBranchName(TString name) { fBranchNameHit = name; }
        void SetTrackletBranchName(TString name) { fBranchNameTracklet = name; }

    private:
        int StepInitArray();
        int StepNewTrack();
        int StepRemoveTrack();
        int StepInitTrack();
        int StepInitTrackAddHit();
        int StepContinuum();
        int StepContinuumAddHit();
        int StepExtrapolation();
        int StepExtrapolationAddHit();
        int StepConfirmation();
        int StepFinalizeTrack();
        int StepNextPhase();
        int StepEndEvent();

        void ReturnBadHitsToPadPlane();

        double CorrelateHitWithTrackCandidate(LKHelixTrack *track, LKHit *hit);
        double CorrelateHitWithTrack(LKHelixTrack *track, LKHit *hit, Double_t scale=1);

        int CheckParentTrackID(LKHit *hit);
        bool CheckTrackQuality(LKHelixTrack *track);
        double CheckTrackContinuity(LKHelixTrack *track);
        bool CheckHitDistInAlphaIsLargerThanQuarterPi(LKHelixTrack *track, Double_t dLength);

        bool BuildAndConfirmTrack(LKHelixTrack *track, bool &tailToHead);
        bool AutoBuildByExtrapolation(LKHelixTrack *track, bool &buildHead, Double_t &extrapolationLength);
        bool AutoBuildAtPosition(LKHelixTrack *track, TVector3 p, bool &tailToHead, Double_t &extrapolationLength, Double_t scale=1);

    private:
        LKDetector *fTpc = nullptr;
        LKDetectorPlane *fPadPlane = nullptr;
        TClonesArray *fHitArray = nullptr;
        TClonesArray *fTrackArray = nullptr;

        TString fBranchNameHit = "Hit";
        TString fBranchNameTracklet = "Tracklet";

        bool fPersistency = true;

        LKHitArray *fTrackHits = nullptr;
        LKHitArray *fCandHits = nullptr;
        LKHitArray *fGoodHits = nullptr;
        LKHitArray *fBadHits  = nullptr;

        Double_t fDefaultScale = 2.5;
        Double_t fTrackWCutLL  = 5.;  ///< Track width cut low limit
        Double_t fTrackWCutHL  = 15.; ///< Track width cut high limit
        Double_t fTrackHCutLL  = 5.;  ///< Track height cut low limit
        Double_t fTrackHCutHL  = 15.; ///< Track height cut high limit
        LKVector3::Axis fReferenceAxis = LKVector3::kZ;

        Int_t fPhaseIndex = 0;

        Int_t fMinHitsToFitInitTrack = 7; ///< try track fit if track has more than this number of hits in track
        Int_t fCutMinNumHitsInitTrack = 10; ///
        Int_t fCutMaxNumHitsInitTrack = 15; ///< return hits if track has more than this number of hits within initialization stage
        Int_t fCutMinNumHitsFinalTrack = 15; ///< remove track if track has smaller than this number of hits within initialization stage
        Double_t fCutMinHelixRadius = 30.; ///< helix radius cut for initialization stage
        Double_t fTrackLengthCutScale = 2.5; ///< track length cut for initialization stage is [this_var] * track -> GetRMSR()
        Double_t fCutdkInExpectedTrackPath = 4.; // the correlation distance cut through helix axis during the track path between two hits

        LKHelixTrack *fCurrentTrack = nullptr;

        Int_t fNextStep = StepNo::kStepInitArray;
        Int_t fNumCandHits;
        Int_t fNumGoodHits;
        Int_t fNumBadHits;

        TCanvas *fCvsCurrentTrack = nullptr;
        TGraphErrors *fGraphCurrentTrackPoint = nullptr;

#ifdef DEBUG_STEP
        Int_t dCountStep = 0;
#endif

    ClassDef(LTHelixTrackFindingTask, 2)
};

#endif
