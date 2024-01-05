#include "LKRun.h"
#include "LTHelixTrackFindingTask.h"

#include <iostream>

ClassImp(LTHelixTrackFindingTask)

bool LTHelixTrackFindingTask::Init()
{
    fPadPlane = fRun -> GetDetectorPlane();

    fHitArray = (TClonesArray *) fRun -> GetBranchA(fBranchNameHit);
    fTrackArray = fRun -> RegisterBranchA("Track","LKHelixTrack",100);

    fTrackHits = new LKHitArray();
    fCandHits = new LKHitArray();
    fGoodHits = new LKHitArray();
    fBadHits  = new LKHitArray();

    //fTrackHits = fRun -> RegisterBranchA("TrackHit");
    //fCandHits  = fRun -> RegisterBranchA("CandHit");
    //fGoodHits  = fRun -> RegisterBranchA("GoodHit");
    //fBadHits   = fRun -> RegisterBranchA("BadHit");

    if (fPar -> CheckPar("LTHelixTrackFindingTask/defaultScale" )) fDefaultScale  = fPar -> GetParDouble("LTHelixTrackFindingTask/defaultScale");
    if (fPar -> CheckPar("LTHelixTrackFindingTask/trackWCutLL"  )) fTrackWCutLL   = fPar -> GetParDouble("LTHelixTrackFindingTask/trackWCutLL");
    if (fPar -> CheckPar("LTHelixTrackFindingTask/trackWCutHL"  )) fTrackWCutHL   = fPar -> GetParDouble("LTHelixTrackFindingTask/trackWCutHL");
    if (fPar -> CheckPar("LTHelixTrackFindingTask/trackHCutLL"  )) fTrackHCutLL   = fPar -> GetParDouble("LTHelixTrackFindingTask/trackHCutLL");
    if (fPar -> CheckPar("LTHelixTrackFindingTask/trackHCutHL"  )) fTrackHCutHL   = fPar -> GetParDouble("LTHelixTrackFindingTask/trackHCutHL");
    if (fPar -> CheckPar("LTHelixTrackFindingTask/tpcBFieldAxis")) fReferenceAxis = fPar -> GetParAxis  ("LTHelixTrackFindingTask/tpcBFieldAxis");

    fNextStep = StepNo::kStepInitArray;

    return true;
}

void LTHelixTrackFindingTask::Exec(Option_t*)
{
#ifdef DEBUG_STEP
    dCountStep = 0;
#endif
    fNextStep = StepNo::kStepInitArray;
    while (ExecStep()) {}
}

bool LTHelixTrackFindingTask::ExecStep()
{
         if (fNextStep==kStepEndOfEvent)          return false; 
    else if (fNextStep==kStepInitArray)           fNextStep = StepInitArray();
    else if (fNextStep==kStepNewTrack)            fNextStep = StepNewTrack();
    else if (fNextStep==kStepRemoveTrack)         fNextStep = StepRemoveTrack();
    else if (fNextStep==kStepInitTrack)           fNextStep = StepInitTrack();
    else if (fNextStep==kStepInitTrackAddHit)     fNextStep = StepInitTrackAddHit();
    else if (fNextStep==kStepContinuum)           fNextStep = StepContinuum();
    else if (fNextStep==kStepContinuumAddHit)     fNextStep = StepContinuumAddHit();
    else if (fNextStep==kStepExtrapolation)       fNextStep = StepExtrapolation();
    else if (fNextStep==kStepExtrapolationAddHit) fNextStep = StepExtrapolationAddHit();
    else if (fNextStep==kStepConfirmation)        fNextStep = StepConfirmation();
    else if (fNextStep==kStepFinalizeTrack)       fNextStep = StepFinalizeTrack();
    else if (fNextStep==kStepNextPhase)           fNextStep = StepNextPhase();
    else if (fNextStep==kStepEndEvent)            fNextStep = StepEndEvent();

#ifdef DEBUG_STEP
    if (DEBUG_STEP>0 && dCountStep==DEBUG_STEP) {
        lk_strong << "The step is forced to break by setting step cut " << DEBUG_STEP << endl;
        return false;
    }
         if (fNextStep==kStepInitArray)           lk_debug << "[" << dCountStep++ << "] NextStep is StepInitArray" << endl;
    else if (fNextStep==kStepNewTrack)            lk_debug << "[" << dCountStep++ << "] NextStep is StepNewTrack" << endl;
    else if (fNextStep==kStepRemoveTrack)         lk_debug << "[" << dCountStep++ << "] NextStep is StepRemoveTrack" << endl;
    else if (fNextStep==kStepInitTrack)           lk_debug << "[" << dCountStep++ << "] NextStep is StepInitTrack" << endl;
    else if (fNextStep==kStepInitTrackAddHit)     lk_debug << "[" << dCountStep++ << "] NextStep is StepInitTrackAddHit" << endl;
    else if (fNextStep==kStepContinuum)           lk_debug << "[" << dCountStep++ << "] NextStep is StepContinuum" << endl;
    else if (fNextStep==kStepContinuumAddHit)     lk_debug << "[" << dCountStep++ << "] NextStep is StepContinuumAddHit" << endl;
    else if (fNextStep==kStepExtrapolation)       lk_debug << "[" << dCountStep++ << "] NextStep is StepExtrapolation" << endl;
    else if (fNextStep==kStepExtrapolationAddHit) lk_debug << "[" << dCountStep++ << "] NextStep is StepExtrapolationAddHit" << endl;
    else if (fNextStep==kStepConfirmation)        lk_debug << "[" << dCountStep++ << "] NextStep is StepConfirmation" << endl;
    else if (fNextStep==kStepFinalizeTrack)       lk_debug << "[" << dCountStep++ << "] NextStep is StepFinalizeTrack" << endl;
    else if (fNextStep==kStepNextPhase)           lk_debug << "[" << dCountStep++ << "] NextStep is StepNextPhase" << endl;
    else if (fNextStep==kStepEndEvent)            lk_debug << "[" << dCountStep++ << "] NextStep is StepEndEvent" << endl;

    if (fCurrentTrack != nullptr)
        fCurrentTrack -> Print();
#endif
    return true;
}

bool LTHelixTrackFindingTask::ExecStepUptoTrackNum(Int_t numTracks)
{
    if (fNextStep==kStepEndOfEvent)
        return false;

    while (ExecStep()) {
        if (fNextStep==kStepNewTrack && fTrackArray -> GetEntriesFast() >= numTracks)
            break;
    }

    return true;
}

int LTHelixTrackFindingTask::StepInitArray()
{
#ifdef DEBUG_STEP
    lk_debug << "StepInitArray" << endl;
#endif
    fCurrentTrack = nullptr;

    fPhaseIndex = 0;

    fPadPlane -> ResetHitMap();
    fPadPlane -> SetHitArray(fHitArray);

    fTrackArray -> Clear("C");
    fTrackHits -> Clear();
    fCandHits -> Clear();
    fGoodHits -> Clear();
    fBadHits -> Clear(); 

    return kStepNewTrack;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LTHelixTrackFindingTask::StepNewTrack()
{
#ifdef DEBUG_STEP
    lk_debug << "StepNewTrack" << endl;
#endif

    fTrackHits -> Clear();
    fCandHits -> Clear();
    fGoodHits -> Clear();
    ReturnBadHitsToPadPlane();

    LKHit *hit = fPadPlane -> PullOutNextFreeHit();
    if (hit == nullptr) {
#ifdef DEBUG_STEP
        lk_debug << "no more hits" << endl;
#endif
        return kStepNextPhase;
    }

    Int_t idx = fTrackArray -> GetEntries();
    fCurrentTrack = new ((*fTrackArray)[idx]) LKHelixTrack(idx);
    fCurrentTrack -> SetReferenceAxis(fReferenceAxis);
    fCurrentTrack -> AddHit(hit);
    fGoodHits -> AddHit(hit);

#ifdef DEBUG_STEP
    fCurrentTrack -> Print(">");
#endif

    return kStepInitTrack;
}

int LTHelixTrackFindingTask::StepRemoveTrack()
{
#ifdef DEBUG_STEP
    lk_debug << "StepRemoveTrack" << endl;
#endif
    fGoodHits -> Clear();
    fTrackHits -> Clear();
    fCandHits -> Clear();

#ifdef DEBUG_STEP
    lk_debug << "Return bad hits" << endl;
#endif
    ReturnBadHitsToPadPlane();

    auto trackHits = fCurrentTrack -> GetHitArray();
    Int_t numTrackHits = trackHits -> GetNumHits();
#ifdef DEBUG_STEP
    lk_debug << "Return track hits" << endl;
    lk_debug << numTrackHits << endl;
#endif
    for (Int_t iTrackHit = 0; iTrackHit < numTrackHits; ++iTrackHit) {
        auto trackHit = (LKHit *) trackHits -> GetHit(iTrackHit);
#ifdef DEBUG_STEP
        lk_debug << iTrackHit << endl;
        trackHit -> Print();
#endif
        trackHit -> AddTrackCand(-1);
        fPadPlane -> AddHit(trackHit);
    }
    fTrackArray -> Remove(fCurrentTrack);
#ifdef DEBUG_STEP
    lk_debug << "return" << endl;
#endif

    fCurrentTrack = nullptr;

    return kStepNewTrack;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LTHelixTrackFindingTask::StepInitTrack()
{
#ifdef DEBUG_STEP
    lk_debug << "StepInitTrack" << endl;
#endif
    fCandHits -> Clear();
    fPadPlane -> PullOutNeighborHits(fGoodHits, fCandHits);
    fGoodHits -> MoveHitsTo(fTrackHits);
#ifdef DEBUG_STEP
    lk_debug << "After pull hits" << endl;
    lk_debug << "good hits>>" << endl;
    fGoodHits -> Print();
    lk_debug << "cand hits>>" << endl;
    fCandHits -> Print();
#endif
    fNumCandHits = fCandHits -> GetEntriesFast();
    if (fNumCandHits == 0)
        return kStepRemoveTrack;

    fCandHits -> SortByDistanceTo(fCurrentTrack -> GetMean(), true);
#ifdef DEBUG_STEP
    lk_debug << "After sort" << endl;
#endif
    return kStepInitTrackAddHit;
}

int LTHelixTrackFindingTask::StepInitTrackAddHit()
{
#ifdef DEBUG_STEP
    lk_debug << "StepInitTrackAddHit, fNumCandHits=" << fNumCandHits << endl;
#endif
    auto candHit = (LKHit *) fCandHits -> GetLastHit();
    fCandHits -> RemoveLastHit();

    Double_t quality;
    if (fCurrentTrack -> IsHelix())
        quality = CorrelateHitWithTrack(fCurrentTrack, candHit);
    else
        quality = CorrelateHitWithTrackCandidate(fCurrentTrack, candHit);

    if (quality > 0) {
        fGoodHits -> AddHit(candHit);
        fCurrentTrack -> AddHit(candHit);
        fCurrentTrack -> FitPlane(); // XXX should comment out

        auto numHitsInTrack = fCurrentTrack -> GetNumHits();
#ifdef DEBUG_STEP
        lk_debug << "Quality is " << quality << ", hit is added, numHitsInTrack=" << numHitsInTrack << endl;
#endif

        if (numHitsInTrack > fCutMaxNumHitsInitTrack) {
#ifdef DEBUG_STEP
            lk_debug << "track has too many hits. remove track" << endl;
#endif
            Int_t numCandHits2 = fCandHits -> GetEntriesFast();
            for (Int_t iCand = 0; iCand < numCandHits2; ++iCand)
                fPadPlane -> AddHit((LKHit *) fCandHits -> GetHit(iCand));
            fCandHits -> Clear("C");
            return kStepRemoveTrack;
        }

        if (numHitsInTrack >= fMinHitsToFitInitTrack) {
#ifdef DEBUG_STEP
            lk_debug << "fit track since track has enough number of hits" << endl;
#endif
            //fCurrentTrack -> SetIsHelix();
            fCurrentTrack -> Fit();
#ifdef DEBUG_STEP
            lk_debug << "After fit, " << endl;
            fCurrentTrack -> Print();
            lk_debug << "After fit, track radius : " << fCurrentTrack -> GetHelixRadius() << " >? " << fCutMinHelixRadius << endl;
            lk_debug << "RMS  " << fCurrentTrack -> GetRMS () << endl;
            lk_debug << "RMSR " << fCurrentTrack -> GetRMSR() << endl;
            lk_debug << "RMST " << fCurrentTrack -> GetRMST() << endl;
            lk_debug << "After fit, track length : " << fCurrentTrack -> TrackLength() << " >? " << fTrackLengthCutScale << " * " << fCurrentTrack -> GetRMST() << endl;
            lk_debug << "... " << fTrackLengthCutScale * fCurrentTrack -> GetRMST() << endl;
#endif
            if (numHitsInTrack > fCutMinNumHitsInitTrack &&
                    fCurrentTrack -> GetHelixRadius() > fCutMinHelixRadius &&
                    fCurrentTrack -> TrackLength() > fTrackLengthCutScale * fCurrentTrack -> GetRMST()) {
#ifdef DEBUG_STEP
                lk_debug << "Track is good" << endl;
#endif
                return kStepContinuum;
            }
            else {
#ifdef DEBUG_STEP
                lk_debug << "Track need more hits" << endl;
#endif
                fCurrentTrack -> FitPlane();
            }
        }

#ifdef DEBUG_STEP
        lk_debug << "Quality is " << quality << ", track fit, numHitsInTrack=" << numHitsInTrack  << endl;
#endif
        //fCurrentTrack -> FitPlane(); // XXX should comment in
    }
    else {
#ifdef DEBUG_STEP
        lk_debug << "Quality is " << quality << ", bad hit" << endl;
#endif
        fBadHits -> AddHit(candHit);
    }

    fNumCandHits = fCandHits -> GetEntriesFast();
    if (fNumCandHits==0)
        return kStepInitTrack;

    return kStepInitTrackAddHit;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LTHelixTrackFindingTask::StepContinuum()
{ 
#ifdef DEBUG_STEP
    lk_debug << "StepContinuum" << endl;
    lk_debug << "goodhits" << endl;
    fGoodHits -> PrintHits(1);
#endif
    fPadPlane -> PullOutNeighborHits(fGoodHits, fCandHits);
    fGoodHits -> MoveHitsTo(fTrackHits);

    fNumCandHits = fCandHits -> GetEntries();
    if (fNumCandHits==0)
        return kStepExtrapolation;

    fCandHits -> SortByCharge(false);
    return kStepContinuumAddHit;
}

int LTHelixTrackFindingTask::StepContinuumAddHit()
{
#ifdef DEBUG_STEP
    lk_debug << "StepContinuumAddHit" << endl;
#endif
    for (Int_t iHit = 0; iHit < fNumCandHits; iHit++) {
        LKHit *candHit = (LKHit *) fCandHits -> GetLastHit();
        fCandHits -> RemoveLastHit();

        Double_t quality = 0; 
        if (CheckParentTrackID(candHit) == -2)
            quality = CorrelateHitWithTrack(fCurrentTrack, candHit);

        if (quality > 0) {
            fGoodHits -> AddHit(candHit);
            fCurrentTrack -> AddHit(candHit);
            fCurrentTrack -> Fit();
#ifdef DEBUG_STEP
            lk_debug << iHit << "(/" << fNumCandHits << ") is good, quality is " << quality << endl;
#endif
        } else {
            fBadHits -> AddHit(candHit);
#ifdef DEBUG_STEP
            lk_debug << iHit << "(/" << fNumCandHits << ") is bad" << endl;
#endif
        }
    }

    return kStepContinuum;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LTHelixTrackFindingTask::StepExtrapolation()
{ 
#ifdef DEBUG_STEP
    lk_debug << "StepExtrapolation" << endl;
#endif
    ReturnBadHitsToPadPlane();

    return kStepExtrapolationAddHit;
}

int LTHelixTrackFindingTask::StepExtrapolationAddHit()
{
#ifdef DEBUG_STEP
    lk_debug << "StepExtrapolationAddHit" << endl;
#endif
    Int_t countIteration = 0;
    bool buildFromHeadOrTail = true;
    Double_t extrapolationLength = 0;
    while (AutoBuildByExtrapolation(fCurrentTrack, buildFromHeadOrTail, extrapolationLength)) {
        if (++countIteration > 200)
            break;
    }

    countIteration = 0;
    buildFromHeadOrTail = !buildFromHeadOrTail;
    extrapolationLength = 0;
    while (AutoBuildByExtrapolation(fCurrentTrack, buildFromHeadOrTail, extrapolationLength)) {
        if (++countIteration > 200)
            break;
    }

    ReturnBadHitsToPadPlane();

    bool isGood = CheckTrackQuality(fCurrentTrack);

    if (isGood) {
        return kStepConfirmation;
#ifdef DEBUG_STEP
        lk_debug << "Track is good! " << endl;
#endif
    }

#ifdef DEBUG_STEP
    lk_debug << "Track is bad! " << endl;
#endif

    return kStepRemoveTrack;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LTHelixTrackFindingTask::StepConfirmation()
{
#ifdef DEBUG_STEP
    lk_debug << "StepConfirmation" << endl;
#endif

    bool tailToHead = false;
    LKVector3 pTail(fCurrentTrack -> PositionAtTail(), fReferenceAxis);
    LKVector3 pHead(fCurrentTrack -> PositionAtHead(), fReferenceAxis);

    if (pHead.K() > pHead.K())
        tailToHead = true;

    ReturnBadHitsToPadPlane();
#ifdef DEBUG_STEP
    lk_debug << "First Try" << endl;
    fCurrentTrack -> Print(">");
#endif
    if (BuildAndConfirmTrack(fCurrentTrack, tailToHead) == false) {
#ifdef DEBUG_STEP
        lk_debug << "failed in confirm!" << endl;
        fCurrentTrack -> Print(">");
#endif
        return kStepRemoveTrack;
    }

    tailToHead = !tailToHead; 

    ReturnBadHitsToPadPlane();
#ifdef DEBUG_STEP
    lk_debug << "Second Try" << endl;
#endif
    if (BuildAndConfirmTrack(fCurrentTrack, tailToHead) == false)
        return kStepRemoveTrack;

#ifdef DEBUG_STEP
    lk_debug << "End" << endl;
#endif

    ReturnBadHitsToPadPlane();

    return kStepFinalizeTrack;
}


int LTHelixTrackFindingTask::StepFinalizeTrack()
{
    auto trackHits = fCurrentTrack -> GetHitArray();
    Int_t trackID = fCurrentTrack -> GetTrackID();
    Int_t numTrackHits = trackHits -> GetNumHits();
    for (Int_t iTrackHit = 0; iTrackHit < numTrackHits; ++iTrackHit) {
        auto trackHit = (LKHit *) trackHits -> GetHit(iTrackHit);
        trackHit -> AddTrackCand(trackID);
        fPadPlane -> AddHit(trackHit);
    }

    fGoodHits -> MoveHitsTo(fTrackHits);
    fGoodHits -> Clear();

    return kStepNewTrack;
}

int LTHelixTrackFindingTask::StepNextPhase()
{
    if (fPhaseIndex==0) {
        fPhaseIndex = 1;

        fPadPlane -> ResetEvent();

        fTrackHits -> Clear();
        fCandHits -> Clear();
        fGoodHits -> Clear();
        ReturnBadHitsToPadPlane();

        //fMinHitsToFitInitTrack = 7;
        //fCutMinNumHitsInitTrack = 10;
        fCutMaxNumHitsInitTrack = 25;
        //fCutMinNumHitsFinalTrack = 15;
        //fCutMinHelixRadius = 30.;
        //fTrackLengthCutScale = 2.5;
        //fCutdkInExpectedTrackPath = 4.;

        return kStepNewTrack;
    }
    //else if (fPhaseIndex==1)
    return kStepEndEvent;
}

int LTHelixTrackFindingTask::StepEndEvent()
{
    fTrackArray -> Compress();

    Int_t numTracks = fTrackArray -> GetEntriesFast();
    for (Int_t iTrack = 0; iTrack < numTracks; ++iTrack) {
        LKHelixTrack *track = (LKHelixTrack *) fTrackArray -> At(iTrack);
        track -> SetTrackID(iTrack);
        track -> FinalizeHits();
    }

    lk_info << "Number of found tracks: " << fTrackArray -> GetEntries() << endl;

    return kStepEndOfEvent;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void LTHelixTrackFindingTask::ReturnBadHitsToPadPlane()
{
    fNumBadHits = fBadHits -> GetEntriesFast();
    for (Int_t iBad = 0; iBad < fNumBadHits; ++iBad)
        fPadPlane -> AddHit((LKHit *) fBadHits -> GetHit(iBad));
    fBadHits -> Clear();
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

double LTHelixTrackFindingTask::CorrelateHitWithTrackCandidate(LKHelixTrack *track, LKHit *hit)
{
#ifdef DEBUG_STEP
    lk_debug << "======================================== track hits" << endl;
    fCurrentTrack -> Print(">");
    lk_debug << "======================================== hit" << endl;
    hit -> Print();
    //hit -> GetHitArray() -> Print();
    lk_debug << "========================================" << endl;
#endif

    if (hit -> GetNumTrackCands() != 0)
        return 0;

    Double_t quality = 0;

    auto row = hit -> GetRow();
    auto layer = hit -> GetLayer();

    LKVector3 qosHit(hit -> GetPosition(), fReferenceAxis);

    auto trackHits = track -> GetHitArray();
    bool passCutdk = false;
    Int_t numTrackHits = trackHits -> GetNumHits();
#ifdef DEBUG_STEP
    lk_debug << numTrackHits << endl;
#endif
    for (Int_t iTrackHit = 0; iTrackHit < numTrackHits; ++iTrackHit) {
        auto trackHit = (LKHit *) trackHits -> GetHit(iTrackHit);
        //if (row == trackHit -> GetRow() && layer == trackHit -> GetLayer()) continue; // XXX
        LKVector3 qosTrackHit(trackHit -> GetPosition(), fReferenceAxis);
        auto di = qosHit.I() - qosTrackHit.I();
        auto dj = qosHit.J() - qosTrackHit.J();
        auto distPadCenter = sqrt(di*di + dj*dj);
        Double_t tangentDip = abs(qosTrackHit.K())/sqrt(qosTrackHit.I()*qosTrackHit.I()+qosTrackHit.J()*qosTrackHit.J());
        Double_t dkInExpectedTrackPath = 1.2 * distPadCenter * tangentDip;
        Double_t dkBetweenTwoHits = abs(qosHit.K()-qosTrackHit.K());

        if (dkInExpectedTrackPath < fCutdkInExpectedTrackPath) dkInExpectedTrackPath = fCutdkInExpectedTrackPath;
#ifdef DEBUG_STEP
        lk_debug << "TrackHit-" << iTrackHit << ": dkBetweenTwoHits and dkInExpectedTrackPath is " << dkBetweenTwoHits << " <? " << dkInExpectedTrackPath << endl;
#endif
        if (dkBetweenTwoHits < dkInExpectedTrackPath) {
            passCutdk = true;
        }
    }

    if (!passCutdk) {
        return 0;
    }

    if (track -> IsBad()) {
#ifdef DEBUG_STEP
        lk_debug << "track was bad" << endl;
#endif
        quality = 1;
    }
    else if (track -> IsLine()) {
#ifdef DEBUG_STEP
        lk_debug << "track was line" << endl;
#endif
        LKVector3 perp = track -> PerpLine(hit -> GetPosition());

        Double_t rmsCut = track -> GetRMST();
        if (rmsCut < fTrackHCutLL) rmsCut = fTrackHCutLL;
        if (rmsCut > fTrackHCutHL) rmsCut = fTrackHCutHL;
        rmsCut = 3 * rmsCut;

        if (perp.K() > rmsCut) {
            quality = 0;
        }
        else {
            perp.SetK(0);
            auto magcut = 15.;
            if (perp.Mag() < magcut)
                //if (perp.Mag() < 10*pos1.K()/sqrt(pos1.Mag()))
            {
                quality = 1;
            }
        }
    }
    else if (track -> IsPlane()) {
#ifdef DEBUG_STEP
        lk_debug << "track was plane" << endl;
#endif
        Double_t dist = (track -> PerpPlane(hit -> GetPosition())).Mag();

        Double_t rmsCut = track -> GetRMST();
        if (rmsCut < fTrackHCutLL) rmsCut = fTrackHCutLL;
        if (rmsCut > fTrackHCutHL) rmsCut = fTrackHCutHL;
        rmsCut = 3 * rmsCut;

#ifdef DEBUG_STEP
        lk_debug << "dist is " <<  dist << endl;
        lk_debug << "track rmsCut is " <<  rmsCut << "(" << track -> GetRMST() << ")" << endl;
#endif
        if (dist < rmsCut) {
            quality = 1;
        }
    }
    else {
#ifdef DEBUG_STEP
        lk_debug << "track is ..." << endl;
        track -> Print(">");
#endif
    }

    return quality;
}

double LTHelixTrackFindingTask::CorrelateHitWithTrack(LKHelixTrack *track, LKHit *hit, Double_t rScale)
{
    Double_t scale = rScale * fDefaultScale;
    Double_t trackLength = track -> TrackLength();
    if (trackLength < 500.)
        scale = scale + (500. - trackLength)/500.;

    /*
       auto direction = track->Momentum().Unit();
       Double_t dot = abs(direction.Dot(LKVector3(fReferenceAxis,0,0,1).GetXYZ()));
       auto a = 1.;
       auto b = 1.;
       if (dot > .5) {
       a = (2*(dot-.5)+1); //
       b = (2*(dot-.5)+1); //
       }
       auto trackHCutLL = a*fTrackHCutLL;
       auto trackHCutHL = a*fTrackHCutHL;
       auto trackWCutLL = b*fTrackWCutLL;
       auto trackWCutHL = b*fTrackWCutHL;
     */
    auto trackHCutLL = fTrackHCutLL;
    auto trackHCutHL = fTrackHCutHL;
    auto trackWCutLL = fTrackWCutLL;
    auto trackWCutHL = fTrackWCutHL;

    Double_t rmsWCut = track -> GetRMSR();
    if (rmsWCut < trackWCutLL) rmsWCut = trackWCutLL;
    if (rmsWCut > trackWCutHL) rmsWCut = trackWCutHL;
    rmsWCut = scale * rmsWCut;

    Double_t rmsHCut = track -> GetRMST();
    if (rmsHCut < trackHCutLL) rmsHCut = trackHCutLL;
    if (rmsHCut > trackHCutHL) rmsHCut = trackHCutHL;
    rmsHCut = scale * rmsHCut;

    TVector3 qHead = track -> Map(track -> PositionAtHead());
    TVector3 qTail = track -> Map(track -> PositionAtTail());
    TVector3 q = track -> Map(hit -> GetPosition());

    if (qHead.Z() > qTail.Z()) {
        if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, q.Z() - qHead.Z())) return 0;
        if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, qTail.Z() - q.Z())) return 0;
    } else {
        if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, q.Z() - qTail.Z())) return 0;
        if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, qHead.Z() - q.Z())) return 0;
    }

    Double_t dr = abs(q.X());
    Double_t quality = 0;
    if (dr < rmsWCut && abs(q.Y()) < rmsHCut)
        quality = sqrt((dr-rmsWCut)*(dr-rmsWCut)) / rmsWCut;

    return quality;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

bool LTHelixTrackFindingTask::BuildAndConfirmTrack(LKHelixTrack *track, bool &tailToHead)
{
    track -> SortHits(!tailToHead);
    auto trackHits = track -> GetHitArray();
    Int_t numHits = trackHits -> GetNumHits();

    TVector3 q, m;
    //Double_t lPre = track -> ExtrapolateByMap(trackHits->at(numHits-1)->GetPosition(), q, m);

    Double_t extrapolationLength = 10.;
    for (Int_t iHit = 1; iHit < numHits; iHit++) 
    {
        //auto trackHit = (LKHit *) trackHits -> GetHit(numHits-iHit-1);
        auto trackHit = (LKHit *) trackHits -> GetLastHit();
        //Double_t lCur = track -> ExtrapolateByMap(trackHit->GetPosition(), q, m);

        Double_t quality = CorrelateHitWithTrack(track, trackHit);

        if (quality <= 0) {
            track -> RemoveHit(trackHit);
            trackHit -> RemoveTrackCand(trackHit -> GetTrackID());
#ifdef DEBUG_STEP
            lk_debug << "after remove hit" << endl;
            track -> Print();
#endif
            Int_t helicity = track -> Helicity();
            track -> Fit();
            if (helicity != track -> Helicity())
                tailToHead = !tailToHead;

            continue;
        }

        /*
           Double_t dLength = abs(lCur - lPre);
           extrapolationLength = 10;
           while(dLength > 0 && AutoBuildByInterpolation(track, tailToHead, extrapolationLength, 1)) { dLength -= 10; }
         */
    }

    extrapolationLength = 0;
    while (AutoBuildByExtrapolation(track, tailToHead, extrapolationLength)) {}

    if (track -> GetNumHits() < fCutMinNumHitsFinalTrack)
        return false;

    return true;
}

bool LTHelixTrackFindingTask::AutoBuildByExtrapolation(LKHelixTrack *track, bool &buildHead, Double_t &extrapolationLength)
{
    if (track -> GetNumHits() < fCutMinNumHitsFinalTrack)
        return false;

    TVector3 p;
#ifdef DEBUG_STEP
    track -> Print();
#endif
    if (buildHead) p = track -> ExtrapolateHead(extrapolationLength);
    else           p = track -> ExtrapolateTail(extrapolationLength);

    return AutoBuildAtPosition(track, p, buildHead, extrapolationLength);
}

bool LTHelixTrackFindingTask::AutoBuildAtPosition(LKHelixTrack *track, TVector3 p, bool &tailToHead, Double_t &extrapolationLength, Double_t rScale)
{
    LKVector3 p2(p,fReferenceAxis);
    if (fPadPlane -> IsInBoundary(p2.I(), p2.J()) == false)
        return false;

    Int_t helicity = track -> Helicity();

    Double_t rms = 3*track -> GetRMSR();
    if (rms < 25) 
        rms = 25;

    Int_t range = Int_t(rms/8);
    fPadPlane -> PullOutNeighborHits(p2.I(), p2.J(), range, fCandHits);

    fNumCandHits = fCandHits -> GetEntriesFast();
    Bool_t foundHit = false;

    if (fNumCandHits != 0) 
    {
        fCandHits -> SortByCharge(false);

        for (Int_t iHit = 0; iHit < fNumCandHits; iHit++) {
            LKHit *candHit = (LKHit *) fCandHits -> GetLastHit();
            fCandHits -> RemoveLastHit();
            TVector3 pos = candHit -> GetPosition();

            Double_t quality = 0; 
            if (CheckParentTrackID(candHit) < 0) 
                quality = CorrelateHitWithTrack(track, candHit, rScale);

            if (quality > 0) {
                track -> AddHit(candHit);
                track -> Fit();
                foundHit = true;
            } else
                fBadHits -> AddHit(candHit);
        }
    }

    if (foundHit) {
        extrapolationLength = 10; 
        if (helicity != track -> Helicity())
            tailToHead = !tailToHead;
    }
    else {
        extrapolationLength += 10; 
        if (extrapolationLength > 3 * track -> TrackLength()) {
            return false;
        }
    }

    return true;
}

int LTHelixTrackFindingTask::CheckParentTrackID(LKHit *hit)
{
    vector<Int_t> *candTracks = hit -> GetTrackCandArray();
    Int_t numCandTracks = candTracks -> size();
    if (numCandTracks == 0)
        return -2;

    Int_t trackID = -1;
    for (Int_t iCand = 0; iCand < numCandTracks; ++iCand) {
        Int_t candTrackID = candTracks -> at(iCand);
        if (candTrackID != -1) {
            trackID = candTrackID;
        }
    }

    return trackID;
}

bool LTHelixTrackFindingTask::CheckTrackQuality(LKHelixTrack *track)
{
    /*
       Double_t continuity = CheckTrackContinuity(track);
       if (continuity < .6) {
       if (track -> TrackLength() * continuity < 500)
       return false;
       }
     */

    if (track -> GetHelixRadius() < fCutMinHelixRadius)
        return false;

    return true;
}

double LTHelixTrackFindingTask::CheckTrackContinuity(LKHelixTrack *track)
{
    Int_t numHits = track -> GetNumHits();
    if (numHits < 2)
        return -1;

    track -> SortHits();

    Double_t total = 0;
    Double_t continuous = 0;

    TVector3 qPre, mPre, qCur, mCur; // I need q(position on helix)

    LKVector3 kqPre;
    LKVector3 kqCur;

    auto trackHits = track -> GetHitArray();
    auto pPre = trackHits->GetHit(0)->GetPosition();
    track -> ExtrapolateByMap(pPre, qPre, mPre);
    kqPre = LKVector3(qPre,fReferenceAxis);

    auto axis1 = fPadPlane -> GetAxis1();
    auto axis2 = fPadPlane -> GetAxis2();

    for (auto iHit = 1; iHit < numHits; iHit++)
    {
        auto pCur = trackHits->GetHit(iHit)->GetPosition();
        track -> ExtrapolateByMap(pCur, qCur, mCur);

        kqCur = LKVector3(qCur,fReferenceAxis);

        auto val1 = kqCur.At(axis1) - kqPre.At(axis1);
        auto val2 = kqCur.At(axis2) - kqPre.At(axis2);

        auto length = sqrt(val1*val1 + val2*val2);

        total += length;
        if (length <= 1.2 * fPadPlane -> PadDisplacement())
            continuous += length;


        kqPre = kqCur;
    }

    return continuous/total;
}

bool LTHelixTrackFindingTask::CheckHitDistInAlphaIsLargerThanQuarterPi(LKHelixTrack *track, Double_t dLength)
{
    if (dLength > 0) {
        if (dLength > .5*track -> TrackLength()) {
            if (abs(track -> AlphaAtTravelLength(dLength)) > .5*TMath::Pi()) {
                return true;
            }
        }
    }
    return false;
}
