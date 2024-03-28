void next(int eventID=-1)
{
    if (eventID==0) LKRun::GetRun()->ExecuteNextEvent();
    if (eventID>=0) LKRun::GetRun()->ExecuteEvent(eventID);
    else LKRun::GetRun() -> RunSelectedEvent("@RawData.GetEntries()>0");
}

#include "LKLogger.h"

void run_eve()
{
    lk_logger("log/eve.log");
    auto run = new LKRun();
    run -> AddPar("config_eve.mac");
    run -> AddDetector(new LAMPSTPC());
    run -> Add(new LKPulseShapeAnalysisTask);
    run -> Add(new LTHelixTrackFindingTask);
    run -> Add(new LKEveTask);
    run -> Init();

    next(0);
    //next(15);
    //next(39);
}
