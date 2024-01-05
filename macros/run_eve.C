void next(int eventID=-1) {
    if (eventID>=0)
        LKRun::GetRun()->ExecuteEvent(eventID);
    else
        LKRun::GetRun()->ExecuteNextEvent();
}

#include "LKLogger.h"

void run_eve()
{
    lk_logger("log/eve.log");
    auto run = new LKRun();
    run -> SetTag("eve");
    run -> AddPar("config_lamps.mac");
    run -> AddInputFile("data/lamps_0000.0.all.root");
    run -> Add(new LKPulseShapeAnalysisTask);
    run -> Add(new LTHelixTrackFindingTask);
    run -> Add(new LKEveTask);
    run -> AddDetector(new LAMPSTPC());
    run -> Init();

    next();
    //next(8);
    //next(39);
}
