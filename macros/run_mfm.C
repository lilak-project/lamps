#include "LKLogger.h"

void run_mfm(int runNo=0, int split=0)
{
    lk_logger("log/lamps.log");

    auto run = new LKRun();

    run -> AddPar("config_lamps.mac");
    //run -> AddDetector(new LAMPSTPC());

    run -> SetEventTrigger(new LKMFMConversionTask());
    //run -> Add(new TTEventPreviewTask());
    //run -> Add(new TTPulseAnalysisTask());
    //run -> Add(new TTHTTrackingTask());

    run -> Init();
    run -> Print();
    run -> SetEventCountForMessage(1);
    run -> WriteExitLog("log_end");
    run -> Run();
}
