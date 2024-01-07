#include "LKLogger.h"

void run_pulse_extraction()
{
    lk_logger("log/lamps_pulse_extraction.log");

    auto run = new LKRun();

    run -> AddPar("config_lamps.mac");

    run -> SetTag("pulse_extraction");
    run -> AddPar("config_lamps.mac");
    run -> AddInputFile("data/lamps_0000.0.all.root");
    run -> Add(new LKPulseExtractionTask);
    run -> Init();

    run -> SetEventCountForMessage(1);
    run -> Run();
}
