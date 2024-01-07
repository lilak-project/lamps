void next(int eventID=-1) {
    if (eventID>=0)
        LKRun::GetRun()->ExecuteEvent(eventID);
    else
        LKRun::GetRun()->ExecuteNextEvent();
}

void run_eve()
{
    auto run = new LKRun();
    run -> SetTag("eve");
    run -> AddPar("config_lamps.mac");
    run -> AddInputFile("data/lamps_0000.0.all.root");
    run -> Add(new LKPulseShapeAnalysisTask);
    run -> Add(new LKEveTask);
    run -> AddDetector(new LAMPSTPC());
    run -> Init();

    next();
}
