LKPulseShapeAnalysisTask *fPSA;

void next(int eventID=-1) {
    if (eventID>=0)
        LKRun::GetRun()->ExecuteEvent(eventID);
    else
        LKRun::GetRun()->ExecuteNextEvent();
}

void drawPSA()
{
    auto chana = fPSA -> GetChannelAnalyzer();
}

void run_eve()
{
    auto run = new LKRun();
    run -> SetTag("eve");
    run -> AddPar("config_lamps.mac");
    run -> AddInputFile("data/lamps_0000.0.all.root");
    run -> AddDetector(new LAMPSTPC());
    fPSA = new LKPulseShapeAnalysisTask;
    run -> Add(fPSA);
    run -> Add(new LKEveTask);
    run -> Init();

    run -> SetEventCountForMessage(1);

    //LKWindowManager::GetWindowManager() -> FixCanvasPosition();
    //next(8);
    next();
}
