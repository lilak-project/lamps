void next() { LKRun::GetRun() -> RunSelectedEvent("EventHeader[0].IsGoodEvent()"); }

void run_eve()
{
    auto run = new LKRun();
    run -> SetTag("eve");
    run -> AddPar("config_lamps.mac");
    run -> AddInputFile("data/lamps_0000.0.all.root");
    //run -> AddInputFile("data/texat_0801.all.root");
    //run -> AddInputFile("/home/cens-alpha-00/data/texat/reco/texat_0801.reco.root");
    //run -> AddFriend("~/data/texat/reco/texat_0801.conv.root");
    run -> AddDetector(new LAMPSTPC());
    run -> Add(new LKEveTask);
    run -> Init();

    LKWindowManager::GetWindowManager() -> FixCanvasPosition();
    next();
}
