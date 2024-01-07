int iChannel = 0;
TClonesArray* channelArray;
LKChannelAnalyzer* ana;

void next()
{
    if (iChannel==channelArray->GetEntries())
        return;
    auto channel = (GETChannel*) channelArray -> At(iChannel++);
    ana -> Analyze(channel->GetWaveformY());
    ana -> Draw();
}

void test_psa()
{
    auto run = new LKRun();
    run -> SetTag("psa_test");
    run -> AddPar("config_lamps.mac");
    run -> AddInputFile("data/lamps_0000.0.all.root");
    run -> AddDetector(new LAMPSTPC());
    run -> Init();
    run -> GetEvent(8);

    channelArray = run -> GetBranchA("RawData");

    auto pp = run -> GetDetectorPlane();
    ana = pp -> GetChannelAnalyzer();
    auto cvs = LKWindowManager::GetWindowManager() -> CanvasDefault("cvs",0.8);
    next();
}
