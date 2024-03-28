void run_conv()
{
    auto run = new LKRun();
    //run -> SetTag("conv");

    run -> AddPar("config_conv.mac");
    //run -> AddDetector(new LAMPSTPC());

    run -> SetEventTrigger(new LKMFMConversionTask());
    //run -> Add(new TTEventPreviewTask());
    //run -> Add(new TTPulseAnalysisTask());
    //run -> Add(new TTHTTrackingTask());

    run -> Init();
    run -> Print();
    run -> SetEventCountForMessage(1);
    run -> Run();
}
