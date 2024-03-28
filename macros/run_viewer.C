void run_viewer()
{
    auto run = new LKRun();
    run -> AddPar("config_viewer.mac");
    run -> AddInputFile("lamps_0123.storage0123.root");
    run -> Add(new LKGETChannelViewer);
    run -> Init();
    run -> ExecuteEvent(1);
}

