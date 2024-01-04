void draw_padplane()
{
    auto pp = new LTPadPlane();
    pp -> AddPar("config_lamps.mac");
    pp -> Init();
    pp -> GetPar() -> AddLine("!LTPadPlane/position_map /Users/ejungwoo/lilak/lamps_2023/common/LAMPS_TPC_PAD_info_11");
    LKWindowManager::GetWindowManager() -> CanvasSquare("cvs");
    auto hist = pp -> GetHist();
    hist -> Draw();

    if (0) {
        auto graph = new TGraph();
        graph -> SetMarkerColor(kBlue);
        graph -> SetMarkerStyle(20);
        graph -> SetMarkerSize(0.3);
        double x;
        double y;
        int cobo;      // 0-21
        int slot;      // 0-3
        int aget;      // 0-3
        int channelID; // 0-67
        int section;   // 1-8
        int padID;     // 0-2697
        int layer;     // 0-41
        int asad;      // 1-11
        pp -> GetPar() -> Print();
        TString mapFileName = "../common/LAMPS_TPC_PAD_info_11";
        TString mapFileNameNew = mapFileName + "_update";
        cout << mapFileName << endl;
        cout << mapFileNameNew << endl;
        ifstream fileMap(mapFileName);
        ofstream fileMapNew(mapFileNameNew);
        for (auto i=0; i<21584; ++i) {
            fileMap >> cobo >> slot >> aget >> channelID >> x >> y >> section >> padID >> layer >> asad;
            int bin = hist -> FindBin(x,y) - 1;
            if (bin<0) {
                double r = sqrt(x*x + y*y);
                TVector3 dir(-x,-y,0);
                dir = dir.Unit();
                x = x+3*dir.X();
                y = y+3*dir.Y();
                bin = hist -> FindBin(x,y) - 1;
            }
            fileMapNew << cobo<<" "<<slot<<" "<<aget<<" "<<channelID<<" "<<x<<" "<<y <<" "<<section<<" "<<padID<<" "<<layer<<" "<<asad<<" "<<bin << endl;
        }
        graph -> Draw("samep");
    }
}
