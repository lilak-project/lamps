void draw_padplane()
{
    auto pp = new LTPadPlane();
    pp -> AddPar("config_lamps.mac");//CreateParameterContainer();
    pp -> Init();
    pp -> GetPar() -> AddLine("!LTPadPlane/position_map /Users/ejungwoo/lilak/lamps_2023/common/LAMPS_TPC_PAD_info_11");
    LKWindowManager::GetWindowManager() -> CanvasSquare("cvs");
    auto hist = pp -> GetHist();
    //auto hist = pp -> GetHist(5);
    hist -> Draw();
    //pp -> Draw();

    if (1) {
        auto graph = new TGraph();
        graph -> SetMarkerColor(kBlue);
        graph -> SetMarkerStyle(20);
        graph -> SetMarkerSize(0.3);
        int cobo; // 0-21
        int slot; // 0-3
        int aget; // 0-3
        int channelID; // 0-67
        double x;
        double y;
        int section; // 1-8
        int padID; // 0-2697
        int layer; // 0-41
        int asad; // 1-11
        pp -> GetPar() -> Print();
        //TString mapFileName = pp -> GetPar() -> GetParString("LTPadPlane/position_map");
        TString mapFileName = "../common/LAMPS_TPC_PAD_info_11";
        //pp -> GetPar() -> GetParString("LTPadPlane/position_map");
        TString mapFileNameNew = mapFileName + "_update";
        cout << mapFileName << endl;
        cout << mapFileNameNew << endl;
        ifstream fileMap(mapFileName);
        ofstream fileMapNew(mapFileNameNew);
        for (auto i=0; i<21584; ++i) {
            fileMap >> cobo >> slot >> aget >> channelID >> x >> y >> section >> padID >> layer >> asad;
            int bin = hist -> FindBin(x,y) - 1;
            //cout << bin << endl;
            if (bin<0) {
                double r = sqrt(x*x + y*y);
                TVector3 dir(-x,-y,0);
                dir = dir.Unit();
                x = x+3*dir.X();
                y = y+3*dir.Y();
                bin = hist -> FindBin(x,y) - 1;
            }
            //if (bin<0) {
            //    x = x + 1;//0.5;
            //    bin = hist -> FindBin(x,y) - 5;
            //}
            //if (bin<0) {
            //    graph -> SetPoint(graph->GetN(),x,y);
            //    cout << cobo<<" "<<slot<<" "<<aget<<" "<<channelID<<" "<<x<<" "<<y <<" "<<section<<" "<<padID<<" "<<layer<<" "<<asad<<" "<<bin << endl;
            //} else
                fileMapNew << cobo<<" "<<slot<<" "<<aget<<" "<<channelID<<" "<<x<<" "<<y <<" "<<section<<" "<<padID<<" "<<layer<<" "<<asad<<" "<<bin << endl;
        }
        graph -> Draw("samep");
    }
}
