#include "LKLogger.h"

void read()
{
    auto win = LKWindowManager::GetWindowManager();

    auto run = new LKRun();
    run -> AddInputFile("data/lamps_0000.0.all.root");
    run -> Init();
    auto array = run -> GetBranchA("RawData");
    auto numEvents = run -> GetNumEvents();
    if (numEvents>5)
        numEvents = 5;

    for (auto event=0; event<numEvents; ++event)
    {
        run -> GetEntry(event);
        auto numChannels = array -> GetEntries();
        if (numChannels==0)
            continue;
        e_cout << endl;
        e_info << "event: " << event << " (" << numChannels << ")" << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel) {
            auto channel = (GETChannel*) array -> At(iChannel);
            e_cout << channel -> GetCAAC() << " ";
        }
        e_cout << endl;

        TCanvas *cvs;
        int countNewCanvases = 0;
        int numPadX = 5;
        int numPadY = 4;
        int numPads = numPadX * numPadY;
        int countPad = numPads+1;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            if (countPad==numPads+1) {
                countPad = 1;
                auto name = Form("cvs_%d_%d",event,countNewCanvases++);
                cvs = win -> CanvasFullRatio(name,name,0.8);
                cvs -> Divide(numPadX,numPadY);
            }
            auto channel = (GETChannel*) array -> At(iChannel);
            cvs -> cd(countPad++);
            auto hist = channel -> GetHist(Form("hist_%d_%d",event,iChannel));
            hist -> SetMinimum(200);
            hist -> SetMaximum(600);
            hist -> Draw();
        }
    }
}
