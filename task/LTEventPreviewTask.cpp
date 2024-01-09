#include "LTEventPreviewTask.h"
#include "LKEventHeader.h"
//#include "LTEventHeader.h"
#include "GETChannel.h"

ClassImp(LTEventPreviewTask);

LTEventPreviewTask::LTEventPreviewTask()
{
    fName = "LTEventPreviewTask";
}

bool LTEventPreviewTask::Init()
{
    lk_info << "Initializing LTEventPreviewTask" << std::endl;

    //fDetector = (TexAT2 *) fRun -> GetDetector();
    fPadPlane = (LTPadPlane*) fRun -> GetDetectorPlane();
    //fEventHeaderArray = fRun -> RegisterBranchA("EventHeader", "LTEventHeader", 1);
    fFrameHeaderArray = fRun -> GetBranchA("FrameHeader");
    fChannelArray = fRun -> GetBranchA("RawData");

    return true;
}

void LTEventPreviewTask::Exec(Option_t *option)
{
    auto numChannels = fChannelArray -> GetEntries();
    for(int iChannel=0; iChannel<numChannels; iChannel++)
    {
        auto channel = (GETChannel*) fChannelArray -> At(iChannel);
        auto cobo = channel -> GetCobo();
        auto asad = channel -> GetAsad();
        auto aget = channel -> GetAget();
        auto chan = channel -> GetChan();
        auto pad = fPadPlane -> GetPad(cobo, asad, aget, chan);
        if (pad==nullptr)
            continue;
        channel -> SetChan2(pad->GetPadID());
    }

    auto frameHeader = (LKEventHeader *) fFrameHeaderArray -> At(0);
    //frameHeader -> SetIsGoodEvent(frameHeader->IsGoodEvent());

    lk_info << endl;
}

bool LTEventPreviewTask::EndOfRun()
{
    return true;
}
