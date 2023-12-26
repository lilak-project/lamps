#ifndef LTEVENTPREVIEWTASK_HH
#define LTEVENTPREVIEWTASK_HH

#include "TClonesArray.h"
#include "LKLogger.h"
#include "LKParameterContainer.h"
#include "LKRun.h"
#include "LKTask.h"
#include "LAMPSTPC.h"
#include "LTPadPlane.h"
#include <iostream>

class LTEventPreviewTask : public LKTask
{
    public:
        LTEventPreviewTask();
        virtual ~LTEventPreviewTask() {}

        bool Init();
        void Exec(Option_t *option="");
        bool EndOfRun();

    private:
        //LAMPSTPC *fDetector = nullptr;
        LTPadPlane* fPadPlane = nullptr;

        TClonesArray* fFrameHeaderArray = nullptr;
        TClonesArray* fEventHeaderArray = nullptr;
        TClonesArray* fChannelArray = nullptr;

    ClassDef(LTEventPreviewTask,2);
};

#endif
