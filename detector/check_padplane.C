void check_padplane()
{
    auto pp = new LHPadPlaneYCut();
    pp -> CreateParameterContainer();
    pp -> Init();
    LKWindowManager::GetWindowManager() -> CanvasSquare("cvs");
    auto hist = pp -> GetHist(0);
    hist -> Draw();
}
