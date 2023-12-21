#include <stdlib.h>
#include <bitset>

void TPCToF(int rn=90){//76,89,90,92
  gROOT->Reset();
  gStyle -> SetOptFit();

  cout << Form("Analysis : run%04d",rn) << endl;

  TCanvas *c1 = new TCanvas("c1", "signal", 800, 800);
  c1->Divide(2,2);
  c1->Draw();

  //histogram
  //  TFile *hist=new TFile(Form("./hbook/run%d.hbook",rn),"recreate");

  TH2F* signal0 = new TH2F("signal0","signal", 514,-1,513., 4100,0.,4100.);

  //  TH2F* padxy = new TH2F("padxy","LAMPS TPC : x vs y", 300,-600.,600., 300,-600.,600.);
  //  TH2F* padzy = new TH2F("padzy","LAMPS TPC : z vs y", 530,-10.,520., 300,-600.,600.);
  //  TH2F* padzx = new TH2F("padzx","LAMPS TPC : z vs x", 530,-10.,520., 300,-600.,600.);
  TH2F* padxy = new TH2F("padxy","LAMPS TPC : x vs y", 750,-750.,750., 750,-750.,750.);
  //TH2F* padzy = new TH2F("padzy","LAMPS TPC : z vs y", 530,-10.,520., 750,-750.,750.);
  TH2F* padzy = new TH2F("padzy","LAMPS TPC : z vs y", 1600,-250.,1350., 750,-750.,750.);
  TH2F* padzx = new TH2F("padzx","LAMPS TPC : z vs x", 530,-10.,520., 750,-750.,750.);

  //open raw(merge) file
  ifstream ifile;
  string ifilename = "/home/cens-alpha-00/data/lamps/run_0123.dat.14-12-23_14h29m55s";

  //string ifilename = Form("./irun%d.graw",rn);
  ifile.open(ifilename, ios::in | ios::binary);
  if(ifile.fail()){
    cout << "Cannot open " << ifilename << endl;
    return;
  }

  char buf1[1], buf2[2], buf3[3], buf4[4], buf6[6], buf9[9];
  int n, nItems;
  int ievt, coboIdx, asadIdx;
  int agetIdx, chanIdx, buckIdx, sample;
  int max[22][4][4][68], hpos[22][4][4][68];
  int ped[22][4][4][68], nped[22][4][4][68];
  int nch, nch0, chId[23936];
  double evttime;
  
  //PAD position database
  ifstream dfile;
  dfile.open("./LAMPS_TPC_PAD_info_11");
  int icobo, iasad, iaget, ich;
  double x[22][4][4][68], y[22][4][4][68];
  int section[22][4][4][68], npad[22][4][4][68], nlayer[22][4][4][68], nasad[22][4][4][68];
  for(int i=0; i<22; i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<4; k++){
        for(int l=0; l<68; l++){
          npad[i][j][k][l] = -1;
        }
      }
    }
  }
  for(int i=0; i<2698*8; i++){
    dfile >> icobo;
    dfile >> iasad;
    dfile >> iaget;
    dfile >> ich;
    dfile >> x[icobo][iasad][iaget][ich];
    dfile >> y[icobo][iasad][iaget][ich];

    dfile >> section[icobo][iasad][iaget][ich];
    dfile >> npad[icobo][iasad][iaget][ich];
    dfile >> nlayer[icobo][iasad][iaget][ich];
    dfile >> nasad[icobo][iasad][iaget][ich];
  }
  dfile.close();

  int cobo, asad, aget, ch, ppos;
  double pedestal, ADC;
  int nnch = 0;
  double xx[21584], yy[21584], zz[21584], aa[21584];
  int ss[21584], ll[21584];
  double zpos;

  ifstream ToFfile;
  ToFfile.open(Form("./20231124/skim_data2/run_%d.dat",rn));

  int ToFevt, nBToF, nFToF, chBToF[48], chFToF[48];
  double BToF[4][48], FToF[4][48];
  double pi = acos(-1.);
  double iphi = 2*pi/48.;
  double r = 707.5;
  double xBToF[48], yBToF[48], zBToF[48];
  double xFToF[48], yFToF[48];
  double dtdc, Bzpos;

  double Bpar[3][48];
  ifstream pfile;
  pfile.open("./20231122/BToF_position_par");
  for(int i=0; i<48; i++){
    pfile >> Bpar[0][i];
    pfile >> Bpar[1][i];
    pfile >> Bpar[2][i];
  }
  pfile.close();

  int index=-1;

  /***************** read data *****************/
  while(1){

    if(index==2) goto ToFana;

    signal0->Reset();
    padxy->Reset();
    padzy->Reset();
    padzx->Reset();
    
    nch = 0;
    ifile.seekg(20,ios::cur);
    for(int ncobo=0; ncobo<11+11; ncobo++){
      for(int nnasad=0; nnasad<4; nnasad++){
	nch0 = 0;
	int metacheck=0;
	int nmeta=0;
	//metaType
      again:
	ifile.read(buf1,1);
	n = (long)(buf1[0] & 0xff);
	metacheck=metacheck+1;
	if(metacheck>1000){
	  cout << "End by metaType!"<< endl;
	  sleep(3);
	  ifile.close();
	  //hist->Write();
	  //hist->Close();
	  return;
	}
	nmeta = nmeta + 1;
	if(n!=8){
	  //cout << Form("metaType : 0x%02X : %d : %d",n,n,metacheck) << endl;
	  /*
	  ifile.close();
	  hist->Write();
	  hist->Close();
	  return;
	  */
	  goto again;
	}
	cout << "nmeta : " << nmeta << endl;
	if(nmeta>1){
	  int bb; cin >> bb;
	}

	//frameSize
	ifile.read(buf3,3);
	n = (long)(buf3[2] & 0xff) + ((long)(buf3[1] & 0xff)<<8)
	  + ((long)(buf3[0] & 0xff)<<16);
	cout << Form("frameSize : 0x%06X : %d",n,n) << endl;
	int frameSize = n;

	//dataSource
	ifile.read(buf1,1);
	n = (long)(buf1[0] & 0xff);
	//cout << Form("dataSource : 0x%02X : %d",n,n) << endl;

	//frameType : 1(partial,zero), 2(full)
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//n = (long)(buf2[1] & 0xff);
	//cout << Form("frameType : 0x%04X : %d",n,n) << endl;

	//revision : 5
	ifile.read(buf1,1);
	n = (long)(buf1[0] & 0xff);
	//cout << Form("revision : 0x%02X : %d",n,n) << endl;

	//headerSize
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("headerSize : 0x%04X : %d",n,n) << endl;

	//itemSize
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	cout << Form("itemSize : 0x%04X : %d",n,n) << endl;

	//nItems
	ifile.read(buf4,4);
	n = (long)(buf4[3] & 0xff) + ((long)(buf4[2] & 0xff)<<8) + 
	  ((long)(buf4[1] & 0xff)<<16) + ((long)(buf4[0] & 0xff)<<24);
	cout << Form("nItems : 0x%08X : %d",n,n) << endl;
	nItems = n;

	//eventTime : 10ns -> sec
	ifile.read(buf6,6);
	evttime = (long)(buf6[5] & 0xff)/100000000. 
	  + ((long)(buf6[4] & 0xff)<<8)/100000000.
	  + ((long)(buf6[3] & 0xff)<<16)/100000000.
	  + ((long)(buf6[2] & 0xff)<<24)/100000000.
	  + ((long)(buf6[1] & 0xff)<<32)/100000000.
	  + ((long)(buf6[0] & 0xff)<<40)/100000000.;
	cout << "eventTime : " << evttime << endl;

	//eventIdx
	ifile.read(buf4,4);
	n = (long)(buf4[3] & 0xff) + ((long)(buf4[2] & 0xff)<<8) + 
	  ((long)(buf4[1] & 0xff)<<16) + ((long)(buf4[0] & 0xff)<<24);
	cout << Form("eventIdx : 0x%08X : %d",n,n) << endl;
	ievt = n;

	//coboIdx
	ifile.read(buf1,1);
	n = (long)(buf1[0] & 0xff);
	cout << Form("coboIdx : 0x%02X : %d",n,n) << endl;
	coboIdx = n;

	//asadIdx
	ifile.read(buf1,1);
	n = (long)(buf1[0] & 0xff);
	cout << Form("asadIdx : 0x%02X : %d",n,n) << endl;
	asadIdx = n;

	//readOffset
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("readOffset : 0x%04X : %d",n,n) << endl;

	//status
	ifile.read(buf1,1);
	n = (long)(buf1[0] & 0xff);
	//cout << Form("status : 0x%02X : %d",n,n) << endl;

	//hitPat_0 (ch0->ch68)
	ifile.read(buf9,9);
	cout << "hitPat_0 : ";
	int nn=-1, hitch=-4;
	for(int i=0; i<9; i++){
	  n = (long)(buf9[i] & 0xff);
	  cout << bitset<8>(n) << " ";
	  for(int j=0; j<8; j++){
	    nn = (long)(buf9[i] >> (7-j)) & 1;
	    if(nn==1 && hitch>=0){
	      cout << "(" << hitch << ") ";
	      max[coboIdx][asadIdx][0][hitch] = -100;
	      hpos[coboIdx][asadIdx][0][hitch] = -100;
	      ped[coboIdx][asadIdx][0][hitch] = 0;
	      nped[coboIdx][asadIdx][0][hitch] = 0;
	      chId[nch] = coboIdx*10000 + asadIdx*1000 + 0*100 + hitch;
	      nch = nch + 1;
	      nch0 = nch0 + 1;
	    }
	    hitch = hitch + 1;
	    nn = -1;
	  }
	}
	cout << endl;

	//hitPat_1 (ch0->ch68)
	ifile.read(buf9,9);
	cout << "hitPat_1 : ";
	hitch=-4;
	for(int i=0; i<9; i++){
	  n = (long)(buf9[i] & 0xff);
	  cout << bitset<8>(n) << " ";
	  for(int j=0; j<8; j++){
	    nn = (long)(buf9[i] >> (7-j)) & 1;
	    if(nn==1 && hitch>=0){
	      cout << "(" << hitch << ") ";
	      max[coboIdx][asadIdx][1][hitch] = -100;
	      hpos[coboIdx][asadIdx][1][hitch] = -100;
	      ped[coboIdx][asadIdx][1][hitch] = 0;
	      nped[coboIdx][asadIdx][1][hitch] = 0;
	      chId[nch] = coboIdx*10000 + asadIdx*1000 + 1*100 + hitch;
	      nch = nch + 1;
	      nch0 = nch0 + 1;
	    }
	    hitch = hitch + 1;
	    nn = -1;
	  }
	}
	cout << endl;

	//hitPat_2 (ch0->ch68)
	ifile.read(buf9,9);
	cout << "hitPat_2 : ";
	hitch=-4;
	for(int i=0; i<9; i++){
	  n = (long)(buf9[i] & 0xff);
	  cout << bitset<8>(n) << " ";
	  for(int j=0; j<8; j++){
	    nn = (long)(buf9[i] >> (7-j)) & 1;
	    if(nn==1 && hitch>=0){
	      cout << "(" << hitch << ") ";
	      max[coboIdx][asadIdx][2][hitch] = -100;
	      hpos[coboIdx][asadIdx][2][hitch] = -100;
	      ped[coboIdx][asadIdx][2][hitch] = 0;
	      nped[coboIdx][asadIdx][2][hitch] = 0;
	      chId[nch] = coboIdx*10000 + asadIdx*1000 + 2*100 + hitch;
	      nch = nch + 1;
	      nch0 = nch0 + 1;
	    }
	    hitch = hitch + 1;
	    nn = -1;
	  }
	}
	cout << endl;

	//hitPat_3 (ch0->ch68)
	ifile.read(buf9,9);
	cout << "hitPat_3 : ";
	hitch=-4;
	for(int i=0; i<9; i++){
	  n = (long)(buf9[i] & 0xff);
	  cout << bitset<8>(n) << " ";
	  for(int j=0; j<8; j++){
	    nn = (long)(buf9[i] >> (7-j)) & 1;
	    if(nn==1 && hitch>=0){
	      cout << "(" << hitch << ") ";
	      max[coboIdx][asadIdx][3][hitch] = -100;
	      hpos[coboIdx][asadIdx][3][hitch] = -100;
	      ped[coboIdx][asadIdx][3][hitch] = 0;
	      nped[coboIdx][asadIdx][3][hitch] = 0;
	      chId[nch] = coboIdx*10000 + asadIdx*1000 + 3*100 + hitch;
	      nch = nch + 1;
	      nch0 = nch0 + 1;
	    }
	    hitch = hitch + 1;
	    nn = -1;
	  }
	}
	cout << endl;
	cout << nch0 << " / " << nch << endl;

	//multip_0
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("multip_0 : 0x%04X : %d",n,n) << endl;

	//multip_1
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("multip_1 : 0x%04X : %d",n,n) << endl;

	//multip_2
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("multip_2 : 0x%04X : %d",n,n) << endl;

	//multip_3
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("multip_3 : 0x%04X : %d",n,n) << endl;

	//windowOut
	ifile.read(buf4,4);
	n = (long)(buf4[3] & 0xff) + ((long)(buf4[2] & 0xff)<<8) + 
	  ((long)(buf4[1] & 0xff)<<16) + ((long)(buf4[0] & 0xff)<<24);
	//cout << Form("windowOut : 0x%08X : %d",n,n) << endl;

	//lastCell_0
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("last_Cell_0 : 0x%04X : %d",n,n) << endl;

	//lastCell_1
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("last_Cell_1 : 0x%04X : %d",n,n) << endl;

	//lastCell_2
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("last_Cell_2 : 0x%04X : %d",n,n) << endl;

	//lastCell_3
	ifile.read(buf2,2);
	n = (long)(buf2[1] & 0xff) + ((long)(buf2[0] & 0xff)<<8);
	//cout << Form("last_Cell_3 : 0x%04X : %d",n,n) << endl;

	//skip between header and data
	ifile.seekg(169,ios::cur);
  
	//data : 4byte
	int nloop = (frameSize*256-256)/4;
	if(nloop<nItems){
	  cout << "nloop : " << nloop << " - " << nItems << endl;
	  int bb; cin >> bb;
	}
	for(int ni=0; ni<nloop; ni++){
	  //for(int ni=0; ni<nItems; ni++){
	  //for(int ni=0; ni<nch0*512; ni++){
	  ifile.read(buf4,4);
	  n = (long)(buf4[3] & 0xff) + ((long)(buf4[2] & 0xff)<<8) + 
	    ((long)(buf4[1] & 0xff)<<16) + ((long)(buf4[0] & 0xff)<<24);
	  //cout << Form("data : 0x%08X : %d",n,n) << endl;
	  agetIdx = (n>>30) & 0x3;
	  chanIdx = (n>>23) & 0x7F;
	  buckIdx = (n>>14) & 0x1FF;
	  sample = n & 0xFFF;
	  //cout << Form("%d-%d : 0x%08X : AGET %d : ch %d : bucket %d : sample : %d", 
	  //	       ni,nItems,n,agetIdx, chanIdx, buckIdx, sample) << endl;
	  signal0->Fill(buckIdx,sample);
	  
	  ped[coboIdx][asadIdx][agetIdx][chanIdx] = ped[coboIdx][asadIdx][agetIdx][chanIdx] + sample;
	  nped[coboIdx][asadIdx][agetIdx][chanIdx] = nped[coboIdx][asadIdx][agetIdx][chanIdx] + 1;
	  if(max[coboIdx][asadIdx][agetIdx][chanIdx]<sample){
	    max[coboIdx][asadIdx][agetIdx][chanIdx]=sample;
	    hpos[coboIdx][asadIdx][agetIdx][chanIdx]=buckIdx;
	  }
	}//data
	//cout << Form("data : 0x%08X : %d",n,n) << endl;
	//cout << ncobo << " " << nnasad << endl;
	//int bb; cin >> bb;
      }//for(int nnasad=0; nnasad<4; nnasad++){
    }//for(int ncobo=0; ncobo<11+11; ncobo++){






    // analysis : max, hpos, nch, chId = coboIdx*10000 + asadIdx*1000 + 3*100 + hitch;
    //    int cobo, asad, aget, ch, ppos;
    //    double pedestal, ADC;
    //    int nnch = 0;
    //    double xx[21584], yy[21584], zz[21584], aa[21584];
    //    int ss[21584], ll[21584];
    nnch = 0;
    cout << "=================================" << endl;
    for(int ich=0; ich<nch; ich++){
      cobo = int(chId[ich]/10000);
      asad = int(chId[ich]/1000) - cobo*10;
      aget = int(chId[ich]/100) - cobo*100 - asad*10;
      ch = chId[ich] - cobo*10000 - asad*1000 - aget*100;
      pedestal = ped[cobo][asad][aget][ch]/nped[cobo][asad][aget][ch];
      //ADC = max[cobo][asad][aget][ch];
      ADC = max[cobo][asad][aget][ch] - pedestal;
      if(ADC>50){
	ppos = hpos[cobo][asad][aget][ch];
	zpos = ppos*40.*61./1000.;
	cout << Form("CoBo%d-AsAd%d-AGET%d,ch%d : ADC-%f, hpos-%d : ",
		     cobo, asad, aget, ch, ADC, ppos);
	padxy->Fill(x[cobo][asad][aget][ch],y[cobo][asad][aget][ch],ADC);
	//padzy->Fill(ppos,y[cobo][asad][aget][ch],ADC);
	padzy->Fill(zpos,y[cobo][asad][aget][ch],ADC);
	padzx->Fill(ppos,x[cobo][asad][aget][ch],ADC);
	cout << x[cobo][asad][aget][ch] << " " << y[cobo][asad][aget][ch] << endl;
	xx[nnch] = x[cobo][asad][aget][ch];
	yy[nnch] = y[cobo][asad][aget][ch];
	zz[nnch] = ppos;
	aa[nnch] = ADC;
	ss[nnch] = section[cobo][asad][aget][ch];
	ll[nnch] = nlayer[cobo][asad][aget][ch];
	nnch = nnch + 1;
      }
    }
    cout << ievt << " : " << evttime << "sec : " << nnch << "/" << nch << endl;
    cout << "=================================" << endl;

  ToFana:
    if(index==1) goto drawing;
    //ToF data
    ToFfile >> ToFevt;
    cout << "Event : " << ToFevt << endl;
    ToFfile >> nBToF;
    cout << "BToF : " << nBToF << endl;
    for(int i=0; i<nBToF; i++){
      ToFfile >> chBToF[i];
      // BToF8 is (707.5,0) view from downstream
      xBToF[i] = r*cos((chBToF[i]-8)*iphi);
      yBToF[i] = r*sin((chBToF[i]-8)*iphi);
      cout << chBToF[i] << " : " << xBToF[i] << ", " << yBToF[i] << ", ";
      for(int j=0; j<4; j++){
        ToFfile >> BToF[j][i];
      }
      //double dtdc, Bzpos;
      dtdc = BToF[1][i]-BToF[3][i];
      Bzpos = Bpar[1][chBToF[i]] + dtdc*Bpar[2][chBToF[i]] + 558.4;
      zBToF[i] = Bzpos;
      cout << zBToF[i] << endl;
      cout << dtdc << " " << Bpar[1][chBToF[i]] << " " << Bpar[2][chBToF[i]] << endl;
    }
    ToFfile >> nFToF;
    cout << "FToF : " << nFToF << endl;
    for(int i=0; i<nFToF; i++){
      ToFfile >> chFToF[i];
      xFToF[i] = r/2*cos((chBToF[i]-8)*iphi);
      yFToF[i] = r/2*sin((chBToF[i]-8)*iphi);
      cout << chFToF[i] << " : " << xFToF[i] << ", " << yFToF[i] << endl;
      for(int j=0; j<4; j++){
        ToFfile >> FToF[j][i];
      }
    }
    cout << "=================================" << endl;










  drawing:
    /////////////////////////////// drawing
    //if(nnch>0){
    if(nBToF>=2){
    c1->cd(1);
    signal0->SetStats(0);
    signal0->Draw("colz");
    
    c1->cd(2);
    padxy->SetStats(0);
    padxy->Draw("colz");
    padxy->GetXaxis()->SetTitle("X (mm)");
    padxy->GetYaxis()->SetTitle("Y (mm)");
    TLine *line = new TLine();
    line->DrawLine(-208.557,503.5, 208.557,503.5);
    line->DrawLine(208.557,503.5, 503.5,208.557);
    line->DrawLine(503.5,208.557, 503.5,-208.557);
    line->DrawLine(503.5,-208.557, 208.557,-503.5);
    line->DrawLine(208.557,-503.5, -208.557,-503.5);
    line->DrawLine(-208.557,-503.5, -503.5,-208.557);
    line->DrawLine(-503.5,-208.557, -503.5,208.557);
    line->DrawLine(-503.5,208.557, -208.557,503.5);
    line->DrawLine(-208.557,503.5, -40.182,97.007);
    line->DrawLine(208.557,503.5, 40.182,97.007);
    line->DrawLine(503.5,208.557, 97.007,40.182);
    line->DrawLine(503.5,-208.557, 97.007,-40.182);
    line->DrawLine(208.557,-503.5, 40.182,-97.007);
    line->DrawLine(-208.557,-503.5, -40.182,-97.007);
    line->DrawLine(-503.5,-208.557, -97.007,-40.182);
    line->DrawLine(-503.5,208.557, -97.007,40.182);
    TArc *arc = new TArc();
    arc->SetFillStyle(0);
    arc->DrawArc(0.,0.,105.);
    arc->DrawArc(0.,0.,170./2);
    arc->DrawArc(0.,0.,r);
    for(int i=0; i<nBToF; i++){
      arc->DrawArc(xBToF[i],yBToF[i],10.);
    }
    
    c1->cd(3);
    padzx->SetStats(0);
    padzx->Draw("colz");
    padzx->GetXaxis()->SetTitle("Z (TB)");
    padzx->GetYaxis()->SetTitle("X (mm)");
    line->DrawLine(0.0,503.5, 512.0,503.5);
    line->DrawLine(0.,-503.5, 512.,-503.5);
    line->DrawLine(0.,503.5, 0.,-503.5);
    line->DrawLine(512.,503.5, 512.,-503.5);
    line->DrawLine(0.,105., 512.,105.);
    line->DrawLine(0.,-105., 512.,-105.);
    
    c1->cd(4);
    padzy->SetStats(0);
    padzy->Draw("colz");
    padzy->GetXaxis()->SetTitle("Z (TB)");
    padzy->GetYaxis()->SetTitle("Y (mm)");
    line->DrawLine(0.0,503.5, 1200.0,503.5);
    line->DrawLine(0.,-503.5, 1200.,-503.5);
    line->DrawLine(0.,503.5, 0.,-503.5);
    line->DrawLine(1200.,503.5, 1200.,-503.5);
    line->DrawLine(0.,105., 1200.,105.);
    line->DrawLine(0.,-105., 1200.,-105.);
    for(int i=0; i<nBToF; i++){
      arc->DrawArc(zBToF[i],yBToF[i],10.);
    }
    /*
    line->DrawLine(0.0,503.5, 512.0,503.5);
    line->DrawLine(0.,-503.5, 512.,-503.5);
    line->DrawLine(0.,503.5, 0.,-503.5);
    line->DrawLine(512.,503.5, 512.,-503.5);
    line->DrawLine(0.,105., 512.,105.);
    line->DrawLine(0.,-105., 512.,-105.);
    */    
    
    c1->Modified();
    c1->Update();
    
    //sleep(1);
    

    cout << "quit(0), TPC(1), ToF(2), both(3) : ";
    cin >> index;
    //index=3;
    sleep(1);
    if(index==0){
      break;
    }

    }//if(nnch>0){
  }//while(1){

  ifile.close();
  ToFfile.close();

  //hist->Write();
  //hist->Close();
}

