void readCalibration(TFile* calibrationFile,                        // file with output of modulecalibration
                    TChain* tree,                                   // TChain of events to analyze
                    TList *formulasAnalysis,                        // list of formulas
                    std::vector<Crystal_t> &crystal,               // collection of crystal structures
                    bool UseTriggerChannelCut,                      // include TriggerChannelCuthannel cut in crystalCut
                    bool UseBroadCut,                               // include broadCut in crystalCut
                    bool UseCutNoise,                               // include CutNoise in crystalCut
                    bool UsePhotopeakEnergyCut,                     // include PhotopeakEnergyCut in crystalCut
                    bool UseCutgs)                                   // include CutGs in crystalCut
                    // std::vector<detector_t> &detectorSaturation,    // collection of detector structures
                    // std::vector<int> forbidden_channels)            // channels excluded from timing correction. this can be an empty vector, and all channels will be included. The number of the channel is the number after "ch_" in the name of the "Graph Delay ch_" and "RMS Graph Delay ch_" graphs located inside the "TimeCorrection" folder of each crystal. The number corresponds to the ADC channel number, i.e. the channel number of the V1740 pair of digitizer (it can go from 0 to 63)
{




  calibrationFile->cd("Module 0.0");

  TList *listModule = gDirectory->GetListOfKeys();
  int nKeysMod = listModule->GetEntries();
  std::vector<std::string> keysModName;
  // fill a vector with the leaves names
  std::string mppc_prefix("MPPC");
  for(int i = 0 ; i < nKeysMod ; i++){
    keysModName.push_back(listModule->At(i)->GetName());
  }

  TCut *taggingPhotopeakCut;
  float marginWZgraph = 0.1;      // defaults
  float WrangeBinsForTiming = 10; // defaults
  float length = 15.0;  // defaults
  // taggingPhotopeakCut->SetName(""); //defaults
  int taggingCrystalTimingChannel = -1; //defaults
  int taggingCrystalChannel = -1; //defaults
  float taggingPosition = 0;
  std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
  std::string taggingCrystalTimingChannel_prefix("taggingCrystalTimingChannel");
  std::string taggingCrystalChannel_prefix("taggingCrystalChannel");
  std::string taggingPosition_prefix("taggingPosition");
  std::string length_prefix("length");
  std::string marginWZgraph_prefix("marginWZgraph");
  std::string WrangeBinsForTiming_prefix("WrangeBinsForTiming");
  std::vector<int> *pChannels;
  std::vector<float> *pSaturation;
  std::vector<float> *pPedestal;
  gDirectory->GetObject("channels",pChannels);
  gDirectory->GetObject("saturation",pSaturation);
  gDirectory->GetObject("pedestal",pPedestal);

  std::vector<int> DetChannels = pChannels[0];
  std::vector<float> saturation = pSaturation[0];
  std::vector<float> pedestal = pPedestal[0];

  std::vector<detector_t> detectorSaturation;

  for(unsigned int iSat = 0; iSat < DetChannels.size(); iSat++)
  {
    detector_t tempDetector;
    tempDetector.digitizerChannel = DetChannels[iSat];
    tempDetector.saturation = saturation[iSat];
    tempDetector.pedestal = pedestal[iSat];
    detectorSaturation.push_back(tempDetector);
  }



  std::vector<std::string> MPPCfolders;
  for(unsigned int i = 0 ; i < keysModName.size() ; i++)
  {
    if (!keysModName[i].compare(0, mppc_prefix.size(), mppc_prefix))
    {
      MPPCfolders.push_back(keysModName[i]);
    }
    if(!keysModName[i].compare(0,marginWZgraph_prefix.size(),marginWZgraph_prefix)) // find
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      marginWZgraph = atof(snameCh.str().c_str());
      std::cout << "margin cut set to " << marginWZgraph << std::endl;
    }
    if(!keysModName[i].compare(0,WrangeBinsForTiming_prefix.size(),WrangeBinsForTiming_prefix)) // find
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      // binningForWCut = atof(snameCh.str().c_str());
      // applyBinRestriction = true;
      // std::cout << "Applying bin restriction, binningForWCut set to " << binningForWCut << std::endl;
      WrangeBinsForTiming = atoi(snameCh.str().c_str());
    }
    if(!keysModName[i].compare(0,length_prefix.size(),length_prefix)) // find
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      // binningForWCut = atof(snameCh.str().c_str());
      // applyBinRestriction = true;
      // std::cout << "Applying bin restriction, binningForWCut set to " << binningForWCut << std::endl;
      length = atof(snameCh.str().c_str());
    }
    if(!keysModName[i].compare(0,taggingPhotopeakCut_prefix.size(),taggingPhotopeakCut_prefix)) // find tcut
    {
      taggingPhotopeakCut = (TCut*) gDirectory->Get( keysModName[i].c_str());
    }
    if(!keysModName[i].compare(0,taggingCrystalTimingChannel_prefix.size(),taggingCrystalTimingChannel_prefix)) // find tcut
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      taggingCrystalTimingChannel = atoi(snameCh.str().c_str());
    }

    if(!keysModName[i].compare(0,taggingCrystalChannel_prefix.size(),taggingCrystalChannel_prefix)) // find tcut
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      taggingCrystalChannel = atoi(snameCh.str().c_str());
    }

    if(!keysModName[i].compare(0,taggingPosition_prefix.size(),taggingPosition_prefix)) // find tcut
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      taggingPosition = atof(snameCh.str().c_str());
    }
  }



  for(unsigned int iMppc = 0 ; iMppc < MPPCfolders.size() ; iMppc++)
  {
    // std::cout << MPPCfolders[iMppc] << std::endl;
    gDirectory->cd(MPPCfolders[iMppc].c_str());
    TList *listMppc = gDirectory->GetListOfKeys();
    int nKeysMppc = listMppc->GetEntries();
    std::vector<std::string> keysMppcName;
    // fill a vector with the leaves names
    std::string crystal_prefix("Crystal");
    // std::string det_prefix("digitizerChannel");
    // std::string saturation_prefix("saturation");
    for(int i = 0 ; i < nKeysMppc ; i++){
      keysMppcName.push_back(listMppc->At(i)->GetName());
    }

    std::vector<std::string> CrystalFolders;

    for(unsigned int i = 0 ; i < keysMppcName.size() ; i++)
    {
      if (!keysMppcName[i].compare(0, crystal_prefix.size(), crystal_prefix))
      {
        CrystalFolders.push_back(keysMppcName[i]);
      }
    }

    for(unsigned int iCry = 0 ; iCry < CrystalFolders.size() ; iCry++)
    {
      //  std::cout << CrystalFolders[iCry] << std::endl;
       gDirectory->cd(CrystalFolders[iCry].c_str());

       Crystal_t temp_crystal;
       temp_crystal.number = -1;
       temp_crystal.minAcceptedW = 0;
       temp_crystal.maxAcceptedW = 1;
       temp_crystal.wMinSlicing = 0;
       temp_crystal.wMaxSlicing = 1;
       temp_crystal.wStepSlicing = 1;
       temp_crystal.marginWZgraph = marginWZgraph;
       temp_crystal.WrangeBinsForTiming = WrangeBinsForTiming;
       temp_crystal.length = length;
       temp_crystal.taggingCrystalTimingChannel = taggingCrystalTimingChannel;
       temp_crystal.taggingCrystalChannel = taggingCrystalChannel;
       temp_crystal.taggingPosition = taggingPosition;
       temp_crystal.taggingPhotopeakCut = taggingPhotopeakCut;
       temp_crystal.CrystalCut = NULL;
       temp_crystal.CrystalCutWithoutCutG = NULL;
       temp_crystal.PhotopeakEnergyCut = NULL;
       temp_crystal.calibrationGraph = NULL;
       temp_crystal.simpleCTR = NULL;
       temp_crystal.centralCTR = NULL;
       temp_crystal.allCTR = NULL;
       temp_crystal.poliCorrCTR = NULL;
       temp_crystal.simpleCTR_norm = NULL;
       temp_crystal.centralCTR_norm = NULL;
       temp_crystal.allCTR_norm = NULL;
       temp_crystal.poliCorrCTR_norm = NULL;
       temp_crystal.likeCTR = NULL;
       temp_crystal.likeCTR_norm = NULL;
       temp_crystal.hybridCTR = NULL;
       temp_crystal.hybridCTR_norm = NULL;
       temp_crystal.wz = NULL;
       temp_crystal.accepted = true;
       temp_crystal.tw_correction = NULL;
       temp_crystal.polishedCorrection = false;
       temp_crystal.detectorSaturation = detectorSaturation;
       temp_crystal.FormulaAnalysis = NULL;

       //set taggingPhotopeakCut Formula if it is found in calibration file
       if((temp_crystal.taggingCrystalTimingChannel != -1) && (temp_crystal.taggingCrystalChannel != -1) )
       {
         TTreeFormula* FormulaTagAnalysis = new TTreeFormula("FormulaTagAnalysis",
                                                             taggingPhotopeakCut->GetTitle(),
                                                             tree);
         formulasAnalysis->Add(FormulaTagAnalysis);
         temp_crystal.FormulaTagAnalysis = FormulaTagAnalysis;
       }

       //get crystal number
       temp_crystal.number = atoi((CrystalFolders[iCry].substr(crystal_prefix.size()+1,CrystalFolders[iCry].size()-crystal_prefix.size()-1)).c_str());

       TList *listCry = gDirectory->GetListOfKeys();
       int nKeysCry = listCry->GetEntries();
       std::vector<std::string> keysCryName;
       if(nKeysCry) //if directory not empty
       {

         for(int i = 0 ; i < nKeysCry ; i++){
           keysCryName.push_back(listCry->At(i)->GetName());
         }
         std::string CalibName;
         std::string CutName;
         std::vector<std::string> cutgNames;
         std::string cutG_prefix("cutg");
         std::string calibration_prefix("Calibration");
         // std::string crystalCut_prefix("CrystalCut");
         // std::string crystalCutWithoutCutG_prefix("CrystalCutWithoutCutG");
         // std::string photopeakEnergyCut_prefix("PhotopeakEnergyCut");
         std::string channel_prefix("digitizerChannel");
         std::string w_channels_prefix("channelsNumRelevantForW");
         std::string timing_channel_prefix("timingChannel");
         std::string wz_prefix("w(z)");
         std::string t_channels_for_poli_mean_prefix("tChannelsForPolishedCorrectionMean");
         std::string t_channels_for_poli_FWHM_prefix("tChannelsForPolishedCorrectionFWHM");
         std::string mean_for_poli_prefix("meanForPolishedCorrection");
         std::string fwhm_for_poli_prefix("fwhmForPolishedCorrection");

         std::string TriggerChannelCut_prefix  ("TriggerChannelCut");
         std::string broadCut_prefix           ("broadCut");
         std::string CutNoise_prefix           ("CutNoise");
         std::string PhotopeakEnergyCut_prefix ("PhotopeakEnergyCut");

         std::string lightCentral_prefix("Light collected in trigger crystal");
         std::string lightAll_prefix("Sum spectrum highlighted");
         std::string basicCTR_prefix("Basic CTR histogram");

         bool dirExists      = true;
         bool dirDelayExists = true;
         bool dirRMSExists   = true;
         std::string TimeCorrection_prefix("TimeCorrection");
         std::string Delay_prefix         ("DelayDir");
         std::string RMS_prefix           ("RMSDir");

         for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
         {
           if(!keysCryName[i].compare(0,calibration_prefix.size(),calibration_prefix)) //find calibration graph
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCanvas* C_graph = NULL;
             TGraph *calibGraph = NULL;
             C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
             if(C_graph)
               calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
             if(calibGraph)
               temp_crystal.calibrationGraph = calibGraph;
           }

           if(!keysCryName[i].compare(0,wz_prefix.size(),wz_prefix)) //find calibration graph
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCanvas* C_graph = NULL;
             TGraph *calibGraph = NULL;
             C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
             if(C_graph)
               calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
             if(calibGraph)
               temp_crystal.wz = calibGraph;
           }

           if(!keysCryName[i].compare(0,TriggerChannelCut_prefix.size(),TriggerChannelCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.TriggerChannelCut = cut;
           }

           if(!keysCryName[i].compare(0,broadCut_prefix.size(),broadCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.broadCut = cut;
           }
           if(!keysCryName[i].compare(0,CutNoise_prefix.size(),CutNoise_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.CutNoise = cut;
           }
           if(!keysCryName[i].compare(0,PhotopeakEnergyCut_prefix.size(),PhotopeakEnergyCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.PhotopeakEnergyCut = cut;
           }



           // if(!keysCryName[i].compare(0,photopeakEnergyCut_prefix.size(),photopeakEnergyCut_prefix)) // find tcut
           // {
           //  //  std::cout << keysCryName[i] << std::endl;
           //   TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
           //   if(cut)
           //     temp_crystal.PhotopeakEnergyCut = cut;
           // }

           if(!keysCryName[i].compare(0,cutG_prefix.size(),cutG_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCutG* cutg = (TCutG*) gDirectory->Get(keysCryName[i].c_str());

             temp_crystal.cutg.push_back(cutg);
           }

           if(!keysCryName[i].compare(0,channel_prefix.size(),channel_prefix)) // find detector channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.detectorChannel = atoi(snameCh.str().c_str());
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }

           if(!keysCryName[i].compare(0,lightCentral_prefix.size(),lightCentral_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.lightCentralHisto = aHisto;
           }

           if(!keysCryName[i].compare(0,lightAll_prefix.size(),lightAll_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.lightAllHisto = aHisto;
           }

           if(!keysCryName[i].compare(0,basicCTR_prefix.size(),basicCTR_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.basicCTRhisto = aHisto;
           }


           if(!keysCryName[i].compare(0,timing_channel_prefix.size(),timing_channel_prefix)) // find timing channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.timingChannel = atoi(snameCh.str().c_str());
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }

           if(!keysCryName[i].compare(0,w_channels_prefix.size(),w_channels_prefix)) // find detector channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             std::vector<int> *v;
             gDirectory->GetObject("channelsNumRelevantForW",v);
             // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.relevantForW = v[0];
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }

           if(!keysCryName[i].compare(0,t_channels_for_poli_mean_prefix.size(),t_channels_for_poli_mean_prefix))
           {
             std::vector<int> *v;
             gDirectory->GetObject("tChannelsForPolishedCorrectionMean",v);
             temp_crystal.polishedCorrection = true;
             temp_crystal.tChannelsForPolishedCorrectionMean = v[0];
           }

           if(!keysCryName[i].compare(0,t_channels_for_poli_FWHM_prefix.size(),t_channels_for_poli_FWHM_prefix))
           {
             std::vector<int> *v;
             gDirectory->GetObject("tChannelsForPolishedCorrectionFWHM",v);
             temp_crystal.polishedCorrection = true;
             temp_crystal.tChannelsForPolishedCorrectionFWHM = v[0];
           }

           if(!keysCryName[i].compare(0,mean_for_poli_prefix.size(),mean_for_poli_prefix))
           {
             std::vector<double> *v;
             gDirectory->GetObject("meanForPolishedCorrection",v);
             temp_crystal.meanForPolishedCorrection = v[0];
           }

           if(!keysCryName[i].compare(0,fwhm_for_poli_prefix.size(),fwhm_for_poli_prefix))
           {
             std::vector<double> *v;
             gDirectory->GetObject("fwhmForPolishedCorrection",v);
             temp_crystal.fwhmForPolishedCorrection = v[0];
           }

           // if(!keysCryName[i].compare(0,TimeCorrection_prefix.size(),TimeCorrection_prefix))
           // {
           //   dirExists = true;
           // }
           // if(!keysCryName[i].compare(0,Delay_prefix.size(),Delay_prefix))
           // {
           //   dirDelayExists = true;
           // }
           // if(!keysCryName[i].compare(0,RMS_prefix.size(),RMS_prefix))
           // {
           //   dirRMSExists = true;
           // }



         }


         std::stringstream sname;
         // dirExists =
         if(dirExists)
         {
           gDirectory->cd("TimeCorrection");
           sname.str("");
           // sname << "Central correction - Crystal " << temp_crystal.number;
           // temp_crystal.centralCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           // sname.str("");
           // sname << "Full correction - Crystal " << temp_crystal.number;
           // temp_crystal.allCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           // sname.str("");
           //
           // if(likelihood)
           // {
           //   sname.str("");
           //   sname << "Likelihood correction - Crystal " << temp_crystal.number;
           //   temp_crystal.likeCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           //   sname.str("");
           // }
           // if(hybridCorrection)
           // {
           //   sname.str("");
           //   sname << "Hybrid correction - Crystal " << temp_crystal.number;
           //   temp_crystal.hybridCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           //   sname.str("");
           // }

           temp_crystal.path = gDirectory->GetPath();
           TList *listTcorr = gDirectory->GetListOfKeys();
           int nKeysTcorr = listTcorr->GetEntries();
           std::vector<std::string> keysTcorrName;
           if(nKeysTcorr) //if directory not empty
           {
             for(int i = 0 ; i < nKeysTcorr ; i++){
               keysTcorrName.push_back(listTcorr->At(i)->GetName());
             }

             std::string deltaWGraph_prefix = "DeltaW Graph";
             std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
             std::string graph_delay_prefix = "Graph Delay ch_";
             std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
             std::string delay_timing_ch_prefix = "delayTimingChannels";
             for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
             {
               if(!keysTcorrName[i].compare(0,deltaWGraph_prefix.size(),deltaWGraph_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   temp_crystal.tw_correction = calibGraph;
                   //fit with straight line
                   TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                   calibGraph->Fit(line,"Q");
                   temp_crystal.tw_correction_line = line;

                 }

               }

               if(!keysTcorrName[i].compare(0,rms_deltaWGraph_prefix.size(),rms_deltaWGraph_prefix))
               {
                TGraphErrors *calibGraph = NULL;
                calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                if(calibGraph)
                {
                  temp_crystal.rms_tw_correction = calibGraph;
                  //fit with straight line
                  TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                  calibGraph->Fit(line,"Q");
                  temp_crystal.rms_tw_correction_line = line;

                }

               }

               if(!keysTcorrName[i].compare(0,graph_delay_prefix.size(),graph_delay_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   // -- check if the channel is not excluded by the user
                   // extract next 2 characters
                   //
                   std::string str2 = keysTcorrName[i].substr (graph_delay_prefix.size(),6);     // take a string with next 6 characters after the prefix
                   std::size_t found = str2.find_first_of("_");                                  // find next "_"
                   std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                   // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                   int current_ch = atoi(str3.c_str());                                          // transform in int
                   // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output

                   bool acceptCh = true;

                   // forbidden_channels part is removed
                   // for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                   // {
                   //   if(current_ch == forbidden_channels[iForb])
                   //   {
                   //     acceptCh = false;
                   //   }
                   // }

                   if(acceptCh)               // add graph if the ch is accepted
                   {
                     temp_crystal.delay.push_back(calibGraph);
                     //fit with straight line
                     TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                     calibGraph->Fit(line,"Q");
                     temp_crystal.delay_line.push_back(line);

                   }



                 }

               }

               if(!keysTcorrName[i].compare(0,rms_graph_delay_prefix.size(),rms_graph_delay_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());


                 if(calibGraph)
                 {
                   // -- check if the channel is not excluded by the user
                   // extract next 2 characters
                   //
                   std::string str2 = keysTcorrName[i].substr (rms_graph_delay_prefix.size(),6);     // take a string with next 6 characters after the prefix
                   std::size_t found = str2.find_first_of("_");                                  // find next "_"
                   std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                   int current_ch = atoi(str3.c_str());                                          // transform in int
                   // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output


                   bool acceptCh = true;

                   // forbidden_channels part is removed
                   // for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                   // {
                   //   if(current_ch == forbidden_channels[iForb])
                   //   {
                   //     acceptCh = false;
                   //   }
                   // }

                   if(acceptCh)               // add graph if the ch is accepted
                   {
                     temp_crystal.rms_delay.push_back(calibGraph);
                     //fit with straight line
                     TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                     calibGraph->Fit(line,"Q");
                     temp_crystal.rms_delay_line.push_back(line);
                   }
                 }

               }

               if(!keysTcorrName[i].compare(0,delay_timing_ch_prefix.size(),delay_timing_ch_prefix))
               {

                 // // -- check if the channel is not excluded by the user
                 // // extract next 2 characters
                 // //
                 // std::string str2 = keysTcorrName[i].substr (delay_timing_ch_prefix.size(),6);     // take a string with next 6 characters after the prefix
                 // std::size_t found = str2.find_first_of("_");                                  // find next "_"
                 // std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                 // // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                 // int current_ch = atoi(str3.c_str());                                          // transform in int
                 // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output
                 //
                 // bool acceptCh = true;
                 // for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                 // {
                 //   if(current_ch == forbidden_channels[iForb])
                 //   {
                 //     acceptCh = false;
                 //   }
                 // }
                 //
                 // if(acceptCh)               // add graph if the ch is accepted
                 // {
                   //  std::cout << keysCryName[i] << std::endl;
                   std::stringstream snameCh;
                   std::vector<int> *v;
                   gDirectory->GetObject("delayTimingChannels",v);
                   // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
                   //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                  //  istringstream()
                   temp_crystal.delayTimingChannels = v[0];
                   //  std::cout <<temp_crystal.detectorChannel << std::endl;
                  //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                            //  << temp_crystal.detectorChannel << std::endl;
                 // }

               }
             }
           }
           gDirectory->cd("..");
         }

         // bool dirDelayExists;
         sname.str("");
         // dirDelayExists =
         if(dirDelayExists)
         {
           gDirectory->cd("DelayDir");
           TList *listTcorr = gDirectory->GetListOfKeys();
           int nKeysTcorr = listTcorr->GetEntries();
           std::vector<std::string> keysTcorrName;
           if(nKeysTcorr) //if directory not empty
           {
             for(int i = 0 ; i < nKeysTcorr ; i++){
               keysTcorrName.push_back(listTcorr->At(i)->GetName());
             }

             // std::string deltaWGraph_prefix = "DeltaW Graph";
             std::string graph_prefix = "Delay graph t";
             // std::string graph_delay_prefix = "Graph Delay ch_";
             // std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
             // std::string delay_timing_ch_prefix = "delayTimingChannels";
             for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
             {
               if(!keysTcorrName[i].compare(0,graph_prefix.size(),graph_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   graphs_t temp_graph;
                   std::string rmsName = calibGraph->GetName();

                   std::size_t found = rmsName.find_last_of(graph_prefix);
                   std::string tChannelStringFromRms   = rmsName.substr(found  +1);

                   int ch   = atoi( tChannelStringFromRms.c_str() );

                   temp_graph.timingChannel = ch;
                   temp_graph.graph = calibGraph;

                   temp_crystal.delay_graphs.push_back(temp_graph);
                 }
               }
             }
           }
           gDirectory->cd("..");
         }



         // bool dirRMSExists;
         sname.str("");
         // dirRMSExists =
         if(dirRMSExists)
         {
           gDirectory->cd("RMSdir");
           TList *listTcorr = gDirectory->GetListOfKeys();
           int nKeysTcorr = listTcorr->GetEntries();
           std::vector<std::string> keysTcorrName;
           if(nKeysTcorr) //if directory not empty
           {
             for(int i = 0 ; i < nKeysTcorr ; i++){
               keysTcorrName.push_back(listTcorr->At(i)->GetName());
             }

             // std::string deltaWGraph_prefix = "DeltaW Graph";
             std::string graph_prefix = "RMS graph t";
             // std::string graph_delay_prefix = "Graph Delay ch_";
             // std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
             // std::string delay_timing_ch_prefix = "delayTimingChannels";
             for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
             {
               if(!keysTcorrName[i].compare(0,graph_prefix.size(),graph_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   graphs_t temp_graph;
                   std::string rmsName = calibGraph->GetName();

                   std::size_t found = rmsName.find_last_of(graph_prefix);
                   std::string tChannelStringFromRms   = rmsName.substr(found  +1);

                   int rmsCh   = atoi( tChannelStringFromRms.c_str() );

                   std::cout << rmsCh << std::endl;

                   temp_graph.timingChannel = rmsCh;
                   temp_graph.graph = calibGraph;

                   temp_crystal.rms_graphs.push_back(temp_graph);

                 }
               }
             }
           }
           gDirectory->cd("..");
         }


         sname.str("");


         // build and save the crystal cut
         TCut crystalCut ;

         if(UseTriggerChannelCut)       crystalCut += temp_crystal.TriggerChannelCut->GetTitle();
         if(UseBroadCut)                crystalCut += temp_crystal.broadCut->GetTitle();
         if(UseCutNoise)                crystalCut += temp_crystal.CutNoise->GetTitle();
         if(UsePhotopeakEnergyCut)      crystalCut += temp_crystal.PhotopeakEnergyCut->GetTitle();
         if(UseCutgs)
         {
           crystalCut += temp_crystal.cutg[0]->GetName();
           crystalCut += temp_crystal.cutg[1]->GetName();
         }

         TTreeFormula* FormulaAnalysis = new TTreeFormula("FormulaAnalysis",crystalCut,tree);
         temp_crystal.FormulaAnalysis = FormulaAnalysis;

         // TList* formulasAnalysis = new TList();
         formulasAnalysis->Add(FormulaAnalysis);

         if(temp_crystal.FormulaAnalysis)
         {
           crystal.push_back(temp_crystal);
         }
       }
       gDirectory->cd("..");
    }
    calibrationFile->cd("Module 0.0");
  }
}



int setWandZcuts(std::vector<Crystal_t> &crystal)
{
  //------------------------------------------------//
  //------------------------------------------------//
  //------------------------------------------------//

  // the w of the crystal is sliced in modulecalibration, by doing
  // zMin = marginWZgraph           --> calc beginW using w(z) tgraph
  // zMax = length - marginWZgraph  --> calc endW using w(z) tgraph
  // then the range  (endW - beginW) is divided by WrangeBinsForTiming
  // so each bin is long (endW - beginW)/WrangeBinsForTiming
  // starts from
  // wmin = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming)
  // and ends in
  // wmax = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming)
  // then on this bin we calc whatever we are calculating, and the values are
  // associated to the mean w of the win, simply
  // wmean = (wmax + wmin) / 2.0
  // as a consequence all the tgraphs can be used only from the first wmean to the las wmean,
  // whatever is outside this range has to be thrown away
  // this has to be written in the crystal as a limit on the acceptance of w

  std::cout << std::endl;
  std::cout << "|--------------------------------------------|" << std::endl;
  std::cout << "|                 w and z cuts               |" << std::endl;
  std::cout << "|--------------------------------------------|" << std::endl;
  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {
    if(crystal[iCry].accepted)
    {
      float beginW = crystal[iCry].wz->Eval(crystal[iCry].length - crystal[iCry].marginWZgraph);
      float endW = crystal[iCry].wz->Eval(crystal[iCry].marginWZgraph);


      // first wmean
      Float_t wmin,wmin2;
      Float_t wmax,wmax2;
      Float_t wmean,wmean2;

      int iBin = 0;

      wmin = beginW + ((iBin*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmax = beginW + (((iBin+1)*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmean = (wmax + wmin) / 2.0;
      crystal[iCry].wStepSlicing = (wmax - wmin);
      crystal[iCry].minAcceptedW = wmean;
      crystal[iCry].wMinSlicing = wmin;

      iBin = crystal[iCry].WrangeBinsForTiming -1;
      wmin2 = beginW + ((iBin*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmax2 = beginW + (((iBin+1)*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmean2 = (wmax2 + wmin2) / 2.0;
      crystal[iCry].maxAcceptedW = wmean2;
      crystal[iCry].wMaxSlicing = wmax2;


      //DEBUG
      std::cout << "Crystal = " << crystal[iCry].number << std::endl;
      std::cout << "Step lenght in w, for slicing = " << crystal[iCry].wStepSlicing << std::endl;
      std::cout << "Steps in w, for slicing = " << crystal[iCry].WrangeBinsForTiming << std::endl;
      std::cout << "min accepted w in slicing = " << crystal[iCry].wMinSlicing
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].wMinSlicing)
                << std::endl;
      std::cout << "min accepted w in Eval = " << crystal[iCry].minAcceptedW
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].minAcceptedW)
                << std::endl;
      std::cout << "max accepted w in slicing = " << crystal[iCry].wMaxSlicing
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].wMaxSlicing)
                << std::endl;
      std::cout << "max accepted w in Eval = " << crystal[iCry].maxAcceptedW
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].maxAcceptedW)
                << std::endl;

      //
      // fix the delay and rms graphs
      // rms are always one more than delay (because there is no delay for central channel, but there is a weight, the rms of the basic ctr), so run on rms vector
      for(unsigned int iRms = 0; iRms < crystal[iCry].rms_graphs.size(); iRms++ )
      {
        correction_graphs_t temp_corr;
        if(crystal[iCry].rms_graphs[iRms].timingChannel == crystal[iCry].timingChannel)
        {
          temp_corr.isMainChannel = true;
          temp_corr.delay = NULL;
        }
        else
        {
          temp_corr.isMainChannel = false;
          for(unsigned int iDel = 0; iDel < crystal[iCry].delay_graphs.size(); iDel++)
          {
            if(crystal[iCry].rms_graphs[iRms].timingChannel == crystal[iCry].delay_graphs[iDel].timingChannel)
            {
              temp_corr.delay = crystal[iCry].delay_graphs[iDel].graph;
            }

          }

        }
        temp_corr.timingChannel = crystal[iCry].rms_graphs[iRms].timingChannel;
        temp_corr.rms = crystal[iCry].rms_graphs[iRms].graph;
        crystal[iCry].correction_graphs.push_back(temp_corr);
      }
      // std::cout << crystal[iCry].delay_graphs.size() << std::endl;
      // std::cout << crystal[iCry].rms_graphs.size() << std::end
      std::cout << "Correction graphs " << crystal[iCry].correction_graphs.size() << std::endl;

      // set the polished correction vectors
      // look in the mean vectors, see if there is correspondance in the rms one, if yes add to the vector
      for(unsigned int iPoli = 0; iPoli < crystal[iCry].tChannelsForPolishedCorrectionMean.size(); iPoli++)
      {
        for(unsigned int jPoli = 0; jPoli < crystal[iCry].tChannelsForPolishedCorrectionFWHM.size(); jPoli++)
        {
          if(crystal[iCry].tChannelsForPolishedCorrectionMean[iPoli] == crystal[iCry].tChannelsForPolishedCorrectionFWHM[jPoli])
          {
            polished_correction_t tempPoli;
            tempPoli.timingChannel = crystal[iCry].tChannelsForPolishedCorrectionMean[iPoli];
            tempPoli.mean          = crystal[iCry].meanForPolishedCorrection[iPoli];
            tempPoli.rms           = crystal[iCry].fwhmForPolishedCorrection[jPoli];
            crystal[iCry].polished_correction.push_back(tempPoli);
          }
        }
      }


    }
  }
  std::cout << "|--------------------------------------------|" << std::endl;
  std::cout << std::endl;

  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {
    if(crystal[iCry].accepted)
    {
      std::cout << "crystal " << crystal[iCry].number << std::endl;
      for(unsigned int iCor = 0; iCor < crystal[iCry].correction_graphs.size(); iCor++)
      {
        std::cout << crystal[iCry].correction_graphs[iCor].isMainChannel << " "
                  << crystal[iCry].correction_graphs[iCor].timingChannel << " ";
        if(crystal[iCry].correction_graphs[iCor].isMainChannel != true)
        {
          std::cout << crystal[iCry].correction_graphs[iCor].delay->GetName() << " ";
        }

        std::cout << crystal[iCry].correction_graphs[iCor].rms->GetName() << " "
                  << std::endl;
      }
      for(unsigned int iPoli = 0; iPoli < crystal[iCry].polished_correction.size(); iPoli++)
      {
        std::cout << crystal[iCry].polished_correction[iPoli].timingChannel << " "
                  << crystal[iCry].polished_correction[iPoli].mean << " "
                  << crystal[iCry].polished_correction[iPoli].rms << std::endl;
      }



    }
  }
  // return 0;
}
