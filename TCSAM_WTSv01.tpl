//Bering Sea Tanner crab model
//********
//--Rec devs indexing changed so that rec_devs(y) enter n-at-size in year y. Previously, rec_devs(y) entered
//      population in year y+1.
//--Changed fishery-related devs so last year is endyr-1, not endyr. This is consistent with year being "survey year". Assessment in 2013
//      uses 2013 survey but last possible fishery year is 2012/13 (so 2012 by survey year).
//--20140419: changed to read model control file name from model configuration, rather than hard-wired name
//--20140425: changed output of selTCFR to retFcn in a number of places
//--20140506: added parameter jitter functionality
//--20140523: updated to use wtsADMB library.
//--20140602: updated to use writeParameter and jitterParameter functions in wtsADMB library.
//--20140805: decreased lower bound on log_sel50_dev_3 (directed fishery log-scale selectivity devs) from -0.5 to -1.0
//--20140814: changed bounds on log_sel50_dev_3 (directed fishery log-scale selectivity devs) from [-1.0,0.5] to [-5,5]
//--20140821: put prior on log_sel50_dev3, implemented decrease with in nll weighting on fishing-related devs 
//--20140821: added flag (doPenRed) to control file to reduce penalties on fishing-related devs with phase. 
//--20140823: added FmRKF phase and log_sel50_dev_3 bounds to control file as inputs, 
//            changed jitter factor on devs to 0.1*fac
//--20140830: corrected penalty reduction for pF_DevsGTF, corrected objfOut for log_sel50_dev_3 and pFmTCF
//            by multiplying nll by llw (was applied correctly in objfun calc)
//--20140915: converted descending limb z50 selectivity parameter for males in the snow crab fishery
//            from arithmetic (with no constraint to be > ascending limb z50) to ln-scale offset from 
//            ascending limb selectivity parameter
//--20150317: updated ModeData.cpp, FisheryData.cpp, this file to use newest version of wtsADMB rFunctions.
//--20150318: added option to assume lognormal error distributions for fishery catch data
//--20150319: added variables 'smlValFsh' and smlValSrv to centralize setting small constants used in lognormal NLL calcs
//--20150319: added some additional output to Jack's .R output file to improve plotting in R
//--20150531: Adding gmacs fishing mortality model as option (optFM). 
//            Revising so fc's are fishery capture rates, fm's are fishery mortality rates.
//            Fits are STILL to retained + discard MORTALITY (same as TCSAM2013), not catches
//--20150615: Corrected use of wts::jitterParameter() functions for change in wtsADMB library.
//--20150817: 1. Replaced hard-wired prior mean, variance (0.88, 0.05) used for srv3_q and srv3_qFem priors
//               with srv3_qPriorMean, Var and srv3_qFemPriorMean, Var, with values read in from the
//               control file.
//            2. Changed q_prior_switch to srv3_qPriorWgt and srv3_qFemPriorWgt to have separate
//               switches for male, female q's. In addition, these also function as likelihood multipliers (weights).
//            3. Changed ...femQ... to ...qFem...
//            4. Now reading phsM from control file to control estimation phase for natural mortality components
//            5. Now reading initial values for Mmult_imat, Mmultm, Mmultf, mat_big from control file
//--20150823: 1. Added switch to control asymptotic behavior of logistic selectivity functions
//--20150901: 1. Corrected when asymptotic behavior was applied to survey selectivities. Previous
//              implementation inadvertently set Q = 1.
//--20160225: 1. Reconfigured population dynamics in PROCEDURE_SECTION into separate function (runPopMod).
//            2. Added call to runPopMod() in PRELIMINARY_CALCS section and output to init files to have
//               output with initial parameter settings available.
//--20160225: 1. Added options for ln-scale female offset parameters to fishing mortality/capture rates. 
//--20160229: 1. Turned off setting TCF devs to 0.00001 in initialization section.
//            2. Added model "version".
//            3. Removed initial value for TCF devs (was 0.05)
//--20160316: 1. Started adding options for normalizing extended size comps (optPrNatZ_GTF)
//            2. Started renaming lots of variables
//--20160317: 1. Completed options for normalizing GTF
//            2. Added verModelControlFile (20160317) to check on model control file compatibility
//--20160321: 1. Added zLegal and iZLegal to calculate legal crab numbers/biomass
//            2. Revised calculations in MiscOutput to vectorize calculations, correct some formulas
//            3. Added output for z-scores from likelihood calculations to TCSAM_WTS.final.R
//--20160323: 1. Changed TCSAM_WTS.oldstyle.R to TCSAM_OLDSTYLE.R.
//            2. Added optTCFMfit to model control file to set option to fit 
//                  total (0) or discard (1) mortality for TCF males
//--20160324: 1. Added recruitment estimation info to control file.
//            2. Incremented versions to 20160324
//IMPORTANT: 2013-09 assessment model had RKC params for 1992+ discard mortality TURNED OFF. 
//           THE ESTIMATION PHASE FOR RKC DISCARD MORTALITY IS NOW SET IN THE CONTROLLER FILE!
//
//
//********
//to run mcmc 
//scmysr2003bayes -nox -mcmc 1000000 -mcsave 200
// then have to run  scmysr2003bayes -mceval to get output
//whatever is in sd report file will have a distribution and output will go to eval.csv
//for whatever have written to post later in program in the mcmc function part
//
// ===============================================================================
// ===============================================================================
GLOBALS_SECTION
    #include <math.h>
    #include <time.h>
    #include <admodel.h>
    #include "wtsADMB.hpp" 
    #include "ModelConstants.hpp"
    #include "ModelConfiguration.hpp"
    #include "FisheryData.hpp"
    #include "ModelData.hpp"
    
    adstring version = "20160324";//model version
    int verModelControlFile = 20160324;//model control file version
    
    double zLegal = 128;//current (2015/16) legal size
    int iZLegal = 0;    //index into size bins for legal size

    //model objects
    ModelConfiguration*  ptrMC;      //ptr to model configuration object
    ModelDatasets* ptrMDS;           //ptr to model datasets object

    //file streams
    ofstream mcmc;
    ofstream R_out;
    ofstream CheckFile;
    ofstream post("eval.csv");
    ofstream echo;                 //stream to echo model inputs
    
    //misc. flags
    int usePin = 0;//flag that a pin file is being used
    
    //jitter stuff
    int jitter = 0; //jitter is off
    int iSeed = -1;//default random number generator seed
    random_number_generator rng(iSeed);//random number generator
    
    //strings
    adstring fnConfigFile;//configuration file
    
    int reclag = 5;       //default lag from fertilization to recruitment (yrs)
    
    double convLBStoG   = 0.00220462262;//conversion from lbs to g
    double convLBStoMT  = 2204.62262;   //conversion from lbs to mt
    double convMLBStoMT = 2.20462262;   //conversion from 10^6 lbs to mt
    
    const int nSXs = 2;//number of sexes
    const int nSCs = 2;//number of  shell 
    const int nMSs = 2;//number of  maturity states
    
    const int nFsh = 4;//number of fisheries
    const int iTCF = 1;//index for TCF in fisheries-related arrays
    const int iSCF = 2;//index for SCF
    const int iRKF = 3;//index for RKF
    const int iGTF = 4;//index for GTF

    //debug flags
    int debugModelConfig     = 0;
    int debugModelDatasets   = 0;
    int debugModelParams     = 0;
    
    int debugDATA_SECTION    = 0;
    
    //strings to help writing to R
    adstring dmX;
    adstring dmS;
    adstring dmM;
    adstring dmY;
    adstring dmYm1;  
    adstring dmZ;
    
    int NUM_LEN_LIKE = 12; //number of likelihood components from size compositions
    int NUM_FOUT     = 38; //total number of likelihood components
    
    adstring_array strFOUT(1,NUM_FOUT);//names for objfOut components
    
    time_t start,finish;
    long hour,minute,second;
    double elapsed_time;

// ===============================================================================
// ===============================================================================
DATA_SECTION
  
 LOCAL_CALCS
    echo.open("EchoOut.dat", ios::trunc);
    echo<<"#Starting TCSAM2013 Code"<<endl;
    echo<<"#model version = "<<version<<endl;
    echo<<"#Starting DATA_SECTION"<<endl;
    cout<<"#Starting TCSAM2013 Code"<<endl;
    cout<<"#model version = "<<version<<endl;
    cout<<"#Starting DATA_SECTION"<<endl;
 END_CALCS
 
 LOCAL_CALCS
    int kf = 1;
    strFOUT(kf++) = "recruitment penalty";
    strFOUT(kf++) = "sex ratio penalty";
    strFOUT(kf++) = "immatures natural mortality penalty";
    strFOUT(kf++) = "mature male natural mortality penalty";
    strFOUT(kf++) = "mature female natural mortality penalty";

    strFOUT(kf++) = "survey q penalty";
    strFOUT(kf++) = "female survey q penalty";
    strFOUT(kf++) = "prior on female growth parameter a";
    strFOUT(kf++) = "prior on female growth parameter b";
    strFOUT(kf++) = "prior on male growth parameter a";
    strFOUT(kf++) = "prior on male growth parameter b";
    strFOUT(kf++) = "smoothing penalty on female maturity curve";
    strFOUT(kf++) = "smoothing penalty on male maturity curve";
    strFOUT(kf++) = "1st difference penalty on changes in male size at 50% selectivity in directed fishery";
    strFOUT(kf++) = "penalty on F-devs in directed fishery";
    strFOUT(kf++) = "penalty on F-devs in snow crab fishery";
    strFOUT(kf++) = "penalty on F-devs in BBRKC fishery";
    strFOUT(kf++) = "penalty on F-devs in groundfish fishery";

    strFOUT(kf++) = "likelihood for  directed fishery: retained males";
    strFOUT(kf++) = "likelihood for  directed fishery: total males";
    strFOUT(kf++) = "likelihood for  directed fishery: discarded females";
    strFOUT(kf++) = "likelihood for  snow crab fishery: discarded males";
    strFOUT(kf++) = "likelihood for  snow crab fishery: discarded females";
    strFOUT(kf++) = "likelihood for  BBRKC fishery: discarded males";
    strFOUT(kf++) = "likelihood for  BBRKC fishery: discarded females";
    strFOUT(kf++) = "likelihood for  groundfish fishery";
    strFOUT(kf++) = "likelihood for  survey: immature males";
    strFOUT(kf++) = "likelihood for  survey: mature males";
    strFOUT(kf++) = "likelihood for  survey: immature females";
    strFOUT(kf++) = "likelihood for  survey: mature females";

    strFOUT(kf++) = "likelihood for survey: mature survey biomass";
    strFOUT(kf++) = "likelihood for directed fishery: male retained catch biomass";
    strFOUT(kf++) = "likelihood for directed fishery: male total catch biomass";
    strFOUT(kf++) = "likelihood for directed fishery: female catch biomass";
    strFOUT(kf++) = "likelihood for snow crab fishery: total catch biomass";
    strFOUT(kf++) = "likelihood for BBRKC fishery: total catch biomass";
    strFOUT(kf++) = "likelihood for groundfish fishery: total catch biomass";

    strFOUT(kf++) = "penalty on sel50 devs for TCF";
 END_CALCS
 
 LOCAL_CALCS  
    //process command line options
    int on = 0;
    int flg = 0;
    //pin file use
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin=1;
        echo<<"#pin file used"<<endl;
        flg = 1;
    }
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-bpin"))>-1) {
        usePin=1;
        echo<<"#pin file used"<<endl;
        flg = 1;
    }
    //recruitment lag
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-lag"))>-1) {
        reclag=atoi(ad_comm::argv[on+1]);
        echo<<"#assumed lag for recruitment changed to: "<<reclag<<endl;
    }
    //configFile
    fnConfigFile = "TCSAM2013_ModelConfig.dat";//default model config filename
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-configFile"))>-1) {
        fnConfigFile = ad_comm::argv[on+1];
        echo<<"#config file changed to '"<<fnConfigFile<<"'"<<endl;
    }
    //debugDATA_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugDATA_SECTION"))>-1) {
        debugDATA_SECTION=1;
        echo<<"#debugDATA_SECTION turned ON"<<endl;
        flg = 1;
    }
    //debugModelDatasets
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelDatasets"))>-1) {
        debugModelDatasets=1;
        echo<<"#debugModelDatasets turned ON"<<endl;
        flg = 1;
    }
    //jitter
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-jitter"))>-1) {
        jitter=1;
        iSeed=(long)start;
        rng.reinitialize(iSeed);
        echo<<"#jitter turned ON"<<endl;
        echo<<iSeed<<"  #iSeed"<<endl;
        flg = 1;
    }
    //iSeed
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-iSeed"))>-1) {
        iSeed=atoi(ad_comm::argv[on+1]);
        rng.reinitialize(iSeed);
        echo<<iSeed<<"  #iSeed"<<endl;
        flg = 1;
    }
 END_CALCS  
    int asmtYr; //assessment year
    int mnYr;   //min model year
    int mxYr;   //max model year  (generally assessment year and last survey year. final fishery year is endyr-1.)
    int nZBs; //number of model size bins
 LOCAL_CALCS
    echo<<"#-----------------------------------"<<endl;
    echo<<"#-----------------------------------"<<endl;
    echo<<"#Reading configuration file '"<<fnConfigFile<<"'"<<endl;
    ad_comm::change_datafile_name(fnConfigFile);
    ptrMC = new ModelConfiguration();
    ptrMC->read(*(ad_comm::global_datafile));
    echo<<"#------------------ModelConfiguration-----------------"<<endl;
    echo<<(*ptrMC);
    echo<<"#-----------------------------------"<<endl;
    echo<<"#----finished model configuration---"<<endl;
    echo<<"#-----------------------------------"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelConfiguration-----------------"<<endl;
        cout<<(*ptrMC);
        cout<<"#-----------------------------------"<<endl;
        cout<<"#----finished model configuration---"<<endl;
        cout<<"#-----------------------------------"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    
    asmtYr  = ptrMC->asmtYr;
    mnYr    = ptrMC->mnYr;
    mxYr    = ptrMC->mxYr; //final pop. model numbers at size given for July 1, mxYr
    nZBs    = ptrMC->nZBs;
 END_CALCS   
    vector zBins(1,nZBs)
 LOCAL_CALCS
    zBins  = ptrMC->zBins;
    //determine starting bin for legal-sized crab
    for (int z=1;z<=nZBs;z++){if (zBins(z)<zLegal) iZLegal++;}
    echo<<"iZLegal = "<<iZLegal<<", where zLegal = "<<zLegal<<endl;
 END_CALCS   
    
    //read model data
 LOCAL_CALCS
    echo<<"#-----------------------------------"<<endl;
    echo<<"#Reading datasets file '"<<ptrMC->fnMDS<<"'"<<endl;
    if (debugModelDatasets) {
        ModelDatasets::debug=1;
        BioData::debug=1;
        TrawlSurveyData::debug=1;
        RetainedFisheryData::debug=1;
        DiscardFisheryData::debug=1;
        GroundfishTrawlFisheryData::debug=1;
    }
    ptrMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrMDS->read(*(ad_comm::global_datafile));
    if (debugModelDatasets) {
        cout<<"enter 0 to exit debug mode (1 to continue in debug mode) : ";
        cin>>debugModelDatasets;
        if (debugModelDatasets<0) exit(1);
        ModelDatasets::debug=debugModelDatasets;
        BioData::debug=debugModelDatasets;
        TrawlSurveyData::debug=debugModelDatasets;
        RetainedFisheryData::debug=debugModelDatasets;
        DiscardFisheryData::debug=debugModelDatasets;
        GroundfishTrawlFisheryData::debug=debugModelDatasets;
    }
    echo<<"#------------------ModeDatasets-----------------"<<endl;
    echo<<(*ptrMDS);
    echo<<"#-----------------------------------"<<endl;
    echo<<"#----finished model datasets---"<<endl;
    echo<<"#-----------------------------------"<<endl;
    echo<<"#------------------ModeDatasets-----------------"<<endl;
    ptrMDS->writeToR(echo,adstring("model.data")); echo<<endl;
    echo<<"#-----------------------------------"<<endl;
    echo<<"#----finished model datasets---"<<endl;
    echo<<"#-----------------------------------"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelDatasets-----------------"<<endl;
        cout<<(*ptrMDS);
        cout<<"#-----------------------------------"<<endl;
        cout<<"#----finished model datasets---"<<endl;
        cout<<"#-----------------------------------"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
 END_CALCS
    
    number smlValSrv       //small value for survey-based bulk lognormal NLL calc.s (2014: was 0.000001)
    !!smlValSrv = 0.000001;
    number smlValFsh       //small value for fishery-based bulk lognormal NLL calc.s (new 2015)
    !!smlValFsh = 0.0001;
 
    number p_const
    !!p_const=0.001;
    int call_no;
    !! call_no = 0;
    
    number spmo;   // spawning month
    !! spmo=0.;    // spmo=deviation in fraction of year from time of fishery to mating 
                   // so, if spmo=0, mating and fishery are concurrent. if spmo=2/12, mating occurs 2 months after fishery
                   
    //old data file inputs
    !!CheckFile<<"#model version = "<<version<<endl;
    int styr; // start year of the model
    int endyr;// end year of the model (generally assessment year and last year of survey data. final fishery year is endyr-1)
    !!styr  = ptrMC->mnYr;
    !!endyr = ptrMC->mxYr;
    !!CheckFile<<"styr  = "<<styr<<endl;
    !!CheckFile<<"endyr = "<<endyr<<endl;
    
    // Data stuff only from here
    //retained size freq.s in the directed tanner crab fishery
    int nObsRetZCsTCF                                   // number of years of directed fishery retained size comps
    !!nObsRetZCsTCF = ptrMDS->pTCFR->nyNatZ;
    ivector yrsObsRetZCsTCF_n(1,nObsRetZCsTCF)           // years when have directed fishery retained size comps
    matrix ssRetZCsTCF_sn(1,ALL_SHELL,1,nObsRetZCsTCF)   // sample sizes for directed fishery retained size comps (by shell condition,year)
 LOCAL_CALCS
    CheckFile<<"nZBs          = "<<nZBs<<endl;
    CheckFile<<"nObsRetZCsTCF = "<<nObsRetZCsTCF<<endl;
    yrsObsRetZCsTCF_n = ptrMDS->pTCFR->yrsNatZ;
    ssRetZCsTCF_sn    = ptrMDS->pTCFR->ssNatZ_sy;
    CheckFile<<"yrsObsRetZCsTCF_n  = "<<yrsObsRetZCsTCF_n<<endl;
    CheckFile<<"ssRetZCsTCF_sn "<<endl;
    for (int s=1;s<=ALL_SHELL;s++) CheckFile<<"sc = "<<s<<tb<<ssRetZCsTCF_sn(s)<<endl;
 END_CALCS   
 
    // ========================
    // tanner crab bycatch catch info
    int nObsDscTCF                          // number of years of directed fishery female and male discard catch data
    !!nObsDscTCF = ptrMDS->pTCFD->nyCatch;
    ivector yrsObsDscTCF_n(1,nObsDscTCF)      // years which have directed fishery discard catch data (has 0's when no fishery)
 LOCAL_CALCS
    CheckFile<<"nObsDscTCF = "<<nObsDscTCF<<endl;
    yrsObsDscTCF_n = ptrMDS->pTCFD->yrsCatch;
    CheckFile<<"yrsObsDscTCF_n  = "<<yrsObsDscTCF_n<<endl;
 END_CALCS   
    
    // tanner crab bycatch size frequency info
    int nObsZCsTCFF                                   // number of years of directed fishery female discard size comps
    !!nObsZCsTCFF = ptrMDS->pTCFD->nyNatZ;
    ivector yrsObsZCsTCFF_n(1,nObsZCsTCFF)             // years which have directed fishery female discard size comps
    vector ssZCsTCFF_n(1,nObsZCsTCFF)             // sample sizes for directed fishery female discard size comps 
    int nObsZCsTCFM                                   // number of years of directed fishery male discard size comps
    !!nObsZCsTCFM = ptrMDS->pTCFD->nyNatZ;
    ivector yrsObsZCsTCFM_n(1,nObsZCsTCFM)             // years which have directed fishery male discard size comps
    matrix ssTotZCsTCFM_sn(1,ALL_SHELL,1,nObsZCsTCFM) // sample sizes for directed fishery male discard length comps (by shell condition, year)
 LOCAL_CALCS
    CheckFile<<"nObsZCsTCFF   = "<<nObsZCsTCFF<<endl;
    yrsObsZCsTCFF_n = ptrMDS->pTCFD->yrsNatZ;
    ssZCsTCFF_n  = ptrMDS->pTCFD->ssNatZ_xsy(FEMALE,ALL_SHELL);
    CheckFile<<"yrsObsZCsTCFF_n = "<<yrsObsZCsTCFF_n<<endl;
    CheckFile<<"ssZCsTCFF_n  = "<<ssZCsTCFF_n<<endl;
    
    CheckFile<<"nObsZCsTCFM   = "<<nObsZCsTCFM<<endl;
    yrsObsZCsTCFM_n = ptrMDS->pTCFD->yrsNatZ;
    for (int s=1;s<=ALL_SHELL;s++) ssTotZCsTCFM_sn(s) = ptrMDS->pTCFD->ssNatZ_xsy(MALE,s);
    CheckFile<<"yrsObsZCsTCFM_n = "<<yrsObsZCsTCFM_n<<endl;
    CheckFile<<"ssTotZCsTCFM_sn "<<endl;
    for (int s=1;s<=ALL_SHELL;s++) CheckFile<<"sc = "<<s<<tb<<ssTotZCsTCFM_sn(s)<<endl;
 END_CALCS   
    
    // snow crab bycatch size frequency info
    int nObsZCsSCF                                  // number of years of snow crab fishery female bycatch length data
    !!nObsZCsSCF = ptrMDS->pSCF->nyNatZ;
    ivector yrsObsZCsSCF_n(1,nObsZCsSCF)             // years which have snow crab fishery bycatch length data
    vector ssZCsSCFF_n(1,nObsZCsSCF)              // sample sizes for snow crab fishery female bycatch size comps 
    matrix ssZCsSCFM_sn(1,ALL_SHELL,1,nObsZCsSCF) // sample sizes for snow crab fishery   male bycatch size comps (by shell condition, year)
 LOCAL_CALCS
    CheckFile<<"nObsZCsSCF     = "<<nObsZCsSCF<<endl;
    yrsObsZCsSCF_n = ptrMDS->pSCF->yrsNatZ;
    ssZCsSCFF_n = ptrMDS->pSCF->ssNatZ_xsy(FEMALE,ALL_SHELL);
    for (int s=1;s<=ALL_SHELL;s++) ssZCsSCFM_sn(s) = ptrMDS->pSCF->ssNatZ_xsy(MALE,s);
    CheckFile<<"yrsObsZCsSCF_n = "<<yrsObsZCsSCF_n<<endl;
    CheckFile<<"ssZCsSCFF_n = "<<ssZCsSCFF_n<<endl;
    CheckFile<<"ssZCsSCFM_sn "<<endl;
    for (int s=1;s<=ALL_SHELL;s++) CheckFile<<"sc = "<<s<<tb<<ssZCsSCFM_sn(s)<<endl;
 END_CALCS   
    
    //red king crab bycatch size frequency info
    int nObsZCsRKF                                  // number of years of BBRKC fishery bycatch size comps
    !!nObsZCsRKF = ptrMDS->pRKF->nyNatZ;
    ivector yrsObsZCsRKF_n(1,nObsZCsRKF)             // years which have BBRKC fishery bycatch size comps
    vector ssZCsRKFF_n(1,nObsZCsRKF)              // sample sizes for BBRKC fishery female bycatch size comps 
    matrix ssZCsRKFM_sn(1,ALL_SHELL,1,nObsZCsRKF) // sample sizes for BBRKC fishery   male bycatch size comps (by year and new/old shell)
 LOCAL_CALCS
    CheckFile<<"nObsZCsRKF     = "<<nObsZCsRKF<<endl;
    yrsObsZCsRKF_n = ptrMDS->pRKF->yrsNatZ;
    ssZCsRKFF_n = ptrMDS->pRKF->ssNatZ_xsy(FEMALE,ALL_SHELL);
    for (int s=1;s<=ALL_SHELL;s++) ssZCsRKFM_sn(s) = ptrMDS->pRKF->ssNatZ_xsy(MALE,s);
    CheckFile<<"yrsObsZCsRKF_n = "<<yrsObsZCsRKF_n<<endl;
    CheckFile<<"ssZCsRKFF_n = "<<ssZCsRKFF_n<<endl;
    CheckFile<<"ssZCsRKFM_sn "<<endl;
    for (int s=1;s<=ALL_SHELL;s++) CheckFile<<"sc = "<<s<<tb<<ssZCsRKFM_sn(s)<<endl;
 END_CALCS   
    
    //snow crab bycatch weight
    int nObsDscSCF                          // number of years of male and female bycatch weight data
    !!nObsDscSCF = ptrMDS->pSCF->nyCatch;
    ivector yrsObsDscSCF(1,nObsDscSCF) // years which have male and female bycatch weight data
 LOCAL_CALCS
    CheckFile<<"nObsDscSCF = "<<nObsDscSCF<<endl;
    yrsObsDscSCF = ptrMDS->pSCF->yrsCatch;
    CheckFile<<"yrsObsDscSCF "<<yrsObsDscSCF<<endl;
 END_CALCS
    
    //BBRKC bycatch weight
    int nObsDscRKF                          // number of years of male and female bycatch weight data
    !!nObsDscRKF = ptrMDS->pRKF->nyCatch;
    ivector yrsObsDscRKF(1,nObsDscRKF) // years which have male and female bycatch weight data
 LOCAL_CALCS
    CheckFile<<"nObsDscRKF = "<<nObsDscRKF<<endl;
    yrsObsDscRKF = ptrMDS->pRKF->yrsCatch;
    CheckFile<<"yrsObsDscRKF "<<yrsObsDscRKF<<endl;
 END_CALCS
    
    //groundfish trawl bycatch 
    int nObsDscGTF                        // number of years of trawl bycatch
    !!nObsDscGTF = ptrMDS->pGTF->nyCatch;
    ivector yrsObsDscGTF(1,nObsDscGTF)    // years which have trawl bycatch data
 LOCAL_CALCS
    CheckFile<<"nObsDscGTF  = "<<nObsDscGTF<<endl;
    yrsObsDscGTF = ptrMDS->pGTF->yrsCatch;
    CheckFile<<"yrsObsDscGTF "<<yrsObsDscGTF<<endl;
 END_CALCS
    
    //groundfish trawl bycatch size freq.s
    int nObsZCsGTF                             // number of years of trawl bycatch length comps
    !!nObsZCsGTF = ptrMDS->pGTF->nyNatZ;
    ivector yrsObsZCsGTF(1,nObsZCsGTF)         // years which have trawl bycatch length data
    matrix ssZCsGTF_xn(1,nSXs,1,nObsZCsGTF) // sample sizes for trawl bycatch length comps (by year and sex)
    vector ssObsZCsGTF_n(1,nObsZCsGTF)         // combined sample sizes for trawl bycatch length comps (by year)
 LOCAL_CALCS
    CheckFile<<"nObsZCsGTF     = "<<nObsZCsGTF<<endl;
    yrsObsZCsGTF = ptrMDS->pGTF->yrsNatZ;
    for (int x=1;x<=nSXs;x++) ssZCsGTF_xn(x) = ptrMDS->pGTF->ssNatZ_xsy(x,ALL_SHELL);
    ssObsZCsGTF_n.initialize();
    for (int n=1;n<=nObsZCsGTF;n++) {
        for (int x=1;x<=nSXs;x++) ssObsZCsGTF_n(n) += ssZCsGTF_xn(x,n);
    }
    CheckFile<<"yrsObsZCsGTF = "<<yrsObsZCsGTF<<endl;
    CheckFile<<"ssZCsGTF_xn  = "<<endl<<ssZCsGTF_xn<<endl;
    CheckFile<<"ssObsZcsGTF_x   = "<<endl<<ssObsZCsGTF_n<<endl;
 END_CALCS   
    
    //survey data
    int nObsSrvBio                        // number of years of survey biomass data
    !!nObsSrvBio = ptrMDS->pTSD->nyAbund;
    ivector yrsObsSrvBio_n(1,nObsSrvBio)  // years which have survey biomass estimates
    int nObsZCsSrv                                              // number of years of survey size comps
    !!nObsZCsSrv = ptrMDS->pTSD->nyNatZ;
    ivector yrsObsZCsSrv_n(1,nObsZCsSrv)                        // years which have survey size comps
    4darray ssObsZCsSrv_msxn(1,nMSs,1,nSCs,1,nSXs,1,nObsZCsSrv) // sample sizes for size comps by maturity, shell condition,sex,year
 LOCAL_CALCS
    CheckFile<<"nObsSrvBio = "<<nObsSrvBio<<endl;
    yrsObsSrvBio_n = ptrMDS->pTSD->yrsAbund;
    yrsObsZCsSrv_n = ptrMDS->pTSD->yrsNatZ;
    for (int m=1;m<=nMSs;m++) {
        for (int s=1;s<=nSCs;s++) {
            for (int x=1;x<=nSXs;x++) ssObsZCsSrv_msxn(m,s,x) = ptrMDS->pTSD->ssNatZ_xsmy(x,s,m);
        }       
    }
    CheckFile<<"yrsObsZCsSrv_n = "<<endl<<tb<<yrsObsZCsSrv_n<<endl;
    CheckFile<<"ssObsZCsSrv_msxn "<<endl;
    CheckFile<<"IMM,NS,FEMALE = "<<ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,FEMALE)<<endl;
    CheckFile<<"IMM,NS,  MALE = "<<ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,  MALE)<<endl;
    CheckFile<<"IMM,OS,FEMALE = "<<ssObsZCsSrv_msxn(IMMATURE,OLD_SHELL,FEMALE)<<endl;
    CheckFile<<"IMM,OS,  MALE = "<<ssObsZCsSrv_msxn(IMMATURE,OLD_SHELL,  MALE)<<endl;
    CheckFile<<"MAT,NS,FEMALE = "<<ssObsZCsSrv_msxn(  MATURE,NEW_SHELL,FEMALE)<<endl;
    CheckFile<<"MAT,NS,  MALE = "<<ssObsZCsSrv_msxn(  MATURE,NEW_SHELL,  MALE)<<endl;
    CheckFile<<"MAT,OS,FEMALE = "<<ssObsZCsSrv_msxn(  MATURE,OLD_SHELL,FEMALE)<<endl;
    CheckFile<<"MAT,OS,  MALE = "<<ssObsZCsSrv_msxn(  MATURE,OLD_SHELL,  MALE)<<endl;
 END_CALCS   
    
//     !! CheckFile<<"yrsObsSrvBio_n "<<yrsObsSrvBio_n<<endl;
//     !! CheckFile <<"nObsZCsSrv "<<nObsZCsSrv<<endl;
//     !! CheckFile <<"yrsObsZCsSrv_n "<<yrsObsZCsSrv_n<<endl;
//     !! CheckFile <<"ssObsZCsSrv_msxn "<<ssObsZCsSrv_msxn<<endl;
    
    // for length data
    // first index,1 immat, 2 mature,1 new shell, 2 old shell, then female 1 male 2
    5darray obsSrvNatZs_msxnz(1,nMSs,1,nSCs,1,nSXs,1,nObsZCsSrv,1,nZBs)  // immat,mat,new, old survey length data,female,male,year then bin
 LOCAL_CALCS
    for (int m=1;m<=nMSs;m++) {
        for (int s=1;s<=nSCs;s++) {
            for (int x=1;x<=nSXs;x++) {
                for (int y=1;y<=nObsZCsSrv;y++) obsSrvNatZs_msxnz(m,s,x,y) = ptrMDS->pTSD->nAtZ_xsmyz(x,s,m,y);
                CheckFile<<"obsSrvNatZs_msxnz("<<m<<cc<<s<<cc<<x<<") = "<<endl<<obsSrvNatZs_msxnz(m,s,x)<<endl;
            }
        }
    }
 END_CALCS
    
    //fishery data: retained catch size freqs (males x shell condition)    
    3darray obsRetZCsTCF_snz(NEW_SHELL,OLD_SHELL,1,nObsRetZCsTCF,1,nZBs)         // MALE retained length data by shell condition in directed fishery
 LOCAL_CALCS
    for (int s=1;s<=nSCs;s++) obsRetZCsTCF_snz(s) = ptrMDS->pTCFR->nAtZ_syz(s);
    CheckFile<<" obs fishery retained male new_shell length "<<endl;
    CheckFile<<obsRetZCsTCF_snz(NEW_SHELL)<<endl;
    CheckFile<<" obs fishery retained male old_shell length "<<endl;
    CheckFile<<obsRetZCsTCF_snz(OLD_SHELL)<<endl;  
//     CheckFile<<" obs fishery retained male all_shell length "<<endl;
//     CheckFile<<obsRetZCsTCF_snz(ALL_SHELL)<<endl;  
 END_CALCS
 
    //fishery data: directed fishery catch size freqs (females, (males x shell condition))    
    matrix obsZCsTCFF_nz(1,nObsZCsTCFF,1,nZBs)              // FEMALE size comps in directed fishery
    3darray obsZCsTCFM_snz(1,ALL_SHELL,1,nObsZCsTCFM,1,nZBs) // MALE size comps by shell condition in directed fishery
 LOCAL_CALCS
    obsZCsTCFF_nz = ptrMDS->pTCFD->nAtZ_xsyz(FEMALE,ALL_SHELL);
    CheckFile<<" obsZCsTCFF_nz discarded female all_shell length "<<endl;
    CheckFile<<obsZCsTCFF_nz<<endl;
    
    for (int s=1;s<=ALL_SHELL;s++) obsZCsTCFM_snz(s) = ptrMDS->pTCFD->nAtZ_xsyz(MALE,s);
    CheckFile<<" obsZCsTCFM_snz discarded male new_shell length "<<endl;
    CheckFile<<obsZCsTCFM_snz(NEW_SHELL)<<endl;
    CheckFile<<" obsZCsTCFM_snz discarded male old_shell length "<<endl;
    CheckFile<<obsZCsTCFM_snz(OLD_SHELL)<<endl;  
//     CheckFile<<" obs fish_discmd discarded male all_shell length "<<endl;
//     CheckFile<<obsZCsTCFM_snz(ALL_SHELL)<<endl;  
 END_CALCS
    
    //fishery data: snow crab fishery bycatch size freqs (females, (males x shell condition))    
    matrix obsZCsSCFF_nz(1,nObsZCsSCF,1,nZBs)               // FEMALE bycatch size comps in snow crab fishery
    3darray obsZCsSCFM_snz(1,ALL_SHELL,1,nObsZCsSCF,1,nZBs) //   MALE bycatch size comps in snow crab fishery by shell condition
 LOCAL_CALCS
    obsZCsSCFF_nz = ptrMDS->pSCF->nAtZ_xsyz(FEMALE,ALL_SHELL);
    CheckFile<<" obsZCsSCFF_nz discarded female all_shell length "<<endl;
    CheckFile<<obsZCsSCFF_nz<<endl;
    
    for (int s=1;s<=ALL_SHELL;s++) obsZCsSCFM_snz(s) = ptrMDS->pSCF->nAtZ_xsyz(MALE,s);
    CheckFile<<" obsZCsSCFM_snz discarded male new_shell length "<<endl;
    CheckFile<<obsZCsSCFM_snz(NEW_SHELL)<<endl;
    CheckFile<<" obsZCsSCFM_snz discarded male old_shell length "<<endl;
    CheckFile<<obsZCsSCFM_snz(OLD_SHELL)<<endl;  
//     CheckFile<<" obsZCsSCFM_snz discarded male all_shell length "<<endl;
//     CheckFile<<obsZCsSCFM_snz(ALL_SHELL)<<endl;  
 END_CALCS
    
    //fishery data: BBRKC fishery bycatch size freqs (females, (males x shell condition))    
    matrix obsZCsRKFF_nz(1,nObsZCsRKF,1,nZBs)               // FEMALE bycatch size comps in BBRKC fishery
    3darray obsZCsRKFM_snz(1,ALL_SHELL,1,nObsZCsRKF,1,nZBs) //   MALE bycatch size comps in BBRKC fishery by shell condition
 LOCAL_CALCS
    obsZCsRKFF_nz = ptrMDS->pRKF->nAtZ_xsyz(FEMALE,ALL_SHELL);
    CheckFile<<" obsZCsRKFF_nz discarded female all_shell length "<<endl;
    CheckFile<<obsZCsRKFF_nz<<endl;
    
    for (int s=1;s<=ALL_SHELL;s++) obsZCsRKFM_snz(s) = ptrMDS->pRKF->nAtZ_xsyz(MALE,s);
    CheckFile<<" obsZCsRKFM_snz discarded male new_shell length "<<endl;
    CheckFile<<obsZCsRKFM_snz(NEW_SHELL)<<endl;
    CheckFile<<" obsZCsRKFM_snz discarded male old_shell length "<<endl;
    CheckFile<<obsZCsRKFM_snz(OLD_SHELL)<<endl;  
//     CheckFile<<" obsZCsRKFM_snz discarded male all_shell length "<<endl;
//     CheckFile<<obsZCsRKFM_snz(ALL_SHELL)<<endl;  
 END_CALCS
    
    //fishery data: groundfish trawl fishery discard catch size freqs (females, males)    
    3darray obsZCsGTF_xnz(1,nSXs,1,nObsZCsGTF,1,nZBs)             // groundfish trawl discards by sex
 LOCAL_CALCS
    CheckFile<<"nObsZCsGTF = "<<tb<<nObsZCsGTF<<endl;
    for (int x=1;x<=nSXs;x++) obsZCsGTF_xnz(x) = ptrMDS->pGTF->nAtZ_xsyz(x,ALL_SHELL);
    CheckFile<<"obsZCsGTF_xnz(FEMALE)"<<endl;
    CheckFile<<obsZCsGTF_xnz(FEMALE)<<endl;
    CheckFile<<"obsZCsGTF_xnz(MALE)"<<endl;
    CheckFile<<obsZCsGTF_xnz(MALE)<<endl;
 END_CALCS
    
    int nYrsTCF;
    int nlog_sel50_dev_3                      //number of years of directed fishery post 1990
    ivector hasDirectedFishery(styr,endyr-1); //flags indicating if directed fishery is prosecuted (>0)
    vector obsRetCatchNum(styr,endyr-1)       //retained catch, numbers                 (IMPORTANT CHANGE: used to be "1965,endyr")
    vector obsRetCatchBio(styr,endyr-1)       //retained catch, millions of lbs of crab (IMPORTANT CHANGE: used to be "1965,endyr")
 LOCAL_CALCS
    nYrsTCF = 0;
    nlog_sel50_dev_3    = 0;
    hasDirectedFishery  = 0;//set all years to "no directed fishery"
    obsRetCatchNum.initialize();
    obsRetCatchBio.initialize();
    for (int i=1;i<=ptrMDS->pTCFR->nyCatch;i++) {
        int y = ptrMDS->pTCFR->yrsCatch(i);
        if ((styr<=y)&&(y<endyr)){
            nYrsTCF++;
            hasDirectedFishery(y) = 1;
            obsRetCatchNum(y)     = ptrMDS->pTCFR->catch_ty(1,i);
            obsRetCatchBio(y)     = ptrMDS->pTCFR->catch_ty(2,i);
        }
        if (y>1990) nlog_sel50_dev_3++;
    }
    CheckFile<<"number of directed fishery years after 1990: "<<nlog_sel50_dev_3<<endl;
    CheckFile<<"retained numbers: obsRetCatchNum"     <<endl<<obsRetCatchNum<<endl;
    CheckFile<<"retained biomass (mt): obsRetCatchBio"<<endl<<obsRetCatchBio<<endl;
    obsRetCatchBio /= 2.2045; // convert from millions lbs to 1000's tons   
 END_CALCS
 
    matrix obsDscBioTCF_xn(1,nSXs,1,nObsDscTCF)    // observed directed discard catch millions lbs female,male
    !! obsDscBioTCF_xn = ptrMDS->pTCFD->catch_xy;
    !! CheckFile <<"obsDscBioTCF_xn (millions lbs)"<<endl;
    !! CheckFile <<obsDscBioTCF_xn<<endl;
    !! obsDscBioTCF_xn /= 2.2045;                  // in 1000's of tons
    
    matrix obsDscBioSCF_xn(1,nSXs,1,nObsDscSCF)   // observed snow discard million lbs female,male
    !! cout<<wts::getBounds(obsDscBioSCF_xn)<<endl;
    !! cout<<wts::getBounds(ptrMDS->pSCF->catch_xy)<<endl;
    !! obsDscBioSCF_xn = ptrMDS->pSCF->catch_xy;
    !! CheckFile <<"obsDscBioSCF_xn (millions lbs)"<<endl;
    !! CheckFile <<obsDscBioSCF_xn<<endl;
    !! obsDscBioSCF_xn /= 2.2045;                 // in 1000's of tons
    
    matrix obsDscBioRKF_xn(1,nSXs,1,nObsDscRKF)     // observed red king discard millions lbs female,male
    !! obsDscBioRKF_xn = ptrMDS->pRKF->catch_xy;
    !! CheckFile <<"obsDscBioRKF_xn (millions lbs)"<<endl;
    !! CheckFile <<obsDscBioRKF_xn<<endl;
    !! obsDscBioRKF_xn /= 2.2045;                  // in 1000's of tons
        
    vector obsDscBioGTF_n(1,nObsDscGTF)             // trawl bycatch millions lbs sex combined need to apply mort 80%
    !! obsDscBioGTF_n = ptrMDS->pGTF->catch_y;
    !! CheckFile <<"obsDscBioGTF_n (millions lbs)"<<endl<<tb<<obsDscBioGTF_n<<endl;
    !! obsDscBioGTF_n /= 2.2045;                   // convert to 1000s tons
    
    // Survey indices
    vector obsSrvNum_n(1,nObsSrvBio)         // survey numbers (total) in millions of crab
    matrix obsSrvCV_xn(1,nSXs,1,nObsSrvBio)  // survey cv by sex x year
 LOCAL_CALCS
    obsSrvNum_n = ptrMDS->pTSD->abund_y;
    CheckFile<<"obsSrvNum_n = "<<endl<<tb<<obsSrvNum_n<<endl;
    obsSrvNum_n /= 1000.;                        // change to billions of crabs  (why?! later multiplied by 1000000)
    obsSrvCV_xn = ptrMDS->pTSD->cvsAbund_xy;
    CheckFile<<"obsSrvCV_xn = "<<endl<<obsSrvCV_xn<<endl;
 END_CALCS   
    
    matrix wtf_mz(1,nMSs,1,nZBs)  // weight at length juvenile and mature females (from kodiak program) in kg (??)
    vector wtm_z(1,nZBs)          // weight at length males (same as used in kodiak)                    in kg (??)
 LOCAL_CALCS
    wtf_mz(IMMATURE) = ptrMDS->pBio->wAtZ_xmz(FEMALE,IMMATURE);
    wtf_mz(  MATURE) = ptrMDS->pBio->wAtZ_xmz(FEMALE,  MATURE);
    wtm_z = ptrMDS->pBio->wAtZ_xmz(MALE,MATURE);
    CheckFile<<"wtf_mz = "<<endl<<wtf_mz<<endl;
    CheckFile<<"wtm_z  = "<<endl<<tb<<wtm_z<<endl;
 END_CALCS   
 
    vector obsPrMatureM_z(1,nZBs)          // logistic maturity probability curve for new shell immature males
    !!obsPrMatureM_z = ptrMDS->pBio->prMature_xz(MALE);//not given for females
    !!CheckFile<<"obsPrMatureM_z"<<endl<<tb<<obsPrMatureM_z<<endl;
    
    matrix obsAvgMatNS_xz(1,nSXs,1,nZBs)     // probability mature for new immature females, average proportion mature by length for new males
    !!obsAvgMatNS_xz(FEMALE) = ptrMDS->pBio->frMature_xsz(FEMALE,NEW_SHELL);
    !!obsAvgMatNS_xz(  MALE) = ptrMDS->pBio->frMature_xsz(  MALE,NEW_SHELL);
    !!CheckFile<<"obsAvgMatNS_xz(FEMALE)"<<endl<<tb<<obsAvgMatNS_xz(FEMALE)<<endl;
    !!CheckFile<<"obsAvgMatNS_xz(  MALE)"<<endl<<tb<<obsAvgMatNS_xz(  MALE)<<endl;
    
    matrix obsAvgMatOS_xz(1,nSXs,1,nZBs) // average proportion mature by length for old shell females, males
    !!obsAvgMatOS_xz(FEMALE) = ptrMDS->pBio->frMature_xsz(FEMALE,OLD_SHELL);
    !!obsAvgMatOS_xz(  MALE) = ptrMDS->pBio->frMature_xsz(  MALE,OLD_SHELL);
    !!CheckFile<<"obsAvgMatOS_xz(FEMALE)"<<endl<<tb<<obsAvgMatOS_xz(FEMALE)<<endl;
    !!CheckFile<<"obsAvgMatOS_xz(  MALE)"<<endl<<tb<<obsAvgMatOS_xz(  MALE)<<endl;
      
    vector length_bins(1,nZBs)             // Midpoints of length bins
    !!length_bins = ptrMDS->pBio->zBins-0.5;//IMPORTANT: subtract 0.5 to match up with old numbers
    !!CheckFile <<"length_bins"<<endl<<tb<<length_bins<<endl;
    
    vector mdptFshs_y(styr,endyr-1)        // Timing of catches  (IMPORTANT CHANGE: used to be "endyr")
 LOCAL_CALCS
    mdptFshs_y.initialize();
    for (int i=1;i<=ptrMDS->pBio->nyFshSeasons;i++){
        int y = (int) ptrMDS->pBio->mdptFshSeasons_yc(i,1);
        if ((styr<=y)&&(y<endyr)) mdptFshs_y(y) = ptrMDS->pBio->mdptFshSeasons_yc(i,2);
    }
    CheckFile<<"mdptFshs_y"<<endl<<tb<<mdptFshs_y<<endl;
 END_CALCS   
//     init_vector catch_ghl(1979,endyr)             // historical CPUE data                              //fixed index!
//     !! CheckFile<<"catch_ghl"<<endl;
//     !! CheckFile<<catch_ghl<<endl;
    
    vector effSCF_y(1978,endyr-1)       //fixed index          
 LOCAL_CALCS
    effSCF_y.initialize();
     for (int i=1;i<=ptrMDS->pSCF->nyEff;i++){
         int y = (int) ptrMDS->pSCF->effort_yc(i,1);
         if ((1978<=y)&&(y<endyr)) effSCF_y(y) = ptrMDS->pSCF->effort_yc(i,2);
    }
    CheckFile<<"effSCF_y"<<endl<<tb<<effSCF_y<<endl;
 END_CALCS   
 
    vector effRKF_y(1953,endyr-1)       //fixed index          wts: now includes rkceffortjap values
 LOCAL_CALCS
    effRKF_y.initialize();
     for (int i=1;i<=ptrMDS->pRKF->nyEff;i++){
         int y = (int) ptrMDS->pRKF->effort_yc(i,1);
         if ((1953<=y)&&(y<endyr)) effRKF_y(y) = ptrMDS->pRKF->effort_yc(i,2);
    }
    CheckFile<<"effRKF_y"<<endl<<tb<<effRKF_y<<endl;
 END_CALCS
    
    vector effTCF_y(1968,endyr-1)     //fixed index
 LOCAL_CALCS
    effTCF_y.initialize();
     for (int i=1;i<=ptrMDS->pTCFD->nyEff;i++){
         int y = (int) ptrMDS->pTCFD->effort_yc(i,1);
         if ((1968<=y)&&(y<endyr)) effTCF_y(y) = ptrMDS->pTCFD->effort_yc(i,2);
    }
    CheckFile<<"effTCF_y"<<endl<<tb<<effTCF_y<<endl;
 END_CALCS
    !!cout<<"Finished reading data files"<<endl;
    !!cout<<"#--------------------------------------------------"<<endl;
    !!CheckFile<<"Finished reading data files"<<endl;
    !!CheckFile<<"#--------------------------------------------------"<<endl;
    // End of reading normal data file 
        
    //---------------------------------------------------------------------------------------
    // Open control file....
    !! ad_comm::change_datafile_name(ptrMC->fnCtl);
    !!CheckFile<<"--------------------------------------"<<endl;
    !!CheckFile<<"Reading control file ''"<<ptrMC->fnCtl<<"'"<<endl;
    
    init_int inpVerMCF   //input version number for control file
 LOCAL_CALCS
    CheckFile<<"Model ControlFile version = "<<inpVerMCF<<endl;
    if (inpVerMCF!=verModelControlFile){
        cout<<"Model Control File version inconsistent with model."<<endl;
        cout<<"Current version = "<<verModelControlFile<<endl;
        cout<<"Version in file = "<<inpVerMCF<<endl;
        cout<<"Please use correct version."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
 END_CALCS
    init_number multQ                 // Q  mult by pop biomass to get survey biomass
    init_int phsM                     // phase to turn on ordinary M component estimation
    init_vector M_in(1,nSXs)          // base value for natural mortality by sex
    init_vector M_matn_in(1,nSXs)     // base value for natural mortality mature new shell by sex
    init_vector M_mato_in(1,nSXs)     // base value for natural mortality mature old shell by sex
    3darray baseM_msx(1,nMSs,1,nSCs,1,nSXs) //convenience array for all base M's
    init_number Mmult_imat_in         // initial value for natural mortality multiplier
    init_number Mmultm_in             // initial value for natural mortality multiplier
    init_number Mmultf_in             // initial value for natural mortality multiplier
    3darray inpMfac_msx(1,nMSs,1,nSCs,1,nSXs) //convenience array for all initial M multipliers
    init_vector mat_big_in(1,nSXs)    // natural mortality multiplier by sex during "kill 'em" period
    init_int phase_moltingp           // phase to estimate molting prob for mature males
    init_int phase_fishsel            // phase to estimate dome shape parameters for fishery selectivities
    init_int survsel_phase            // switch for which survey selectivty to use for 1989 to present - positive estimated negative fixed at somerton and otto
    init_int survsel1_phase           // switch for fixing all survey selectivities to somerton and otto - <0 fix, >0 estimate
    init_int phase_logistic_sel       // phase to estimate selectivities using logistic function
    init_vector sel_som(1,5)          // parameters for somerton-otto selectivity curve
 LOCAL_CALCS
    baseM_msx(IMMATURE,NEW_SHELL,FEMALE) = M_in(FEMALE);
    baseM_msx(IMMATURE,NEW_SHELL,  MALE) = M_in(  MALE);
    baseM_msx(IMMATURE,OLD_SHELL,FEMALE) = M_in(FEMALE);
    baseM_msx(IMMATURE,OLD_SHELL,  MALE) = M_in(  MALE);
    baseM_msx(  MATURE,NEW_SHELL,FEMALE) = M_matn_in(FEMALE);
    baseM_msx(  MATURE,NEW_SHELL,  MALE) = M_matn_in(  MALE);
    baseM_msx(  MATURE,OLD_SHELL,FEMALE) = M_mato_in(FEMALE);
    baseM_msx(  MATURE,OLD_SHELL,  MALE) = M_mato_in(  MALE);
    CheckFile<<"multQ              = "<<multQ<<endl;
    CheckFile<<"phsM               = "<<phsM<<endl;
    CheckFile<<"M_in               = "<<M_in<<endl;
    CheckFile<<"M_matn_in          = "<<M_matn_in<<endl;
    CheckFile<<"M_mato_in          = "<<M_mato_in<<endl;
    CheckFile<<"Mmult_imat_in      = "<<Mmult_imat_in<<endl;
    CheckFile<<"Mmultm_in          = "<<Mmultm_in<<endl;
    CheckFile<<"Mmultf_in          = "<<Mmultf_in<<endl;
    CheckFile<<"mat_big_in         = "<<mat_big_in<<endl;
    CheckFile<<"phase_moltingp     = "<<phase_moltingp<<endl;
    CheckFile<<"phase_fishsel      = "<<phase_fishsel<<endl;
    CheckFile<<"survsel_phase      = "<<survsel_phase<<endl;
    CheckFile<<"survsel1_phase     = "<<survsel1_phase<<endl;
    CheckFile<<"phase_logistic_sel = "<<phase_logistic_sel<<endl;
    CheckFile<<"sel_som            = "<<sel_som<<endl;
 END_CALCS
            
    init_vector wt_like(1,8)              // weights for selectivity likelihoods 1 fishery female, 2 survey female, 3 fishery male, 4 survey male
    init_vector like_wght(1,NUM_LEN_LIKE) // likelihood weights for size composition data
    init_number like_wght_CatchBio        // likelihood weight for fishery catch biomass fits [was like_wght_(6)]
    init_number like_wght_fbio            // likelihood weight for female biomass fit [was like_wght(5)]
    init_number like_wght_mbio            // likelihood weight for male biomass fit
    init_number like_wght_rec         // ??
    init_number like_wght_recf        // ??
    init_number like_wght_sexr        // ??
    init_number like_wght_sel50       // ??
    !!CheckFile<<"wt_like            = "<<wt_like<<endl;
    !!CheckFile<<"like_wght          = "<<like_wght<<endl;
    !!CheckFile<<"like_wght_CatchBio = "<<like_wght_CatchBio<<endl;
    !!CheckFile<<"like_wght_fbio     = "<<like_wght_fbio<<endl;
    !!CheckFile<<"like_wght_mbio     = "<<like_wght_mbio<<endl;
    !!CheckFile<<"like_wght_rec      = "<<like_wght_rec<<endl;
    !!CheckFile<<"like_wght_recf     = "<<like_wght_recf<<endl;
    !!CheckFile<<"like_wght_sexr     = "<<like_wght_sexr<<endl;
    !!CheckFile<<"like_wght_sel50    = "<<like_wght_sel50<<endl;
    
    init_number hm_pot                                                // fraction of pot discards that die (.5)
    init_number hm_trawl                                               // fraction of trawl discards that die(.8)
    !!CheckFile<<"hm_pot  = "<<hm_pot<<endl;
    !!CheckFile<<"hm_trawl = "<<hm_trawl<<endl;
    
    init_number mate_ratio                                            // mating ratio (USED)
    init_number fraction_new_error                                    // accounts for shell error
    init_number maturity_switch                                       // Set > 0 for logistic maturity instead of fractions by year (males)
    init_int nages                                                    // number of ages to track for mature old shell 
    init_number wght_total_catch                                      // weight for total catch biomass (WHY here!) (UNUSED except for output)
    init_number wght_female_potcatch                                  // weight for female pot bycatch              (UNUSED except for output)
    init_number wt_lmlike
    !!CheckFile<<"fraction_new_error   = "<<fraction_new_error<<endl;
    !!CheckFile<<"mate_ratio           = "<<mate_ratio<<endl;
    !!CheckFile<<"maturity_switch      = "<<maturity_switch<<endl;
    !!CheckFile<<"nages                = "<<nages<<endl;
    !!CheckFile<<"wght_total_catch     = "<<wght_total_catch<<endl;
    !!CheckFile<<"wght_female_potcatch = "<<wght_female_potcatch<<endl;
    !!CheckFile<<"wt_lmlike            = "<<wt_lmlike<<endl;
    
    init_number srv3_qPriorWgt        //likelihood multiplier for survey q prior (off=0)
    init_number srv3_qPriorMean       //prior mean
    init_number srv3_qPriorStD        //prior variance
    init_number srv3_qFemPriorWgt     //likelihood multiplier for survey female q prior (off=0)
    init_number srv3_qFemPriorMean    //prior mean
    init_number srv3_qFemPriorStD     //prior variance
    init_int mort_switch                                            // extra mort on 1, off 0
    init_int mort_switch2                                           // apply mort_switch correctly (1), or as in 2012 (0)
    init_int lyr_mort                                               // start yr extra mort
    init_int uyr_mort                                               //end yr extra mort
    !!CheckFile<<"srv3_qPriorWgt  = "<<srv3_qPriorWgt<<endl;
    !!CheckFile<<"srv3_qPriorMean = "<<srv3_qPriorMean<<endl;
    !!CheckFile<<"srv3_qPriorStD  = "<<srv3_qPriorStD<<endl;
    !!CheckFile<<"srv3_qFemPriorWgt  = "<<srv3_qFemPriorWgt<<endl;
    !!CheckFile<<"srv3_qFemPriorMean = "<<srv3_qFemPriorMean<<endl;
    !!CheckFile<<"srv3_qFemPriorStD  = "<<srv3_qFemPriorStD<<endl;
    !!CheckFile<<"mort_switch    = "<<mort_switch<<endl;
    !!CheckFile<<"mort_switch2   = "<<mort_switch2<<endl;
    !!CheckFile<<"lyr_mort       = "<<lyr_mort<<endl;
    !!CheckFile<<"uyr_mort       ="<<uyr_mort<<endl;
    
    init_int doPenRed           //flag (0/1) to reduce penalties on fishing-related devs by phase
    init_int phsFmRKF           //phase to turn on RKF fishing estimation (originally -4 or 5)
    init_number llw_sel50_dev_3 //llw for penalty on log_sel50_dev_3 (originally 0))
    init_number bnd_sel50_dev_3 //upper/lower bounds on log_sel50_dev_3 deviations (originally 0.5)
    !!CheckFile<<"doPenRed = "<<doPenRed<<endl;    
    !!CheckFile<<"phsFmRKF = "<<phsFmRKF<<endl;    
    !!CheckFile<<"llw_sel50_dev_3 = "<<llw_sel50_dev_3<<endl;    
    !!CheckFile<<"bnd_sel50_dev_3 = "<<bnd_sel50_dev_3<<endl;  
    
    //new 20150318
    init_int optFshNLLs  //flag indicating error model for fishery catch data (0=norm2, 1=lognormal)
    !!CheckFile<<"optFshNLLs = "<<optFshNLLs<<endl;    
    //the following obsErr are sd for normal error, cv for lognormal error
    init_number obsErrTCFR  //assumed observation error level for retained catch data
    init_number obsErrTCFD  //assumed observation error level for discard catch data in TCF
    init_number obsErrSCF   //assumed observation error level for discard catch data in SCF
    init_number obsErrRKF   //assumed observation error level for discard catch data in RKF
    init_number obsErrGTF   //assumed observation error level for discard catch data in GTF
    !!CheckFile<<"obsErrTCFR = "<<obsErrTCFR<<endl;   
    !!CheckFile<<"obsErrTCFD = "<<obsErrTCFD<<endl;   
    !!CheckFile<<"obsErrSCF  = "<<obsErrSCF<<endl;   
    !!CheckFile<<"obsErrRKF  = "<<obsErrRKF<<endl;   
    !!CheckFile<<"obsErrGTF  = "<<obsErrGTF<<endl;   
    
    //new 20150601
    init_int optFM;//original or gmacs-like fishing mortality
 LOCAL_CALCS
    if (optFM!=1) optFM = 0;//if gmacs not specified, then use original
    CheckFile<<"#---Fishing mortality model options"<<endl;
    CheckFile<<"optFM = "<<optFM<<endl;
 END_CALCS
    
    //new 20150823: revised 20151005
    init_int optFshSel;//option to force asymptote=1 for logistic fishery selectivity functions
    init_int optSrvSel;//option to force asymptote=1 for logistic survey selectivity functions
 LOCAL_CALCS
    if (optFshSel!=1) optFshSel = 0;//if not 1 (on), use old approach to normalization
    if (optSrvSel!=1) optSrvSel = 0;//if not 1 (on), use old approach to normalization
    CheckFile<<"#---Options for asymptotic behavior of selectivity functions"<<endl;
    CheckFile<<"optFshSel = "<<optFshSel<<endl;
    CheckFile<<"optSrvSel = "<<optSrvSel<<endl;
 END_CALCS
    
    //new 20160226: added female fishing mortality offsets
    init_int phsTCFF;//phase to turn on estimation of female fishing mortality in TCF
    init_int phsSCFF;//phase to turn on estimation of female fishing mortality in SCF
    init_int phsRKFF;//phase to turn on estimation of female fishing mortality in RKF
    init_int phsGTFF;//phase to turn on estimation of female fishing mortality in GTF
 LOCAL_CALCS
    CheckFile<<"#---Phases for female fishing mortality"<<endl;
    CheckFile<<"phsTCFF = "<<phsTCFF<<endl;
    CheckFile<<"phsSCFF = "<<phsSCFF<<endl;
    CheckFile<<"phsRKFF = "<<phsRKFF<<endl;
    CheckFile<<"phsGTFF = "<<phsGTFF<<endl;
 END_CALCS
            
    //new 20160316: options for calculating size comps
    init_int optPrNatZ_GTF; //
 LOCAL_CALCS
    CheckFile<<"#---Phases for female fishing mortality"<<endl;
    CheckFile<<"optPrNatZ_GTF = "<<optPrNatZ_GTF<<endl;
 END_CALCS
            
    //new 20160324: set info for estimating "historic" recruitment
    init_int    phsMnLnRecHist;  //phase to start estimating "historic" ln-scale mean recruitment (pMnLnRecHist)
    init_int    phsRecDevsHist;  //phase to start estimating "historic" rec devs (pRecDevsHist)
    init_number inpMnLnRecHist;  //initial value for "historic" ln-scale mean recruitment (pMnLnRecHist)
    init_int    mnYrRecDevsHist; //year to start estimating "historic" rec devs
 LOCAL_CALCS
    if (mnYrRecDevsHist<=0) mnYrRecDevsHist = mnYr;
    CheckFile<<"#---'Historic' recruitment info"<<endl;
    CheckFile<<"phsMnLnRecHist  = "<<phsMnLnRecHist <<tb<<"#---initial phase to estimate 'historic' mean ln-scale recruitment"<<endl;
    CheckFile<<"phsRecDevsHist  = "<<phsRecDevsHist <<tb<<"#---initial phase to estimate 'historic' recruitment deviations"<<endl;
    CheckFile<<"inpMnLnRecHist  = "<<inpMnLnRecHist <<tb<<"#---initial value for 'historic' mean ln-scale recruitment"<<endl;
    CheckFile<<"mnYrRecDevsHist = "<<mnYrRecDevsHist<<tb<<"#---Start year for 'historic' recruitment deviations"<<endl;
 END_CALCS
            
    //new 20160324: set info for estimating "current" recruitment
    init_int    phsMnLnRec;  //phase to start estimating "current" ln-scale mean recruitment (pMnLnRec)
    init_int    phsRecDevs;  //phase to start estimating "current" rec devs (pRecDevs)
    init_number inpMnLnRec;  //initial value for "current" ln-scale mean recruitment (pMnLnRec)
    init_int    mnYrRecCurr; //year to start estimating "current" recruitment
 LOCAL_CALCS
    CheckFile<<"#---'Current' recruitment info"<<endl;
    CheckFile<<"phsMnLnRec  = "<<phsMnLnRec <<tb<<"#---initial phase to estimate 'current' mean ln-scale recruitment"<<endl;
    CheckFile<<"phsRecDevs  = "<<phsRecDevs <<tb<<"#---initial phase to estimate 'current' recruitment deviations"<<endl;
    CheckFile<<"inpMnLnRec  = "<<inpMnLnRec <<tb<<"#---initial value for 'current' mean ln-scale recruitment"<<endl;
    CheckFile<<"mnYrRecCurr = "<<mnYrRecCurr<<tb<<"#---Start year for 'current' recruitment"<<endl;
 END_CALCS
            
    //new 20160323: options for fitting male mortality in TCF
    init_int optTCFMfit;
 LOCAL_CALCS
    CheckFile<<"#---Option for fitting male mortality in TCF"<<endl;
    CheckFile<<"optTCFMfit = "<<optTCFMfit<<endl;
 END_CALCS
 
    //Finished reading control file
    !!CheckFile<<"Finished reading control file."<<endl;
    !!CheckFile<<"--------------------------------------"<<endl;
    // the rest are working variables 
    
    matrix obsDscBioMortTCF(1,nSXs,1,nObsDscTCF)   // observed discard MORTALITY in TCF (1000's tons)
    matrix obsDscBioMortSCF(1,nSXs,1,nObsDscSCF)   // observed discard MORTALITY in SCF (1000's tons)
    matrix obsDscBioMortRKF(1,nSXs,1,nObsDscRKF)   // observed discard MORTALITY in RKF (1000's tons)
    vector obsDscBioMortGTF(1,nObsDscGTF)    // observed discard MORTALITY in GTF (1000's tons)
    vector obsTotBioMortTCFM(1,nObsDscTCF)   // observed total male MORTALITY in TCF (1000's tons)
 LOCAL_CALCS
    obsSrvNum_n=obsSrvNum_n*1000000;       // survey numbers read in are millions of crab (but /1000 above so obsSrvNum_n in thousands of crab, now)
    wtf_mz=wtf_mz*0.001;                     // change weights from kg to tons
    wtm_z =wtm_z*0.001;                     // change weights from kg  to tons
    //apply discard mortality to observed discards
    obsDscBioMortTCF = hm_pot*obsDscBioTCF_xn;//apply discard mortality to catches
    obsDscBioMortSCF = hm_pot*obsDscBioSCF_xn;
    obsDscBioMortRKF = hm_pot*obsDscBioRKF_xn;
    obsDscBioMortGTF = hm_trawl*obsDscBioGTF_n;
    // Calculate TOTAL male mortality in Tanner crab fishery
    obsTotBioMortTCFM.initialize();
    obsTotBioMortTCFM = obsDscBioMortTCF(MALE)+obsRetCatchBio(yrsObsDscTCF_n);
//    for (int i=1;i<=nObsDscTCF;i++) {
//        obsTotBioMortTCFM(i) = obsDscBioMortTCF(MALE,i)+obsRetCatchBio(yrsObsDscTCF_n(i));
//    }
    CheckFile<<"yrsObsDscTCF_n"<<endl<<yrsObsDscTCF_n<<endl;
    CheckFile<<"obsDscBioMortTCF(  MALE)"<<endl<<obsDscBioMortTCF(  MALE)<<endl;
    CheckFile<<"obsDscBioMortTCF(FEMALE)"<<endl<<obsDscBioMortTCF(FEMALE)<<endl;
    CheckFile<<"obsRetBioMortTCFM"<<endl<<obsRetCatchBio(yrsObsDscTCF_n)<<endl;
    CheckFile<<"obsTotBioMortTCFM"<<endl<<obsTotBioMortTCFM<<endl;
    CheckFile<<"yrsObsDscSCF"<<endl<<yrsObsDscSCF<<endl;
    CheckFile<<"obsDscBioMortSCF(  MALE)"<<endl<<obsDscBioMortSCF(  MALE)<<endl;
    CheckFile<<"obsDscBioMortSCF(FEMALE)"<<endl<<obsDscBioMortSCF(FEMALE)<<endl;
    CheckFile<<"yrsObsDscRKF"<<endl<<yrsObsDscRKF<<endl;
    CheckFile<<"obsDscBioMortRKF(  MALE)"<<endl<<obsDscBioMortRKF(  MALE)<<endl;
    CheckFile<<"obsDscBioMortRKF(FEMALE)"<<endl<<obsDscBioMortRKF(FEMALE)<<endl;
    CheckFile<<"yrsObsDscGTF"<<endl<<yrsObsDscGTF<<endl;
    CheckFile<<"obsDscBioMortGTF"<<endl<<obsDscBioMortGTF<<endl;
    
 END_CALCS
    
    vector obsTotSrvNum(1,nObsZCsSrv)                       // Total survey numbers
    
    3darray obsPrNatZ_TCFR_sn(1,nSCs,1,nObsRetZCsTCF,1,nZBs)   // length-frequency of retained catch?  
    3darray obsPrNatZ_TCFM_snz(1,nSCs,1,nObsRetZCsTCF,1,nZBs)   // length-frequency of total male catch in directed fishery
    3darray obs_p_fish_discm(1,nSCs,1,nObsZCsTCFM,1,nZBs)// males discards                             what's 1,2?
    matrix obsPrNatZ_TCFF_nz(1,nObsZCsTCFF,1,nZBs)        // female discards
    3darray obsPrNatZ_GTF_xnz(1,nSXs,1,nObsZCsGTF,1,nZBs)         // bycatch in trawl fishery
    
    5darray obsPrNatZ_Srv_msxnz(1,nMSs,1,nSCs,1,nSXs,1,nObsZCsSrv,1,nZBs)// Survey length frequency by maturity state, shell condition, sex
    matrix obsSrvMatBio_xy(1,nSXs,styr,endyr)                            // Survey mature biomass (by sex)
    3darray obsSrvMatNum_sxy(1,nSCs,1,nSXs,styr,endyr)                   // Survey mature numbers (by shell condition, sex)
    
    3darray obsPrNatZ_SCF_xnz(1,nSXs,1,nObsZCsSCF,1,nZBs)  
    3darray obsPrNatZ_RKF_xnz(1,nSXs,1,nObsZCsRKF,1,nZBs)    
    
    vector obsSrvLegalNum_n(1,nObsZCsSrv)                // Legal male numbers, in survey
    vector obsSrvLegalBio_n(1,nObsZCsSrv)                // Legal male biomass, in survey
    
    3darray obsSrvNum_xyz(1,nSXs,styr,endyr,1,nZBs)  // Survey numbers, by sex
    vector obsSrvNum_y(styr,endyr)                   // Total survey numbers
    matrix obsSrvBio_xy(1,nSXs,styr,endyr)           // Survey biomass, by sex
    vector obsSrvBio_y(styr,endyr)                   // Total survey biomass
    
    number avgwt2
    number avgwtall
    vector avgwt(styr,endyr)
    
    int phs_mat_big;
    !!if (mort_switch==1) phs_mat_big = phsM+1; else phs_mat_big = -8;
    
 LOCAL_CALCS
    dmX   = "x=c("+qt+STR_FEMALE+qt    +cc+ qt+STR_MALE+qt+")";
    dmS   = "s=c("+qt+STR_NEW_SHELL+qt +cc+ qt+STR_OLD_SHELL+qt+")";
    dmM   = "m=c("+qt+STR_IMMATURE+qt  +cc+ qt+STR_MATURE+qt+")";
    dmY   = "y="+str(styr)+":"+str(endyr);
    dmYm1 = "y="+str(styr)+":"+str(endyr-1);  
    dmZ   = "z=c("+wts::to_qcsv(length_bins)+")";
 END_CALCS
            
    !!CheckFile<<"End of DATA_SECTION-------------------------"<<endl;
    !!CheckFile<<"--------------------------------------------"<<endl<<endl;

// =======================================================================
// =======================================================================
INITIALIZATION_SECTION
    pMnLnRec 11.4
    pAvgLnF_TCF -0.7
    //  log_avg_fmortdf -1.0
    pAvgLnF_GTF  -4.0
    pAvgLnF_SCF -3.0
    pAvgLnF_RKF -5.25       //to initialize same as TCSAM_WTS for 2013-09
    log_avg_sel50_mn  4.87  //this is 130.3 mm
//    pF_DevsTCF 0.00001                           //wts: dev.s should be mean 0!
    matestf -1.0
    matestm -1.0
    mat_big  1.0      //<-NEW by wts!!
    //  selGTFF_slpA 0.05
    //  selGTFF_z50A 85.0
    //  selGTFM_slpA 0.07
    //  selGTFM_z50A 65.0
    
    //  moltp_af 2.0
    //  moltp_bf 150.
    //  moltp_am 0.02
    //  moltp_bm 300.
    //  moltp_ammat 0.05
    //  moltp_bmmat 105.0
    
    //  af 15.75
    //  bf 1.01
    //  am2 15.75
    //  bm2 1.07
    
    //  srv1_slope 0.07
    //  srv1_sel50 60.0
    //  srv1_q 1.0
    //  srv2_q 1.0
    //  srv3_q 1.0
    //  srv1_sel95 100
    //  srv1_sel50 60 
    //  srv2_sel95 100
    //  srv2_sel50  60
    //  srv3_sel95 100
    //  srv3_sel50  60
    //  fish_fit_sel50_mn 95.1
    //  log_sel50_dev_mn  0.0000

    selSCFM_z50A1 80.0
    selSCFM_z50A2 80.0
    selSCFM_z50A3 80.0
    selSCFM_lnZ50D1 4.0
    selSCFM_lnZ50D2 4.0
    selSCFM_lnZ50D3 4.0

    pAvgLnF_TCFF 0.0
    pAvgLnF_SCFF 0.0
    pAvgLnF_RKFF 0.0
    pAvgLnF_GTFF 0.0
 
// =======================================================================
// =======================================================================
PARAMETER_SECTION
    
    init_bounded_number af1(0.4,0.7,8)                       // Female growth-increment
    init_bounded_number bf1(0.6,1.2,8)                       // Female growth-increment
    init_bounded_number am1(0.3,0.6,8)                       // Male growth-increment
    init_bounded_number bm1(0.7,1.2,8)                       // Male growth-increment
    
    init_bounded_vector growth_beta(1,nSXs,0.75000,0.75001,-2)  // Growth beta                                //this is NOT estimated (why?)
    init_bounded_number Mmult_imat(0.2,2.0,phsM)                   // natural mortality multiplier for females and males
    init_bounded_number Mmultm(0.1,1.9,phsM)                       // natural mortality multiplier for mature new and old shell male
    init_bounded_number Mmultf(0.1,1.9,phsM)                       // natural mortality multiplier for mature new and old shell female
    init_bounded_vector mat_big(1,nSXs,0.1,10.0,phs_mat_big)    // mult. on 1980-1980 M for mature males and females                     
    init_bounded_number alpha1_rec(11.49,11.51,-8)              // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    init_bounded_number beta_rec(3.99,4.01,-8)                  // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    
    init_bounded_number moltp_af(0.04,3.0,-6)                    // paramters for logistic function molting   //this is NOT estimated (why?)
    init_bounded_number moltp_bf(130.,300.,-6)                   // female                                    //this is NOT estimated (why?)
    init_bounded_number moltp_am(0.04,3.0,-5)                    // paramters for logistic function molting   //this is NOT estimated (why?)
    init_bounded_number moltp_bm(130.0,300.0,-5)                 // immature males                            //this is NOT estimated (why?)
    init_bounded_number moltp_ammat(.0025,3.0,phase_moltingp)    // logistic molting prob for mature males
    init_bounded_number moltp_bmmat(1,120,phase_moltingp)        // logistic molting prob for mature males
    
    init_number pMnLnRec(phsMnLnRec)                                                         // Mean ln-scale "current" recruitment 1974+ (males, females are equal)
    init_bounded_dev_vector pRecDevs(mnYrRecCurr,endyr,-15,15,phsRecDevs)                    // "current" recruitment devs
    init_number pMnLnRecHist(phsMnLnRecHist)                                                 // Mean ln-scale "historic" recruitment
    init_bounded_dev_vector pRecDevsHist(mnYrRecDevsHist,mnYrRecCurr-1,-15,15,phsRecDevsHist)// "historic" ln-scale recruitment deviations
    vector rec_y(styr,endyr)                                                                 //arithmetic-scale recruitments (1000's ??)
    
    //20150601: changed ...Fm... to ...F_... because of ambiguity as to whether
    //they represent fishing mortality (original model) of capture (gmacs) rates
    init_number pAvgLnF_TCF(1)                                        //log-scale mean directed fishing mortality
    init_bounded_dev_vector pF_DevsTCF(1,nYrsTCF,-15,15,2)//log-scale directed fishing mortality devs IMPORTANT CHANGE: USED TO BE "1966,endyr-12"
    init_number pAvgLnF_GTF(2)                                        // fishing mortality (trawl)
    init_bounded_dev_vector pF_DevsGTF(1973,endyr-1,-15,15,3)         // trawl fishery f-devs       (IMPORTANT CHANGE: used to be "endyr") 1973 seems OK
    init_number pAvgLnF_SCF(3)                                        // fishing mortality snow crab fishery discards
    init_bounded_dev_vector pF_DevsSCF(1992,endyr-1,-15,15,4)         // snow crab fishery f-devs   (IMPORTANT CHANGE: used to be "endyr")  1992 is OK
    //wts:following lines were used in 2013 assessment (i.e., RKF Fm NOT estimated)
//    init_bounded_number pAvgLnF_RKF(-5.25,-5.25,-4)                   // fishing mortality red king crab fishery discards //this is NOT estimated (why?)
//    init_bounded_dev_vector pF_DevsRKF(1,nObsDscRKF,-15,15,-5) //this is NOT estimated (why?)  IMPORTANT CHANGEA: was nObsDscRKF-1.  why -1 in "nObsDscRKF-1"
//    //wts:following lines turn ON Fm for RKF
//    init_bounded_number pAvgLnF_RKF(-10,5,5)                   // fishing mortality red king crab fishery discards //IMPORTANT CHANGE: now estimated
//    init_bounded_dev_vector pF_DevsRKF(1,nObsDscRKF,-15,15,6) //this is NOT estimated (why?)  IMPORTANT CHANGE: now estimated!
    //wts:implemented 20140823 (estimation phase now set in control file)
    init_bounded_number pAvgLnF_RKF(-10,5,phsFmRKF)                           // fishing mortality red king crab fishery discards
    init_bounded_dev_vector pF_DevsRKF(1,nObsDscRKF,-15,15,phsFmRKF+1) // 
    
    // Selectivity pattern for males (directed fishery)
    // Set -phase so not estimated if using @3 selectivity periods
    init_bounded_number fish_slope_mn(0.1,0.4,-phase_logistic_sel)                //this is NOT estimated (why?)
    init_bounded_number log_avg_sel50_mn(4,5.0,-phase_logistic_sel)               //this is NOT estimated (why?)
    init_bounded_dev_vector log_sel50_dev_mn(1,nObsRetZCsTCF,-5,5,-phase_logistic_sel)//this is NOT estimated (why?)
    vector fish_sel50_mn(styr,endyr-1)  //(IMPORTANT CHANGE: used to be "endyr")
    
    // Retention function
    // init_bounded_number fish_fit_slope_mn(.250,1.001,phase_logistic_sel)
    // init_bounded_number fish_fit_sel50_mn(85.0,160.0,phase_logistic_sel)
    // 1981 - 1992
    init_bounded_number fish_fit_slope_mn1(00.250,001.001,phase_logistic_sel)
    init_bounded_number fish_fit_sel50_mn1(85.000,160.000,phase_logistic_sel)
    // 2005-endyr  
    init_bounded_number fish_fit_slope_mn2(00.250,002.001,phase_logistic_sel)
    init_bounded_number fish_fit_sel50_mn2(85.000,160.000,phase_logistic_sel)
    
    // Directed fishery selectivity pattern for period-1: 1993-1996
    init_bounded_number fish_slope_1(00.05,000.75,phase_logistic_sel)      
    init_bounded_number fish_sel50_1(50.00,170.00,-phase_logistic_sel)//this is NOT estimated! (why?)
    
    // Directed fishery selectivity pattern changing by year for period-3: 2005-P
    init_bounded_number fish_slope_yr_3(0.1,0.4,phase_logistic_sel)      
    init_bounded_number log_avg_sel50_3(4.0,5.0,phase_logistic_sel)
    init_bounded_dev_vector log_sel50_dev_3(1,nlog_sel50_dev_3,-bnd_sel50_dev_3,bnd_sel50_dev_3,phase_logistic_sel) //Fixed index (why 2000?) (IMPORTANT CHANGE: used to be "endyr-2000")
    
    // for a dome-shaped selex pattern
    init_bounded_number fish_slope_mn2(000.01,002.0,phase_fishsel)
    init_bounded_number fish_sel50_mn2(100.00,160.0,phase_fishsel)
    
    // Female discards
    init_bounded_number selTCFF_slp(00.1,000.4,phase_logistic_sel)
    init_bounded_number selTCFF_z50(80.0,150.0,phase_logistic_sel)
    
    // snow fishery female discards for period-1: 1989-1996
    init_bounded_number selSCFF_slpA1(00.05,000.5,phase_logistic_sel+1) //add 1 to phase
    init_bounded_number selSCFF_z50A1(50.00,150.0,phase_logistic_sel+1)
    
    // snow fishery female discards for period-2: 1997-2004
    init_bounded_number selSCFF_slpA2(00.05,000.5,phase_logistic_sel+1) //add 1 to phase
    init_bounded_number selSCFF_z50A2(50.00,120.0,phase_logistic_sel+1)
    
    // snow fishery female discards for period-3: 2005-P
    init_bounded_number selSCFF_slpA3(00.05,000.5,phase_logistic_sel+1) //add 1 to phase
    init_bounded_number selSCFF_z50A3(50.00,120.0,phase_logistic_sel+1)
    
    // snow fishery male discards for period-1: 1989-1996
    init_bounded_number selSCFM_slpA1(00.1,000.5,phase_logistic_sel+1)  //add 1 to phase
    init_bounded_number selSCFM_z50A1(40.0,140.0,phase_logistic_sel+1)
    init_bounded_number selSCFM_slpD1(00.1,000.5,phase_logistic_sel+1)
    init_bounded_number selSCFM_lnZ50D1(2,4.5,phase_logistic_sel+1)
    
    // snow fishery male discards for period-2: 1997-2004
    init_bounded_number selSCFM_slpA2(00.1,000.5,phase_logistic_sel+1)  //add 1 to phase
    init_bounded_number selSCFM_z50A2(40.0,140.0,phase_logistic_sel+1)
    init_bounded_number selSCFM_slpD2(00.1,000.5,phase_logistic_sel+1)
    init_bounded_number selSCFM_lnZ50D2(2,4.5,phase_logistic_sel+1)
    
    // snow fishery male discards for period-3: 2005-P
    init_bounded_number selSCFM_slpA3(00.1,000.5,phase_logistic_sel+1)  //add 1 to phase
    init_bounded_number selSCFM_z50A3(40.0,140.0,phase_logistic_sel+1)
    init_bounded_number selSCFM_slpD3(00.1,000.5,phase_logistic_sel+1)
    init_bounded_number selSCFM_lnZ50D3(2,4.5,phase_logistic_sel+1)  //was selSCF_Z50_ma2_1
    
    // red king fishery female discards
    init_bounded_number selRKFF_slpA1(00.05,000.5,phase_logistic_sel)     //add 2 to phase
    init_bounded_number selRKFF_z50A1(50.00,150.0,phase_logistic_sel) //add 2 to phase
    init_bounded_number selRKFF_slpA2(00.05,000.5,phase_logistic_sel) //add 2 to phase
    init_bounded_number selRKFF_z50A2(50.00,150.0,phase_logistic_sel) //add 2 to phase
    init_bounded_number selRKFF_slpA3(00.05,000.5,phase_logistic_sel) //add 2 to phase
    init_bounded_number selRKFF_z50A3(50.00,170.0,phase_logistic_sel) //add 2 to phase
    
    // red king fishery male discards
    init_bounded_number selRKFM_slpA1(.01,.50,phase_logistic_sel)          //add 2 to phase
    init_bounded_number selRKFM_z50A1(95.0,150.0,phase_logistic_sel)
    init_bounded_number selRKFM_slpA2(.01,.50,phase_logistic_sel)          //add 2 to phase
    init_bounded_number selRKFM_z50A2(95.0,150.0,phase_logistic_sel)
    init_bounded_number selRKFM_slpA3(.01,.50,phase_logistic_sel)          //add 2 to phase
    init_bounded_number selRKFM_z50A3(95.0,150.0,phase_logistic_sel)
    
    // Trawl fishery selectivity female, 1973-1987
    init_bounded_number selGTFF_slpA1(0.01,0.5,phase_logistic_sel)
    init_bounded_number selGTFF_z50A1(40.0,125.01,phase_logistic_sel)
    // Trawl fishery selectivity female, 1988-1996
    init_bounded_number selGTFF_slpA2(0.005,0.5,phase_logistic_sel)
    init_bounded_number selGTFF_z50A2(40.0,250.01,phase_logistic_sel) 
    // Trawl fishery selectivity female, 1997-P
    init_bounded_number selGTFF_slpA3(0.01,0.5,phase_logistic_sel)
    init_bounded_number selGTFF_z50A3(40.0,150.01,phase_logistic_sel)
    // Trawl fishery selectivity male, 1973-1987
    init_bounded_number selGTFM_slpA1(0.01,0.5,phase_logistic_sel)
    init_bounded_number selGTFM_z50A1(40.0,120.01,phase_logistic_sel)
    // Trawl fishery selectivity male, 1988-1996
    init_bounded_number selGTFM_slpA2(0.01,0.5,phase_logistic_sel)
    init_bounded_number selGTFM_z50A2(40.0,120.01,phase_logistic_sel)
    // Trawl fishery selectivity male, 1997-P
    init_bounded_number selGTFM_slpA3(0.01,0.5,phase_logistic_sel)
    init_bounded_number selGTFM_z50A3(40.0,120.01,phase_logistic_sel)
    // Tanner 1968 to 2008 use these to estimate survey selectivities - same for males and females
    // put negative for phase for q's, change to positive when use som and otto because input a negative phase 
    //  init_bounded_number srv1_slope(.01,.4,-1)
    //1969 to 1973  
    //  init_bounded_number srv1_q(0.20,1.001,survsel1_phase)
    //  init_bounded_number srv1_sel95(140,200.01,survsel1_phase)
    //    init_bounded_number srv1_sel50(108.334,108.33401,-survsel1_phase)
    //    init_bounded_number srv1_sel50(20.0,140.01,survsel1_phase)
    //1974 to 1981 
    init_bounded_number srv2_q(0.50,1.001,survsel1_phase)
    //  init_bounded_number srv2_sel95(120.409,120.4091,-survsel1_phase)
    //  init_bounded_number srv2_sel50(62.519,62.5191,-survsel1_phase)
    init_bounded_number srv2_seldiff(0.0,100.0,survsel1_phase)
    init_bounded_number srv2_sel50(0.0,90.0,survsel1_phase)
    //1982-86 net change; 1982 first year of 83-112; burn-in period
    init_bounded_number srv2a_q(0.1,1.00001,-survsel1_phase)
    init_bounded_number srv2a_seldiff(0.0,100.0,-survsel1_phase)
    init_bounded_number srv2a_sel50(10.0,120.01,-survsel_phase)
    //1987-P
    //max of the underbag at 182.5 mm is 0.873
    init_bounded_number srv3_q(0.2,2.0,survsel1_phase)
    init_bounded_number srv3_seldiff(0.0,100.0,survsel1_phase)
    init_bounded_number srv3_sel50(0.0,69.0,survsel_phase)    
    
    init_bounded_vector matestf(1,16,-15.0,0.0,5)
    init_bounded_vector matestm(1,nZBs,-15.0,0.0,5)
    
    init_bounded_number srv2_qFem(0.5,1.001,survsel1_phase)
    init_bounded_number srv2_seldiff_f(0.0,100.0,survsel1_phase)
    
    init_bounded_number srv2_sel50_f(-200.0,100.01,survsel1_phase)
    
    init_bounded_number srv2a_qFem(0.25,1.001,-survsel1_phase)
    init_bounded_number srv2a_seldiff_f(0.0,100.0,-survsel1_phase)
    init_bounded_number srv2a_sel50_f(-200.0,100.01,-survsel1_phase)
    
    init_bounded_number srv3_qFem(0.2,1.0,survsel1_phase)
    init_bounded_number srv3_seldiff_f(0.0,100.0,survsel1_phase)
    init_bounded_number srv3_sel50_f(-50.0,69.0,survsel1_phase)
    
    init_bounded_number prNewShellRecruits(1.0,1.0,-2)       // proportion new shell in recruits  NOT ESTIMATED

    init_bounded_number pAvgLnF_TCFF(-5.0,5.0,phsTCFF)  // female offset to ln-scale mean fishing mortality in directed fishery
    init_bounded_number pAvgLnF_SCFF(-5.0,5.0,phsSCFF)  // female offset to ln-scale mean fishing mortality in snow crab fishery
    init_bounded_number pAvgLnF_RKFF(-5.0,5.0,phsRKFF)  // female offset to ln-scale mean fishing mortality in BBRKC fishery
    init_bounded_number pAvgLnF_GTFF(-5.0,5.0,phsGTFF)  // female offset to ln-scale mean fishing mortality in groundfish trawl fisheries

    ////end of estimated parameters///////////////
    
    3darray retFcn(1,nSCs,styr,endyr-1,1,nZBs)    // Retention curve for males caught in directed fishery    (IMPORTANT CHANGE: used to be "endyr")
    3darray selTCFR(1,nSCs,styr,endyr-1,1,nZBs)   // full selectivity for retained males in directed fishery (IMPORTANT CHANGE: used to be "endyr")
    3darray selTCFM(1,nSCs,styr,endyr-1,1,nZBs)   // selectivity for all males in directed fishery           (IMPORTANT CHANGE: used to be "endy1")
    vector  selTCFF(1,nZBs)                       // selectivity for females in directed fishery             
    3darray selGTF(1,3,1,nSXs,1,nZBs)      // 3D array to accommodate 3 selectivity periods
    3darray selSCF(1,3,1,nSXs,1,nZBs)      // 3D array to accommodate 3 selectivity periods 
    3darray selRKF(1,3,1,nSXs,1,nZBs)      // 3D array to accommodate 3 selectivity periods 
    
    matrix selSrv2(1,nSXs,1,nZBs)  // Survey selectivity 2
    matrix selSrv2a(1,nSXs,1,nZBs) // Survey selectivity 2a
    matrix selSrv3(1,nSXs,1,nZBs)  // Survey selectivity 3
    
    vector modPopNum_y(styr,endyr)   // Total population numbers on July 1, endyr (output)
    number popn_snowm         // Total population numbers (output)
    number popn_snowf         // Total population numbers (output)
    number popn_rkm           // Total population numbers (output)
    number popn_rkf           // Total population numbers (output)
    
    3darray M_msx(1,nMSs,1,nSCs,1,nSXs);//natural mortality rates
    
    vector pred_bio(styr,endyr) // Predicted biomass (determines sdrDepletion) 
    
    vector fspbio(styr,endyr)                       // Predicted female spawning biomass on July 1
    vector mspbio(styr,endyr)                       // Predicted   male spawning biomass on July 1  
//    vector fspbio_srv1(styr,endyr)                  // Predicted female spawning biomass at survey time, as seen by survey
//    vector mspbio_srv1(styr,endyr)                  // Predicted   male spawning biomass at survey time, as seen by survey
    matrix modSrvMatBio_xy(1,nSXs,styr,endyr)          // Predicted mature biomass at survey time, as seen by survey
    3darray modSrvMatNum_xsy(1,nSXs,1,nSCs,styr,endyr) // Predicted survey numbers for mature crab by sex, shell condition (output)
    
    vector numLegalMales_y(styr,endyr)              // Number of legal males at time of survey 
    vector bioLegalMales_y(styr,endyr)              // Biomass of legal males in survey (output)          
    
    matrix modSrvNum_xy(1,nSXs,styr,endyr)             // Survey abundance, by sex          
    matrix modSrvBio_xy(1,nSXs,styr,endyr)             // Survey biomass, by sex      
    vector modSrvNumLegal_y(styr,endyr)                // Survey-selected legal-sized males   
    matrix modSrvNumLegal_sy(1,nSCs,styr,endyr)        // Survey-selected legal-sized males, by shell condition  
    vector modSrvBioLegal_y(styr,endyr)                // Survey-selected legal male biomass  
//    matrix biom_tmp(1,nSXs,styr,endyr);              // Predicted survey indices        
    3darray modSrvNum_xyz(1,nSXs,styr,endyr,1,nZBs)    // Predicted survey abundance by length bin       
    4darray modSrvPrNatZ_NS_mxyz(1,nMSs,1,nSXs,styr,endyr,1,nZBs) // Predicted new shell length-frequency 
    4darray modSrvPrNatZ_OS_mxyz(1,nMSs,1,nSXs,styr,endyr,1,nZBs) // Predicted old shell length-frequency
    
    3darray modPrNatZ_TCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)  // Predicted proportion (total catch)                     (IMPORTANT CHANGE: used to be "endyr")
    3darray modPrNatZ_TCFR_syz(1,nSCs,styr-1,endyr,1,nZBs)  // Predicted retained catch proportions                   (IMPORTANT CHANGE: used to be "endyr")
    matrix  modPrNatZ_TCFF_yz(styr,endyr-1,1,nZBs)          // Predicted female discard proprtions in directed fishery(IMPORTANT CHANGE: used to be "endyr")
    3darray modPrNatZ_GTF_xyz(1,nSXs,styr,endyr-1,1,nZBs)   // Predicted trawl proportions                            (IMPORTANT CHANGE: used to be "endyr")
    3darray modPrNatZ_SCF_xyz(1,nSXs,styr,endyr-1,1,nZBs)   // Predicted snow crab fishery  proportions               (IMPORTANT CHANGE: used to be "endyr")
    3darray modPrNatZ_RKF_xyz(1,nSXs,styr,endyr-1,1,nZBs)   // Predicted red king crab proportions                    (IMPORTANT CHANGE: used to be "endyr")
    
    3darray natlength(1,nSXs,styr,endyr,1,nZBs)               // Total numbers by sex, length, and year
    3darray natlength_iold(1,nSXs,styr,endyr,1,nZBs)          // Immature old-shell numbers by sex, length, and year
    3darray natlength_inew(1,nSXs,styr,endyr,1,nZBs)          // Immature new-shell numbers by sex, length, and year
    3darray natlength_mold(1,nSXs,styr,endyr,1,nZBs)          // Mature old-shell numbers by sex, length, and year
    3darray natlength_mnew(1,nSXs,styr,endyr,1,nZBs)          // Mature new-shell numbers by sex, length, and year
    4darray natlength_mold_age(1,nSXs,styr,endyr,1,nages,1,nZBs)  // Age- and length-structure          
    3darray natlength_old(1,nSXs,styr,endyr,1,nZBs)           // Old-shell numbers by sex, length, and year
    3darray natlength_new(1,nSXs,styr,endyr,1,nZBs)           // New-shell numbers by sex, length, and year
    3darray natlength_i(1,nSXs,styr,endyr,1,nZBs)             // Immature numbers by sex, length, and year
    3darray natlength_mat(1,nSXs,styr,endyr,1,nZBs)           // Mature numbers by sex, length, and year
    3darray natl_new_fishtime(1,nSXs,styr,endyr,1,nZBs)           // Numbers-at-length (new shell)
    3darray natl_old_fishtime(1,nSXs,styr,endyr,1,nZBs)           // Numbers-at-length (old shell)    
    3darray natl_inew_fishtime(1,nSXs,styr,endyr,1,nZBs)          // Numbers-at-length (immature new shell)          
    3darray natl_iold_fishtime(1,nSXs,styr,endyr,1,nZBs)          // Numbers-at-length (immature old shell) should be identically 0
    3darray natl_mnew_fishtime(1,nSXs,styr,endyr,1,nZBs)          // Numbers-at-length (mmature new shell) 
    3darray natl_mold_fishtime(1,nSXs,styr,endyr,1,nZBs)          // Numbers-at-length (mmature old shell)
    
    3darray len_len(1,nSXs,1,nZBs,1,nZBs)                    // length to length growth array
    matrix moltp(1,nSXs,1,nZBs)                               // molting probabilities for female, male by length bin 
    matrix moltp_mat(1,nSXs,1,nZBs)                           // molting probs for mature female, male by length bin
    matrix mean_length(1,nSXs,1,nZBs)                         // Predicted post-moult sizes
    
    vector rec_len(1,nZBs)             // Recruitment size frequency
    
    vector predpop_sexr(styr,endyr)                             // Population sex-ratio - output
    
    //changed endyr to endyr-1
    //20150601: changed fm... to f... because these could be fishing mortality OR capture rates
    //fully-selected fishery capture (gmacs option) OR mortality (original) rates
    //20160225: changing to matrices to incorporate sex-specific rates
    matrix fTCF_xy(1,nSXs,styr,endyr-1)    //directed fishery
    matrix fSCF_xy(1,nSXs,styr,endyr-1)    //snow crab fishery
    matrix fRKF_xy(1,nSXs,styr,endyr-1)    //BBRKC fishery
    matrix fGTF_xy(1,nSXs,styr,endyr-1)    //groundfish fisheries
    
    //ratios for extrapolating effort to fishing mortality
    number brSCF  //snow crab fishery
    number brRKF  //BBRKC fishery
                      
    //wts: 20150601: fc's are CAPTURE (NOT mortality) rates (ONLY calculated if using gmacs calculations)
    3darray fcTCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)     // male capture rate in directed fishery                             
    matrix  fcTCFF_yz(styr,endyr-1,1,nZBs)             // female capture (bycatch) rate in directed fishery
    3darray fcSCF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // capture (bycatch) rates in snow crab fishery
    3darray fcRKF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // capture (bycatch) rates in BBRKC fishery 
    3darray fcGTF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // capture (bycatch) rates in trawl fishery
    
    //wts: 20150601: fd's are DISCARD (NOT mortality) rates (ONLY calculated if using gmacs calculations)
    3darray fdTCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)     // male discard rates in directed fishery                      
    matrix  fdTCFF_yz(styr,endyr-1,1,nZBs)             // female discard rates
    3darray fdSCF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // discard rates in snow crab fishery
    3darray fdRKF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // discard rates in BBRKC fishery 
    3darray fdGTF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // discard rates in trawl fishery
    
    //fishing MORTALITY RATES (changed f... to fm... to clarify: 20150601)
    3darray fmTCFR_syz(1,nSCs,styr,endyr-1,1,nZBs)     // Retained fishing mortality on males in directed fishery 
    3darray fmTCFD_syz(1,nSCs,styr,endyr-1,1,nZBs)     // Discard fishing mortality on males in directed fishery 
    3darray fmTCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)     // Total fishing mortality on males in directed fishery                             
    matrix  fmTCFF_yz(styr,endyr-1,1,nZBs)             // Female discard mortality rates in directed fishery
    3darray fmSCF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // discard mortality rates in snow crab fishery
    3darray fmRKF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // discard mortality rates in BBRKC fishery        
    3darray fmGTF_xyz(1,nSXs,styr,endyr-1,1,nZBs)       // discard mortality rates in trawl fishery
    
    4darray fmTOT_xsyz(1,nSXs,1,nSCs,styr,endyr-1,1,nZBs) //TOTAL FISHING MORTALITY! new 20150601
    4darray S_xsyz(1,nSXs,1,nSCs,styr,endyr-1,1,nZBs)     // Survival during fisheries
        
    //new arrays for model fishery catch abundance
    6darray cpN_fyxmsz(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//captured numbers at size in fishery f (calculated if optFM=1)
    6darray dsN_fyxmsz(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discarded catch numbers (NOT mortality) at size in fishery f (calculated if optFM=1)
    6darray tmN_fyxmsz(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//total catch mortality (numbers) at size in fishery f
    6darray dmN_fyxmsz(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discard catch mortality (numbers) at size in fishery f
    4darray rmN_ymsz(styr,endyr-1,1,nMSs,1,nSCs,1,nZBs);                //retained catch mortality (numbers) at size in fishery f
    
    //new arrays for model fishery catch biomass
    5darray cpB_fyxms(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs);//captured biomass in fishery f (calculated if optFM=1)
    5darray dsB_fyxms(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs);//discarded catch biomass (NOT mortality) at size in fishery f (calculated if optFM=1)
    5darray tmB_fyxms(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs);//total catch mortality (biomass) at size in fishery f
    5darray dmB_fyxms(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs);//discard catch mortality (biomass) at size in fishery f
    3darray rmB_yms(styr,endyr-1,1,nMSs,1,nSCs);                //retained catch mortality (biomass) at size in fishery f
    
    //TCF
    matrix catch_lmale_new(styr,endyr-1,1,nZBs)       // Predicted catch mortality (males, total, new shell)
    matrix catch_lmale_old(styr,endyr-1,1,nZBs)       // Predicted catch mortality (males, total, old shell)
    matrix catch_lmale(styr,endyr-1,1,nZBs)           // Catch mortality-in-numbers (males, total)
    matrix catch_female_d_new(styr,endyr-1,1,nZBs)    //female discard mortality, new shell
    matrix catch_female_d_old(styr,endyr-1,1,nZBs)    //female discard mortality, old shell
    matrix catch_female_d(styr,endyr-1,1,nZBs)
    matrix catch_male_ret_new(styr,endyr-1,1,nZBs)    // Predicted catch mortality (males, retained, new shell)
    matrix catch_male_ret_old(styr,endyr-1,1,nZBs)    // Predicted catch mortality (males, retained, old shell)  
    matrix catch_male_ret(styr,endyr-1,1,nZBs)        // Catch mortality-in-numbers (males, retained)
    vector predRetBioMortTCFM_y(styr,endyr-1)      // Retained catch mortality (males, directed)
    vector predTotBioMortTCFM_y(styr,endyr-1)      // Total catch mortality (males, directed)
    vector predDscBioMortTCFM_y(styr,endyr-1)      // Discard catch mortality (males, directed)
    vector predDscBioMortTCFF_y(styr,endyr-1)      // Discard catch mortality (females, directed)
    
    //SCF
    matrix catch_male_snowd_new(styr,endyr-1,1,nZBs)
    matrix catch_male_snowd_old(styr,endyr-1,1,nZBs)
    matrix catch_female_snowd_new(styr,endyr-1,1,nZBs)
    matrix catch_female_snowd_old(styr,endyr-1,1,nZBs)
    matrix catch_male_snowd(styr,endyr-1,1,nZBs)
    matrix catch_female_snowd(styr,endyr-1,1,nZBs)
    matrix predDscBioMortSCF_xy(1,nSXs,styr,endyr-1)
    
    //BBRKC
    matrix catch_male_rkd_new(styr,endyr-1,1,nZBs)
    matrix catch_male_rkd_old(styr,endyr-1,1,nZBs)
    matrix catch_female_rkd_new(styr,endyr-1,1,nZBs)
    matrix catch_female_rkd_old(styr,endyr-1,1,nZBs)
    matrix catch_male_rkd(styr,endyr-1,1,nZBs)
    matrix catch_female_rkd(styr,endyr-1,1,nZBs)
    matrix predDscBioMortRKF_xy(1,nSXs,styr,endyr-1)
    
    //GTF
    matrix catch_trawl_female(styr,endyr-1,1,nZBs)
    matrix catch_trawl_male(styr,endyr-1,1,nZBs)
    vector predDscBioMortGTF_y(styr,endyr-1)            // Trawl catch mortality (in mass and by sex)
    
    //changed endyr to endyr-1   
    matrix effn_fish_ret(1,nSCs,styr,endyr-1)           // Effective sample sizes
    matrix effn_fish_tot(1,nSCs,styr,endyr-1)           // Effective sample sizes
    
    4darray effn_srv1(1,nMSs,1,nSCs,1,nSXs,styr,endyr) // Survey effective sample sizes
    
    // Offsets
    vector offset(1,NUM_LEN_LIKE)  //<-these are constants. should be in DATA_SECTION!!
                                              
    vector objfOut(1,NUM_FOUT);      // objective function components (weighted likelihoods and penalties)
    vector likeOut(1,NUM_FOUT);      //unweighted likelihood components and penalties
    vector wgtsOut(1,NUM_FOUT);      //weights
    
    // Penalties
    number nat_penalty
    number penal_rec                                                  // Recruitment
    number initsmo_penal                                              // Initial size-structure
    number initnum_penal                                              // Low size penalty
    number fpen                                                       // Penalties (misc)
    number sel_50m_penal                                              // Penalties on selectivity
    number af_penal                                                   // Prior on af 
    number srv3q_penalty
    number am_penal                                                   // Prior on am
    number bf_penal                                                   // Prior on bf
    number bm_penal                                                   // Prior on bm
    number penal_sexr                                                 // Penalty of sex ratio of recruitment            
    
    // Likelihood components
    vector zsRetMortBio_TCFR_y(styr,endyr-1)        // z-scores for fit to retained male catch mortality in directed TCF 
    vector zsTotMortBio_TCFM_n(1,nObsDscTCF)        // z-scores for fit to total male catch mortality in directed TCF 
    vector zsDscMortBio_TCFM_n(1,nObsDscTCF)        // z-scores for fit to discard male catch mortality in directed TCF 
    vector zsDscMortBio_TCFF_n(1,nObsDscTCF)        // z-scores for fit to female bycatch mortality in directed TCF 
    matrix zsDscMortBio_SCF_xn(1,nSXs,1,nObsDscSCF) // z-scores for fit to bycatch mortality, by sex, in SCF 
    matrix zsDscMortBio_RKF_xn(1,nSXs,1,nObsDscRKF) // z-scores for fit to bycatch mortality, by sex, in RKFCF 
    vector zsDscMortBio_GTF_n(1,nObsDscGTF)         // z-scores for fit to total bycatch mortality in GTF 
    number lkRetMortBio_TCFR          // likelihood for retained male catch mortality in directed TCF 
    number lkTotMortBio_TCFM          // likelihood for total male catch mortality in directed TCF
    number lkDscMortBio_TCFM          // likelihood for discard male catch mortality in directed TCF
    number lkDscMortBio_TCFF          // likelihood for female bycatch mortality in directed TCF 
    vector lkDscMortBio_SCF_x(1,nSXs) // likelihood for bycatch in SCF, by sex
    vector lkDscMortBio_RKF_x(1,nSXs) // likelihood for bycatch in RKF, by sex
    number lkDscMortBio_GTF           // likelihood for bycatch in GTF
    vector lkZCs(1,NUM_LEN_LIKE)      // likelihood for fishery size compositions
    
    matrix zsSrvMatBio_xn(1,nSXs,1,nObsSrvBio) // z-scores for fits to survey mature biomass
    vector lkSrvMatBio_x(1,nSXs)               // likelihood for survey mature biomass data
//    number surv_like_nowt                                             // Surveys (output)
    3darray lkSrvZCs_msx(1,nMSs,1,nSCs,1,nSXs) // likelihood for survey size composition components
    
    //IMPORTANT CHANGE: was "endyr".  cannot be calculated in endyr. 
    vector mspnum_new_matetime(styr,endyr-1)
    vector mspnum_old_matetime(styr,endyr-1)
    vector mspnum_matetime(styr,endyr-1)
    vector mspbio_old_matetime(styr,endyr-1)  
    vector mspbio_matetime(styr,endyr-1)    
    
    vector fspnum_new_matetime(styr,endyr-1)
    vector fspnum_old_matetime(styr,endyr-1)
    vector fspnum_matetime(styr,endyr-1)
    vector fspbio_new_matetime(styr,endyr-1)
    vector fspbio_matetime(styr,endyr-1)
    
    vector emspbio_matetime(styr,endyr-1)     // Spawning biomass at mating time stuff
    vector efspbio_new_matetime(styr,endyr-1) //"effective" female spawning biomass at mating time
    vector efspbio_matetime(styr,endyr-1)   
    
    //can be calculated in endyr:
    vector mspbio_fishtime(styr,endyr)
    vector fspbio_fishtime(styr,endyr)
    
    matrix maturity_est(1,nSXs,1,nZBs)                              // Maturity-at-length
    
    // Outputs (not in the likelihood function)
    vector num_males_gt101(styr,endyr)
    vector bio_males_gt101(styr,endyr)
    vector modTotLglBio_TCFM_y(styr,endyr-1)     //(IMPORTANT CHANGE: used to be "endyr")
    vector modTotLglNum_TCFM_y(styr,endyr-1)  //(IMPORTANT CHANGE: used to be "endyr")
    matrix pred_tmp_snow(1,nSXs,1,nObsDscSCF)   // was  pred_tmp(1,4,1,nobs_discardc), combining both snow and BBRKC (wts: 20130806)
    matrix pred_tmp_rkc(1,nSXs,1,nObsDscRKF)     //                                                 
    vector tmpp1(1,nZBs)
    vector tmpp2(1,nZBs)
    vector tmpp3(1,nZBs)
    vector tmpp4(1,nZBs)
    number like_mat
    
    sdreport_vector sdrSpBioF(mnYrRecCurr,endyr)                                // Sd_report stuff
    sdreport_vector sdrSpBioM(mnYrRecCurr,endyr)                       //male spawning biomass (1000's t)
    sdreport_vector sdrLegalMales(mnYrRecCurr,endyr)
    sdreport_vector sdrRecEarly(styr,mnYrRecCurr-1)
    sdreport_vector sdrRecF(mnYrRecCurr,endyr)  //was endyr-1
    sdreport_vector sdrRecM(mnYrRecCurr,endyr)  //was endyr-1
    sdreport_number sdrDepletion
    
    sdreport_matrix sdrNatMortImm(1,nSXs,styr,endyr);//natural mortality by year on immatures
    sdreport_matrix sdrNatMortNS(1,nSXs,styr,endyr);
    sdreport_matrix sdrNatMortOS(1,nSXs,styr,endyr);
    
      
    sdreport_vector sdrMMB(styr+1,endyr-1);          //MMB (= 0 in styr, so skip it) at fertilization time
    sdreport_vector sdrLnRecMMB(styr+1,endyr-reclag);//ln(rec[yr+reclag-1]/MMB[yr])  at fertilization
    sdreport_vector sdrLnRec(styr+1,endyr-reclag);   //rec[yr+reclag-1] at fertilization
    sdreport_vector sdrRec(styr+1,endyr-reclag);     //rec[yr+reclag-1] at fertilization
    
    objective_function_value f
    
 LOCAL_CALCS
    CheckFile << "Phase: Moulting probabilities:       " << phase_moltingp << endl;
    CheckFile << "Phase: Logistic selectivity pattern: " << phase_logistic_sel << endl;
    CheckFile << "Phase: Dome-shaped selectivity:      " << phase_fishsel << endl;
    CheckFile << "Phase: Survey selectivity #1         " << survsel1_phase << endl;
    CheckFile << "Phase: Survey selectivity #2         " << survsel_phase << endl;
    CheckFile << "Maturity switch:                     " << maturity_switch << endl;
    cout << "Phase: Moulting probabilities:       " << phase_moltingp     << endl;
    cout << "Phase: Logistic selectivity pattern: " << phase_logistic_sel << endl;
    cout << "Phase: Dome-shaped selectivity:      " << phase_fishsel      << endl;
    cout << "Phase: Survey selectivity #1         " << survsel1_phase     << endl;
    cout << "Phase: Survey selectivity #2         " << survsel_phase      << endl;
    cout << "Maturity switch:                     " << maturity_switch    << endl;
 END_CALCS
    

//========================================================================
//========================================================================
PRELIMINARY_CALCS_SECTION

    if (!usePin){//set initial values from control file inputs 
        //natural mortality multipliers
        Mmult_imat = Mmult_imat_in;
        Mmultm = Mmultm_in;
        Mmultf = Mmultf_in;
        mat_big = mat_big_in;
        
        //recruitment
        pMnLnRecHist = inpMnLnRecHist;
        pMnLnRec     = inpMnLnRec;
    }
    if (jitter) jitterParameters(ptrMC->jitFrac);
 
    // use logistic maturity curve for new shell males instead of fractions by year if switch>0
    // this would be for initial population not probability of moving to mature
    if (maturity_switch > 0) {
        obsAvgMatNS_xz(MALE) = obsPrMatureM_z;
    }
//     CheckFile<<"catch_disc(1) "<<catch_disc(1)<<endl;
//     CheckFile<<"catch_disc(2) "<<catch_disc(2)<<endl;
    CheckFile<<"catch ret numbers "<<endl<<obsRetCatchNum<<endl;
    
    //Compute proportions and offsets for multinomials
    offset.initialize();
    
    //retained males in directed fishery
    obsPrNatZ_TCFR_sn.initialize();
    for (int i=1; i <= nObsRetZCsTCF; i++) {
        double tot = 0.0;
        for (int shell=NEW_SHELL; shell<=OLD_SHELL; shell++) tot += sum(obsRetZCsTCF_snz(shell,i));
        
        obsPrNatZ_TCFR_sn(NEW_SHELL,i) = obsRetZCsTCF_snz(NEW_SHELL,i)*fraction_new_error;
        obsPrNatZ_TCFR_sn(OLD_SHELL,i) = obsRetZCsTCF_snz(OLD_SHELL,i)+(1.-fraction_new_error)*obsRetZCsTCF_snz(NEW_SHELL,i);
        
        for (int shell=NEW_SHELL; shell<=OLD_SHELL; shell++) {obsPrNatZ_TCFR_sn(shell,i) = obsPrNatZ_TCFR_sn(shell,i)/tot;}
        
        dvector sumPropn = obsPrNatZ_TCFR_sn(NEW_SHELL,i) + obsPrNatZ_TCFR_sn(OLD_SHELL,i);
        offset(1)       -= ssRetZCsTCF_sn(NEW_SHELL,i)*sumPropn*log(p_const+sumPropn);//wts: why scale this by only NEW_SHELL?
    }
    CheckFile<<"obsPrNatZ_TCFR_sn(NEW_SHELL)"<<endl<<obsPrNatZ_TCFR_sn(NEW_SHELL)<<endl;
    CheckFile<<"obsPrNatZ_TCFR_sn(OLD_SHELL)"<<endl<<obsPrNatZ_TCFR_sn(OLD_SHELL)<<endl;
    CheckFile<<"offset(1)"<<endl<<offset(1)<<endl;
    
    //all males in directed fishery (discard + retained)
    obsPrNatZ_TCFM_snz.initialize();
    for (int i=1; i <= nObsZCsTCFM; i++) {
        double tot = 0.0;
        for (int shell=NEW_SHELL; shell<=OLD_SHELL; shell++) tot += sum(obsZCsTCFM_snz(shell,i));
        for (int shell=NEW_SHELL; shell<=OLD_SHELL; shell++) obsPrNatZ_TCFM_snz(shell,i) = obsZCsTCFM_snz(shell,i)/tot;
            
        dvector sumPropn = obsPrNatZ_TCFM_snz(NEW_SHELL,i)+obsPrNatZ_TCFM_snz(OLD_SHELL,i);
        offset(2)       -= ssTotZCsTCFM_sn(NEW_SHELL,i)*sumPropn*log(p_const+sumPropn);//wts: why is this only NEW_SHELL?
    }
    CheckFile<<"obsPrNatZ_TCFM_snz(NEW_SHELL)"<<endl<<obsPrNatZ_TCFM_snz(NEW_SHELL)<<endl;
    CheckFile<<"obsPrNatZ_TCFM_snz(OLD_SHELL)"<<endl<<obsPrNatZ_TCFM_snz(OLD_SHELL)<<endl;
    CheckFile<<"offset(2)"<<endl<<offset(2)<<endl;
    
    // Female discards in directed fishery
    obsPrNatZ_TCFF_nz.initialize();
    for (int i=1; i <= nObsZCsTCFF; i++) {
        double tot = sum(obsZCsTCFF_nz(i));
        obsPrNatZ_TCFF_nz(i) = obsZCsTCFF_nz(i)/tot;
        offset(3) -= ssZCsTCFF_n(i)*obsPrNatZ_TCFF_nz(i)*log(p_const+obsPrNatZ_TCFF_nz(i));
    }  
    CheckFile<<"obsPrNatZ_TCFF_nz"<<endl<<obsPrNatZ_TCFF_nz<<endl;
    CheckFile<<"offset(3)"<<endl<<offset(3)<<endl;
    
    // SCF male discards  
    obsPrNatZ_SCF_xnz.initialize();
    for (int i=1;i<=nObsZCsSCF;i++) {
        double tot = 0;
        for(int shell=1;shell<=ALL_SHELL;shell++){ tot += sum(obsZCsSCFM_snz(shell,i));}
        // new and old shell together
        obsPrNatZ_SCF_xnz(MALE,i) = (obsZCsSCFM_snz(NEW_SHELL,i)+obsZCsSCFM_snz(OLD_SHELL,i)+obsZCsSCFM_snz(ALL_SHELL,i)) / tot;
        offset(4) -= ssZCsSCFM_sn(NEW_SHELL,i)*obsPrNatZ_SCF_xnz(MALE,i)*log(p_const+obsPrNatZ_SCF_xnz(MALE,i));
    }
    CheckFile<<"obsPrNatZ_SCF_xnz(MALE)"<<endl<<obsPrNatZ_SCF_xnz(MALE)<<endl;
    CheckFile<<"offset(4)"<<endl<<offset(4)<<endl;    
    // SCF female discards
    for (int i=1;i<=nObsZCsSCF;i++) {
        double tot = sum(obsZCsSCFF_nz(i));//sum over size bins
        obsPrNatZ_SCF_xnz(FEMALE,i) = obsZCsSCFF_nz(i) / tot;
        offset(5) -= ssZCsSCFF_n(i)*(obsPrNatZ_SCF_xnz(FEMALE,i)*log(p_const+obsPrNatZ_SCF_xnz(FEMALE,i)));
    }
    CheckFile<<"obsPrNatZ_SCF_xnz(FEMALE)"<<endl<<obsPrNatZ_SCF_xnz(FEMALE)<<endl;
    CheckFile<<"offset(5)"<<endl<<offset(5)<<endl;
        
    // RKF male discards    wts: added ALL_SHELL
    obsPrNatZ_RKF_xnz.initialize();
    for (int i=1;i<=nObsZCsRKF;i++) {
        double tot = 0;
        for(int shell=1;shell<=ALL_SHELL;shell++){ tot += sum(obsZCsRKFM_snz(shell,i));}
        obsPrNatZ_RKF_xnz(MALE,i) = (obsZCsRKFM_snz(NEW_SHELL,i)+obsZCsRKFM_snz(OLD_SHELL,i)+obsZCsRKFM_snz(ALL_SHELL,i)) / tot;
        offset(6)-=ssZCsRKFM_sn(NEW_SHELL,i)*obsPrNatZ_RKF_xnz(MALE,i)*log(p_const+obsPrNatZ_RKF_xnz(MALE,i));
    }
    CheckFile<<"obsPrNatZ_RKF_xnz(MALE)"<<endl<<obsPrNatZ_RKF_xnz(MALE)<<endl;
    CheckFile<<"offset(6)"<<endl<<offset(6)<<endl;    
    // RKF female discards
    for (int i=1;i<=nObsZCsRKF;i++) {
        double tot  =  sum(obsZCsRKFF_nz(i));
        obsPrNatZ_RKF_xnz(FEMALE,i) = obsZCsRKFF_nz(i) / tot;
        offset(7) -= ssZCsRKFF_n(i)*obsPrNatZ_RKF_xnz(FEMALE,i)*log(p_const+obsPrNatZ_RKF_xnz(FEMALE,i));
    }
    CheckFile<<"obsPrNatZ_RKF_xnz(FEMALE)"<<endl<<obsPrNatZ_RKF_xnz(FEMALE)<<endl;
    CheckFile<<"offset(7)"<<endl<<offset(7)<<endl;
       
    // Trawl bycatch
    obsPrNatZ_GTF_xnz.initialize();
    for (int n=1;n<=nObsZCsGTF;n++) {
        if (optPrNatZ_GTF==0){
            //old way: 
            //1. normalize extended size comp by total counts
            //2. weight sex-specific components of extended size comp by sex-specific ss
            double tot = 0;
            for(int x=1;x<=nSXs;x++) tot += sum(obsZCsGTF_xnz(x,n));
            for(int x=1;x<=nSXs;x++) {
                obsPrNatZ_GTF_xnz(x,n) = obsZCsGTF_xnz(x,n)/tot;
                offset(8) -= ssZCsGTF_xn(x,n)*obsPrNatZ_GTF_xnz(x,n)*log(p_const+obsPrNatZ_GTF_xnz(x,n));
            }
        } else if (optPrNatZ_GTF==1){
            //new way 1: 
            //1. normalize extended size comp by total counts
            //2. weight extended size comp by ss summed over sexes
            double tot = 0;
            for(int x=1;x<=nSXs;x++) tot += sum(obsZCsGTF_xnz(x,n));
            for(int x=1;x<=nSXs;x++){
                obsPrNatZ_GTF_xnz(x,n) = obsZCsGTF_xnz(x,n)/tot;
                offset(8) -= ssObsZCsGTF_n(n)*obsPrNatZ_GTF_xnz(x,n)*log(p_const+obsPrNatZ_GTF_xnz(x,n));
            }
        } else if (optPrNatZ_GTF==2){
            //new way 2: 
            //1. normalize sex-specific size comp by count
            //2. create extended comp by weighting sex-specific parts by sex-specific ss 
            //3. weight extended size comp by ss summed over sexes
            for(int x=1;x<=nSXs;x++) {
                double tot = sum(obsZCsGTF_xnz(x,n));
                obsPrNatZ_GTF_xnz(x,n) = (obsZCsGTF_xnz(x,n)/tot)*(ssZCsGTF_xn(x,n)/ssObsZCsGTF_n(n));
                offset(8) -= ssObsZCsGTF_n(n)*obsPrNatZ_GTF_xnz(x,n)*log(p_const+obsPrNatZ_GTF_xnz(x,n));
            }
        } else {
            cout<<endl<<endl<<"-------------------------------"<<endl;
            cout<<"Option for normalizing NatZ for GTF not recognized: "<<optPrNatZ_GTF<<endl;
            cout<<"Please fix value in Model Control File."<<endl;
            cout<<"-------------------------------"<<endl<<endl<<endl;
            CheckFile<<endl<<endl<<"-------------------------------"<<endl;
            CheckFile<<"Option for normalizing NatZ for GTF not recognized: "<<optPrNatZ_GTF<<endl;
            CheckFile<<"Please fix value in Model Control File."<<endl;
            CheckFile<<"-------------------------------"<<endl<<endl<<endl;
            exit(-1);
        }
    }
    CheckFile<<"obsPrNatZ_GTF_xnz(FEMALE)"<<endl<<obsPrNatZ_GTF_xnz(FEMALE)<<endl;
    CheckFile<<"obsPrNatZ_GTF_xnz(MALE)"  <<endl<<obsPrNatZ_GTF_xnz(MALE)<<endl;
    CheckFile<<"offset(8)"<<endl<<offset(8)<<endl;
    
    obsTotSrvNum.initialize();
    for (int ll=1; ll<=nObsZCsSrv; ll++)
        for (int mat=1;mat<=nMSs;mat++)
            for (int shell=1; shell<=nSCs; shell++)
                for (int sex=1; sex<=nSXs; sex++)
                    obsTotSrvNum(ll) += sum(obsSrvNatZs_msxnz(mat,shell,sex,ll));//normalization is over sex, maturity, shell condition
    CheckFile<<"obsTotSrvNum"<<endl<<obsTotSrvNum<<endl;
    
    // Survey data  
    d5_array obs_p_srv1_len1(1,nMSs,1,nSCs,1,nSXs,1,nObsZCsSrv,1,nZBs);    // temporary array for survey length frequency by maturity state, shell condition, sex
    obs_p_srv1_len1.initialize();
    for (int m=1; m<=nMSs; m++) { //maturity
        for (int s=1; s<=nSCs; s++) { //shell condition
            for (int x=1; x<=nSXs;x++) { //sex
                for (int n=1; n <= nObsZCsSrv; n++){ //year counter
                    if(m==IMMATURE){
                        obs_p_srv1_len1(m,s,x,n) = obsSrvNatZs_msxnz(m,s,x,n)/obsTotSrvNum(n);
                    } else {
                        if( s==NEW_SHELL){
                            obs_p_srv1_len1(m,s,x,n)=fraction_new_error*obsSrvNatZs_msxnz(m,s,x,n)/obsTotSrvNum(n);//NOTE: should be applied to PREDICTED, not OBSERVED
                        } else {
                            obs_p_srv1_len1(m,s,x,n)=(obsSrvNatZs_msxnz(m,s,x,n)+obsSrvNatZs_msxnz(m,IMMATURE,x,n)*(1.-fraction_new_error))/obsTotSrvNum(n);
                        }
                    }
                } //year counter
            } //sex
        } //shell condition
    } //maturity
    
    // use logistic maturity curve for new and old shell male survey data if switch>0 instead of yearly samples
    // old shell already uses ok maturity curve (AEP only applies to OLD SHELL?)
    if (maturity_switch > 0){
        for(int i=1; i <= nObsZCsSrv; i++){
//             tmps = (obs_p_srv1_len1(1,2,2,i)+obs_p_srv1_len1(2,2,2,i));
//             obs_p_srv1_len1(2,2,2,i) = elem_prod(obsAvgMatOS_xz(2),tmps);
//             obs_p_srv1_len1(1,2,2,i) = elem_prod(1.0-obsAvgMatOS_xz(2),tmps);
            dvector tmps = (obs_p_srv1_len1(IMMATURE,OLD_SHELL,MALE,i)+obs_p_srv1_len1(MATURE,OLD_SHELL,MALE,i));
            obs_p_srv1_len1(  MATURE,OLD_SHELL,MALE,i) = elem_prod(obsAvgMatOS_xz(MALE),tmps);
            obs_p_srv1_len1(IMMATURE,OLD_SHELL,MALE,i) = elem_prod(1.0-obsAvgMatOS_xz(MALE),tmps);
        }
    }
    
    // Store results
    obsPrNatZ_Srv_msxnz(IMMATURE) = obs_p_srv1_len1(IMMATURE);
    obsPrNatZ_Srv_msxnz(  MATURE) = obs_p_srv1_len1(  MATURE);
    CheckFile<<"obsPrNatZ_Srv_msxnz(IMMATURE)"<<endl<<obsPrNatZ_Srv_msxnz(IMMATURE)<<endl;
    CheckFile<<"obsPrNatZ_Srv_msxnz(  MATURE)"<<endl<<obsPrNatZ_Srv_msxnz(  MATURE)<<endl;
    
    // for maturity and shell condition together in survey length comp fits
    {
        int sex;
        for (int i=1; i <= nObsZCsSrv; i++){
            sex = MALE;
            dvector vall = obsPrNatZ_Srv_msxnz(IMMATURE,NEW_SHELL,sex,i)+obsPrNatZ_Srv_msxnz(IMMATURE,OLD_SHELL,sex,i);
            offset(9) -= ssObsZCsSrv_msxn(MATURE,NEW_SHELL,sex,i)*vall*log(vall+p_const);//DON'T THINK CORRECT NSAMPLES IS BEING APPLIED HERE!
            dvector val2 = obsPrNatZ_Srv_msxnz(  MATURE,NEW_SHELL,sex,i)+obsPrNatZ_Srv_msxnz(  MATURE,OLD_SHELL,sex,i);
            offset(10) -= ssObsZCsSrv_msxn(MATURE,OLD_SHELL,sex,i)*val2*log(val2+p_const);
            sex = FEMALE;
            dvector val3 = obsPrNatZ_Srv_msxnz(IMMATURE,NEW_SHELL,sex,i)+obsPrNatZ_Srv_msxnz(IMMATURE,OLD_SHELL,sex,i);
            offset(11) -= ssObsZCsSrv_msxn(MATURE,NEW_SHELL,sex,i)*val3*log(val3+p_const);
            dvector val4 = obsPrNatZ_Srv_msxnz(  MATURE,NEW_SHELL,sex,i)+obsPrNatZ_Srv_msxnz(  MATURE,OLD_SHELL,sex,i);
            offset(12) -= ssObsZCsSrv_msxn(MATURE,OLD_SHELL,sex,i)*val4*log(val4+p_const);
        }
    }
    CheckFile<<"offset( 9) = "<<offset( 9)<< endl;  
    CheckFile<<"offset(10) = "<<offset(10)<< endl;  
    CheckFile<<"offset(11) = "<<offset(11)<< endl;  
    CheckFile<<"offset(12) = "<<offset(12)<< endl;  
                
    obsSrvNum_y.initialize();//wts: now initializing this
    obsSrvNum_xyz.initialize();
    obsSrvBio_y.initialize();
    obsSrvBio_xy.initialize();
    obsSrvMatBio_xy.initialize();
    obsSrvMatNum_sxy.initialize();
    for(int i=1;i<=nObsSrvBio;i++) obsSrvNum_y(yrsObsSrvBio_n(i)) = obsSrvNum_n(i);//<-wts : yrsObsZCsSrv_n(i) used to index obsSrvNum_y below [so yrsObsSrvBio_n=yrsObsZCsSrv_n?]
    CheckFile<<"obsSrvNum_y"<<endl<<obsSrvNum_y<<endl;
    
    // Compute survey biomass
    obsSrvNum_xyz.initialize();
    obsSrvBio_xy.initialize();
    obsSrvBio_y.initialize();
    obsSrvMatBio_xy.initialize();
    obsSrvMatNum_sxy.initialize();
    for (int mat=1;mat<=2;mat++) { //maturity status
        for (int shell=1;shell<=2;shell++) { //shell condition
            for (int sex=1;sex<=2;sex++) {          //sex
                for (int i=1; i <= nObsZCsSrv; i++) {
                    obsSrvNum_xyz(sex,yrsObsZCsSrv_n(i)) += obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i));
                    if (sex==FEMALE){
                        obsSrvBio_xy(sex,yrsObsZCsSrv_n(i)) += obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wtf_mz(mat);
                        obsSrvBio_y(yrsObsZCsSrv_n(i))    += obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wtf_mz(mat);
                    } else {
                        obsSrvBio_xy(sex,yrsObsZCsSrv_n(i)) += obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wtm_z;
                        obsSrvBio_y(yrsObsZCsSrv_n(i))    += obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wtm_z;
                    }
                    //  sum to get mature biomass by sex (AEP index is mature animals only?)
                    if(mat==MATURE) {
                        if(sex==FEMALE){
                            obsSrvMatBio_xy(sex,yrsObsZCsSrv_n(i)) += obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wtf_mz(mat);
                        } else {
                            obsSrvMatBio_xy(sex,yrsObsZCsSrv_n(i)) += obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wtm_z;    
                        }
                        obsSrvMatNum_sxy(shell,sex,yrsObsZCsSrv_n(i)) += sum(obsPrNatZ_Srv_msxnz(mat,shell,sex,i)*obsSrvNum_y(yrsObsZCsSrv_n(i)));
                    }
                }
            }
        }
    }
    CheckFile<<"obsSrvNum_xyz"   <<endl<<obsSrvNum_xyz    <<endl;
    CheckFile<<"obsSrvBio_y"  <<endl<<obsSrvBio_y  <<endl;
    CheckFile<<"obsSrvBio_xy" <<endl<<obsSrvBio_xy <<endl;
    CheckFile<<"obsSrvMatBio_xy"<<endl<<obsSrvMatBio_xy<<endl;
    CheckFile<<"obsSrvMatNum_sxy" <<endl<<obsSrvMatNum_sxy <<endl;
    
    // Number of large males
    obsSrvLegalNum_n.initialize();
    obsSrvLegalBio_n.initialize();
    for(int i=1;i<=nObsZCsSrv;i++) {
        // take 1/2 of the 100-104 bin, 
        obsSrvLegalNum_n(i) = 0.5*obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i),23);            //<--hardwired index
        obsSrvLegalBio_n(i) = obsSrvLegalNum_n(i)*wtm_z(23);                            //<--hardwired index
        for(int j=24;j<=nZBs;j++) {                                          //<--hardwired index
            obsSrvLegalNum_n(i) += obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i),j);
            obsSrvLegalBio_n(i) += obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i),j)*wtm_z(j);
        }
    }
    CheckFile<<"obsSrvLegalNum_n"    <<endl<<obsSrvLegalNum_n<<endl;
    CheckFile<<"obsSrvLegalBio_n"<<endl<<obsSrvLegalBio_n<<endl;

    //make sure these calculations get done at least once!!
    // Compute the moulting probabilities
    get_moltingp();
    //  cout<<"done moltingp"<<endl;
    // estimate growth function
    get_growth1();//only option now
    //  cout<<"done growth"<<endl;
    // Set maturity
    get_maturity();
    
    //run population mode with initial parameter values
    runPopMod();
    
    //evaluate the objective function for initial parameter values
    evaluate_the_objective_function();
    
    //write reports for initial model configuration
    ofstream initRept("TCSAM_WTS.init.rep");
    writeReport(initRept);
    initRept.close();
    ofstream initReptToR("TCSAM_WTS.init.R");
    myWriteToR(initReptToR);
    initReptToR.close();
    
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        openMCMCFile();
        cout<<"MCEVAL is on"<<endl;
    }
    
    CheckFile<<"------------------------------------------------------------------------"<<endl<<endl;
    CheckFile<<"Initial Parameter Settings"<<endl;
    writeParameters(CheckFile,0,0);
    ofstream osInitParams("TCSAM_WTS.init_params.csv");
    writeParameters(osInitParams,0,1);
    osInitParams.close();
    
    CheckFile<<"End of PRELIMINARY_CALCS SECTION----------------------------------------"<<endl;
    CheckFile<<"------------------------------------------------------------------------"<<endl<<endl;
    
    cout<<"End of PRELIMINARY_CALCS SECTION----------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------"<<endl<<endl;
    
// ============================================================================
// ============================================================================
PROCEDURE_SECTION                                          //wts: revised

    //run the population model
    runPopMod();

    //evaluate the model fit
    evaluate_the_objective_function();
//    cout<<"evaluate_the_objective_function "<<endl;
    
//    Misc_output();          //for testing
//    writeReport(CheckFile); //for testing
//    writeToR(R_out);             //for testing
    
    if((call_no==1)||(call_no/1000)*1000==call_no) {
        CheckFile<<"call_no = "<<call_no<<endl;
        writeParameters(CheckFile,0,1); 
        CheckFile <<"Likes = "<< objfOut << endl;
        CheckFile <<"phase = "<< current_phase() << " call = " << call_no << " Total Like = " << f << endl;
    }    
//     //for testing
//     cout <<"Likes = "<< objfOut << endl;
//     cout <<"phase = "<< current_phase() << " call = " << call_no << " Total Like = " << f << endl;

    if (sd_phase()) Misc_output();//calculate output for sd_phase quantities 
    
    if (mceval_phase()) {
        Misc_output();
        myWriteMCMCToR(mcmc);
    }
    
//     //testing!!
//     echo<<"----------writing MCMC to R--------------------------------------------"<<endl;
//     ofstream mcmc("TCSAM_WTS.MCMC.R");
//     myWriteMCMCToR(mcmc);
//     mcmc.close();
//     echo<<"----------writing to R--------------------------------------------"<<endl;
//     ofstream myRout("TCSAM_WTS.MyR.R");
//     myWriteToR(myRout);
//     myRout.close();
//     exit(1);
    
    //for testing
//     CheckFile<<"fTCF_xy       = "<<endl<<tb<<fTCF_xy<<endl;
//     CheckFile<<"fmortdf     = "<<endl<<tb<<fmortdf<<endl;
//     CheckFile<<"fSCF_xy = "<<endl<<tb<<fSCF_xy<<endl;
//     CheckFile<<"fRKF_xy   = "<<endl<<tb<<fRKF_xy<<endl;
//     CheckFile<<"fGTF_xy      = "<<endl<<tb<<fGTF_xy<<endl;
//     CheckFile<<"selTCFM(NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<selTCFM(NEW_SHELL,iy)<<endl;
//     CheckFile<<"selTCFM(OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<selTCFM(OLD_SHELL,iy)<<endl;

//     CheckFile<<"retFcn(NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<retFcn(NEW_SHELL,iy)<<endl;
//     CheckFile<<"retFcn(OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<retFcn(OLD_SHELL,iy)<<endl;
//     CheckFile<<"selSCF(1-3,FEMALE) = "<<endl<<tb<<selSCF(1,FEMALE)<<endl<<tb<<selSCF(2,FEMALE)<<endl<<tb<<selSCF(3,FEMALE)<<endl;
//     CheckFile<<"selSCF(1-3,  MALE) = "<<endl<<tb<<selSCF(1,  MALE)<<endl<<tb<<selSCF(2,  MALE)<<endl<<tb<<selSCF(3,  MALE)<<endl;
//     CheckFile<<"selRKF(1-3,FEMALE) = "<<endl<<tb<<selRKF(1,FEMALE)<<endl<<tb<<selRKF(2,FEMALE)<<endl<<tb<<selRKF(3,FEMALE)<<endl;
//     CheckFile<<"selRKF(1-3,  MALE) = "<<endl<<tb<<selRKF(1,  MALE)<<endl<<tb<<selRKF(2,  MALE)<<endl<<tb<<selRKF(3,  MALE)<<endl;
//     CheckFile<<"selGTF(1-3,FEMALE) = "<<endl<<tb<<selGTF(1,FEMALE)<<endl<<tb<<selGTF(2,FEMALE)<<endl<<tb<<selGTF(3,FEMALE)<<endl;
//     CheckFile<<"selGTF(1-3,  MALE) = "<<endl<<tb<<selGTF(1,  MALE)<<endl<<tb<<selGTF(2,  MALE)<<endl<<tb<<selGTF(3,  MALE)<<endl;

//     CheckFile<<"fmTCFF_yz = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmTCFF_yz(iy)<<endl;
//     CheckFile<<"fmGTF_xyz(FEMALE) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmGTF_xyz(FEMALE,iy)<<endl;
//     CheckFile<<"fmGTF_xyz(  MALE) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmGTF_xyz(  MALE,iy)<<endl;
//     CheckFile<<"fmSCF_xyz(FEMALE) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmSCF_xyz(FEMALE,iy)<<endl;
//     CheckFile<<"fmSCF_xyz(  MALE) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmSCF_xyz(  MALE,iy)<<endl;
//     CheckFile<<"fmRKF_xyz(FEMALE) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmRKF_xyz(FEMALE,iy)<<endl;
//     CheckFile<<"fmRKF_xyz(  MALE) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmRKF_xyz(  MALE,iy)<<endl;
//     CheckFile<<"fmTCFM_syz(NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmTCFM_syz(NEW_SHELL,iy)<<endl;
//     CheckFile<<"fmTCFM_syz(OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmTCFM_syz(OLD_SHELL,iy)<<endl;
//     CheckFile<<"fmTCFR_syz(NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmTCFR_syz(NEW_SHELL,iy)<<endl;
//     CheckFile<<"fmTCFR_syz(OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<fmTCFR_syz(OLD_SHELL,iy)<<endl;
//     CheckFile<<"S_xsyz(FEMALE,NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<S_xsyz(FEMALE,NEW_SHELL,iy)<<endl;
//     CheckFile<<"S_xsyz(FEMALE,OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<S_xsyz(FEMALE,OLD_SHELL,iy)<<endl;
//     CheckFile<<"S_xsyz(  MALE,NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<S_xsyz(  MALE,NEW_SHELL,iy)<<endl;
//     CheckFile<<"S_xsyz(  MALE,OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<S_xsyz(  MALE,OLD_SHELL,iy)<<endl;
//     CheckFile<<"End of PROCEDURE SECTION----------------------------------------"<<endl;
//     CheckFile<<"------------------------------------------------------------------------"<<endl<<endl;
//     if (call_no>0) exit(-1);

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION runPopMod
    // Update growth (if the parameters are being estimated)
    if (active(moltp_af) || active(moltp_bf) || active(moltp_am) || active(moltp_bm) || active(moltp_ammat) || active(moltp_bmmat)) {
        get_moltingp();
    }
    
    // growth estimated in prelimn calcs if growth parameters estimated in the model
    // then will redo growth matrix, otherwise not
    if(active(am1) || active(bm1) || active(af1) || active(bf1) || active(growth_beta)) {
        get_growth1();//only option now
    }
    
    get_maturity();
//     cout<<" maturity "<<endl;
    get_selectivity();
//     cout<<" selectivity "<<endl;
    get_mortality();
//     cout<<" mortality "<<endl;
    get_numbers_at_len();
//     cout<<" n at len "<<endl;
    get_catch_at_len();
//     cout<<" catch at len "<<endl;
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION void writeParameters(ofstream& os,int toR, int willBeActive)                        //wts: new
    os<<"index, phase, idx.mn, idx.mx, min, max, value, name, type"<<endl;
    wts::writeParameter(os,af1,toR,willBeActive);      
    wts::writeParameter(os,bf1,toR,willBeActive);      
    wts::writeParameter(os,am1,toR,willBeActive);      
    wts::writeParameter(os,bm1,toR,willBeActive);      
    
    wts::writeParameter(os,growth_beta,toR,willBeActive);      
    wts::writeParameter(os,Mmult_imat,toR,willBeActive);      
    wts::writeParameter(os,Mmultm,toR,willBeActive);      
    wts::writeParameter(os,Mmultf,toR,willBeActive);      
    wts::writeParameter(os,mat_big,toR,willBeActive);      
    wts::writeParameter(os,alpha1_rec,toR,willBeActive);      
    wts::writeParameter(os,beta_rec,toR,willBeActive);      
    
    wts::writeParameter(os,moltp_af,toR,willBeActive);      
    wts::writeParameter(os,moltp_bf,toR,willBeActive);      
    wts::writeParameter(os,moltp_am,toR,willBeActive);      
    wts::writeParameter(os,moltp_bm,toR,willBeActive);      
    wts::writeParameter(os,moltp_ammat,toR,willBeActive);      
    wts::writeParameter(os,moltp_bmmat,toR,willBeActive);      
    
    wts::writeParameter(os,pMnLnRec,toR,willBeActive);      
    wts::writeParameter(os,pRecDevs,toR,willBeActive);      
    wts::writeParameter(os,pMnLnRecHist,toR,willBeActive); 
    wts::writeParameter(os,pRecDevsHist,toR,willBeActive); 
    
    wts::writeParameter(os,pAvgLnF_TCF,toR,willBeActive);   
    wts::writeParameter(os,pF_DevsTCF,toR,willBeActive);    
    wts::writeParameter(os,pAvgLnF_GTF,toR,willBeActive);   
    wts::writeParameter(os,pF_DevsGTF,toR,willBeActive);    
    wts::writeParameter(os,pAvgLnF_SCF,toR,willBeActive);   
    wts::writeParameter(os,pF_DevsSCF,toR,willBeActive);    
    wts::writeParameter(os,pAvgLnF_RKF,toR,willBeActive);   
    wts::writeParameter(os,pF_DevsRKF,toR,willBeActive);    
    
    wts::writeParameter(os,fish_slope_mn,toR,willBeActive);   
    wts::writeParameter(os,log_avg_sel50_mn,toR,willBeActive);    
    wts::writeParameter(os,log_sel50_dev_mn,toR,willBeActive);   
    
    wts::writeParameter(os,fish_fit_slope_mn1,toR,willBeActive);   
    wts::writeParameter(os,fish_fit_sel50_mn1,toR,willBeActive);    
    wts::writeParameter(os,fish_fit_slope_mn2,toR,willBeActive);   
    wts::writeParameter(os,fish_fit_sel50_mn2,toR,willBeActive);    
    
    wts::writeParameter(os,fish_slope_1,toR,willBeActive);   
    wts::writeParameter(os,fish_sel50_1,toR,willBeActive);    
    
    wts::writeParameter(os,fish_slope_yr_3,toR,willBeActive);   
    wts::writeParameter(os,log_avg_sel50_3,toR,willBeActive);    
    wts::writeParameter(os,log_sel50_dev_3,toR,willBeActive);    
    
    wts::writeParameter(os,fish_slope_mn2,toR,willBeActive);   
    wts::writeParameter(os,fish_sel50_mn2,toR,willBeActive);    
    
    wts::writeParameter(os,selTCFF_slp,toR,willBeActive);   
    wts::writeParameter(os,selTCFF_z50,toR,willBeActive);    
    
    wts::writeParameter(os,selSCFF_slpA1,toR,willBeActive);   
    wts::writeParameter(os,selSCFF_z50A1,toR,willBeActive);    
    
    wts::writeParameter(os,selSCFF_slpA2,toR,willBeActive);   
    wts::writeParameter(os,selSCFF_z50A2,toR,willBeActive);    
    
    wts::writeParameter(os,selSCFF_slpA3,toR,willBeActive);   
    wts::writeParameter(os,selSCFF_z50A3,toR,willBeActive);    
    
    wts::writeParameter(os,selSCFM_slpA1,toR,willBeActive);   
    wts::writeParameter(os,selSCFM_z50A1,toR,willBeActive);    
    wts::writeParameter(os,selSCFM_slpD1,toR,willBeActive);   
    wts::writeParameter(os,selSCFM_lnZ50D1,toR,willBeActive);    
    
    wts::writeParameter(os,selSCFM_slpA2,toR,willBeActive);   
    wts::writeParameter(os,selSCFM_z50A2,toR,willBeActive);    
    wts::writeParameter(os,selSCFM_slpD2,toR,willBeActive);   
    wts::writeParameter(os,selSCFM_lnZ50D2,toR,willBeActive);    
    
    wts::writeParameter(os,selSCFM_slpA3,toR,willBeActive);   
    wts::writeParameter(os,selSCFM_z50A3,toR,willBeActive);    
    wts::writeParameter(os,selSCFM_slpD3,toR,willBeActive);   
    wts::writeParameter(os,selSCFM_lnZ50D3,toR,willBeActive);    
    
    wts::writeParameter(os,selRKFF_slpA1,toR,willBeActive);   
    wts::writeParameter(os,selRKFF_z50A1,toR,willBeActive);    
    wts::writeParameter(os,selRKFF_slpA2,toR,willBeActive);   
    wts::writeParameter(os,selRKFF_z50A2,toR,willBeActive);    
    wts::writeParameter(os,selRKFF_slpA3,toR,willBeActive);   
    wts::writeParameter(os,selRKFF_z50A3,toR,willBeActive);    
    
    wts::writeParameter(os,selRKFM_slpA1,toR,willBeActive);   
    wts::writeParameter(os,selRKFM_z50A1,toR,willBeActive);    
    wts::writeParameter(os,selRKFM_slpA2,toR,willBeActive);   
    wts::writeParameter(os,selRKFM_z50A2,toR,willBeActive);    
    wts::writeParameter(os,selRKFM_slpA3,toR,willBeActive);   
    wts::writeParameter(os,selRKFM_z50A3,toR,willBeActive);    
    
    wts::writeParameter(os,selGTFF_slpA1,toR,willBeActive);   
    wts::writeParameter(os,selGTFF_z50A1,toR,willBeActive);   
    wts::writeParameter(os,selGTFF_slpA2,toR,willBeActive);   
    wts::writeParameter(os,selGTFF_z50A2,toR,willBeActive);   
    wts::writeParameter(os,selGTFF_slpA3,toR,willBeActive);   
    wts::writeParameter(os,selGTFF_z50A3,toR,willBeActive);   
    
    wts::writeParameter(os,selGTFM_slpA1,toR,willBeActive);   
    wts::writeParameter(os,selGTFM_z50A1,toR,willBeActive);   
    wts::writeParameter(os,selGTFM_slpA2,toR,willBeActive);   
    wts::writeParameter(os,selGTFM_z50A2,toR,willBeActive);   
    wts::writeParameter(os,selGTFM_slpA3,toR,willBeActive);   
    wts::writeParameter(os,selGTFM_z50A3,toR,willBeActive);   
    
    wts::writeParameter(os,srv2_q,toR,willBeActive);       
    wts::writeParameter(os,srv2_seldiff,toR,willBeActive); 
    wts::writeParameter(os,srv2_sel50,toR,willBeActive);   

    wts::writeParameter(os,srv2a_q,toR,willBeActive);       
    wts::writeParameter(os,srv2a_seldiff,toR,willBeActive); 
    wts::writeParameter(os,srv2a_sel50,toR,willBeActive);   

    wts::writeParameter(os,srv3_q,toR,willBeActive);       
    wts::writeParameter(os,srv3_seldiff,toR,willBeActive); 
    wts::writeParameter(os,srv3_sel50,toR,willBeActive);   

    wts::writeParameter(os,matestf,toR,willBeActive); 
    wts::writeParameter(os,matestm,toR,willBeActive); 
    
    wts::writeParameter(os,srv2_qFem,toR,willBeActive);      
    wts::writeParameter(os,srv2_seldiff_f,toR,willBeActive); 
    wts::writeParameter(os,srv2_sel50_f,toR,willBeActive);   
    
    wts::writeParameter(os,srv2a_qFem,toR,willBeActive);      
    wts::writeParameter(os,srv2a_seldiff_f,toR,willBeActive); 
    wts::writeParameter(os,srv2a_sel50_f,toR,willBeActive);   
    
    wts::writeParameter(os,srv3_qFem,toR,willBeActive);      
    wts::writeParameter(os,srv3_seldiff_f,toR,willBeActive); 
    wts::writeParameter(os,srv3_sel50_f,toR,willBeActive);
    
    wts::writeParameter(os,prNewShellRecruits,toR,willBeActive);
    
    wts::writeParameter(os,pAvgLnF_TCFF,toR,willBeActive);
    wts::writeParameter(os,pAvgLnF_SCFF,toR,willBeActive);
    wts::writeParameter(os,pAvgLnF_RKFF,toR,willBeActive);
    wts::writeParameter(os,pAvgLnF_GTFF,toR,willBeActive);
    
// ----------------------------------------------------------------------
FUNCTION void jitterParameters(double fac)   //wts: new 2014-05-10
    cout<<"starting jitterParameters"<<endl;
    
    af1 = wts::jitterParameter(af1,fac,rng);     // Female growth-increment
    bf1 = wts::jitterParameter(bf1,fac,rng);     // Female growth-increment
    am1 = wts::jitterParameter(am1,fac,rng);     // Male growth-increment
    bm1 = wts::jitterParameter(bm1,fac,rng);     // Male growth-increment
    
    growth_beta = wts::jitterParameter(growth_beta,fac,rng); // Growth beta                                //this is NOT estimated (why?)
    Mmult_imat  = wts::jitterParameter(Mmult_imat,fac,rng);  // natural mortality multiplier for immature females and males
    Mmultm      = wts::jitterParameter(Mmultm,fac,rng);      // natural mortality multiplier for mature new and old shell male
    Mmultf      = wts::jitterParameter(Mmultf,fac,rng);      // natural mortality multiplier for mature new and old shell female
    mat_big     = wts::jitterParameter(mat_big,fac,rng);     // mult. on 1980-1984 M for mature males and females                     
    alpha1_rec  = wts::jitterParameter(alpha1_rec,fac,rng);  // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    beta_rec    = wts::jitterParameter(beta_rec,fac,rng);    // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    
    moltp_af = wts::jitterParameter(moltp_af,fac,rng);       // paramters for logistic function molting   //this is NOT estimated (why?)
    moltp_bf = wts::jitterParameter(moltp_bf,fac,rng);       // female                                    //this is NOT estimated (why?)
    moltp_am = wts::jitterParameter(moltp_am,fac,rng);       // paramters for logistic function molting   //this is NOT estimated (why?)
    moltp_bm = wts::jitterParameter(moltp_bm,fac,rng);       // immature males                            //this is NOT estimated (why?)
    moltp_ammat = wts::jitterParameter(moltp_ammat,fac,rng);    // logistic molting prob for mature males
    moltp_ammat = wts::jitterParameter(moltp_bmmat,fac,rng);    // logistic molting prob for mature males
    
    pMnLnRec = wts::jitterParameter(pMnLnRec,fac,rng);         // Mean log-scale recruitment mnYrRecCurr+ (males, females are equal)
    pRecDevs = wts::jitterParameter(pRecDevs,0.1*fac,rng);     // Deviations about mean recruitment mnYrRecCurr+ (IMPORTANT CHANGE: used to be "endyr-1")
    pMnLnRecHist = wts::jitterParameter(pMnLnRecHist,fac,rng);    // Mean log-scale recruitment in early phase (pre-mnYrRecCurr)
    pRecDevsHist = wts::jitterParameter(pRecDevsHist,0.1*fac,rng);// Deviations about logscale mean recruitment in early phase (pre-mnYrRecCurr)
    
    pAvgLnF_TCF = wts::jitterParameter(pAvgLnF_TCF,fac,rng);           //log-scale mean directed fishing mortality
    pF_DevsTCF  = wts::jitterParameter(pF_DevsTCF,0.1*fac,rng);//log-scale directed fishing mortality devs IMPORTANT CHANGE: USED TO BE "1966,endyr-12"
    pAvgLnF_GTF = wts::jitterParameter(pAvgLnF_GTF,fac,rng);   // fishing mortality (trawl)
    pF_DevsGTF  = wts::jitterParameter(pF_DevsGTF,0.1*fac,rng);// trawl fishery f-devs       (IMPORTANT CHANGE: used to be "endyr") 1973 seems OK
    pAvgLnF_SCF = wts::jitterParameter(pAvgLnF_SCF,fac,rng);   // fishing mortality snow crab fishery discards
    pF_DevsSCF  = wts::jitterParameter(pF_DevsSCF,0.1*fac,rng);// snow crab fishery f-devs   (IMPORTANT CHANGE: used to be "endyr")  1992 is OK
    pAvgLnF_RKF = wts::jitterParameter(pAvgLnF_RKF,fac,rng);   // fishing mortality red king crab fishery discards //this is NOT estimated (why?)
    pF_DevsRKF  = wts::jitterParameter(pF_DevsRKF,0.1*fac,rng);//this is NOT estimated (why?)  IMPORTANT CHANGEA: was nObsDscRKF-1.  why -1 in "nObsDscRKF-1"
    
    // Selectivity pattern for males (directed fishery)
    fish_slope_mn = wts::jitterParameter(fish_slope_mn,fac,rng);           //this is NOT estimated (why?)
    log_avg_sel50_mn = wts::jitterParameter(log_avg_sel50_mn,fac,rng);        //this is NOT estimated (why?)
    log_sel50_dev_mn = wts::jitterParameter(log_sel50_dev_mn,fac,rng);//this is NOT estimated (why?)
    
    // Retention function
    // 1981 - 1992
    fish_fit_slope_mn1 = wts::jitterParameter(fish_fit_slope_mn1,fac,rng);
    fish_fit_sel50_mn1 = wts::jitterParameter(fish_fit_sel50_mn1,fac,rng);
    // 2005-endyr  
    fish_fit_slope_mn2 = wts::jitterParameter(fish_fit_slope_mn2,fac,rng);
    fish_fit_sel50_mn2 = wts::jitterParameter(fish_fit_sel50_mn2,fac,rng);
    
    // Directed fishery selectivity pattern for period-1: 1993-1996
    fish_slope_1 = wts::jitterParameter(fish_slope_1,fac,rng);      
    fish_sel50_1 = wts::jitterParameter(fish_sel50_1,fac,rng);
    
    // Directed fishery selectivity pattern changing by year for period-3: 2005-P
    fish_slope_yr_3 = wts::jitterParameter(fish_slope_yr_3,fac,rng);      
    log_avg_sel50_3 = wts::jitterParameter(log_avg_sel50_3,fac,rng);
    log_sel50_dev_3 = wts::jitterParameter(log_sel50_dev_3,0.1*fac,rng);
    
    // for a dome-shaped selex pattern
    fish_slope_mn2 = wts::jitterParameter(fish_slope_mn2,fac,rng);
    fish_sel50_mn2 = wts::jitterParameter(fish_sel50_mn2,fac,rng);
    
    // Female discards
    selTCFF_slp = wts::jitterParameter(selTCFF_slp,fac,rng);
    selTCFF_z50 = wts::jitterParameter(selTCFF_z50,fac,rng);
    
    // snow fishery female discards for period-1: 1989-1996
    selSCFF_slpA1 = wts::jitterParameter(selSCFF_slpA1,fac,rng);
    selSCFF_z50A1 = wts::jitterParameter(selSCFF_z50A1,fac,rng);
    
    // snow fishery female discards for period-2: 1997-2004
    selSCFF_slpA2 = wts::jitterParameter(selSCFF_slpA2,fac,rng);
    selSCFF_z50A2 = wts::jitterParameter(selSCFF_z50A2,fac,rng);
    
    // snow fishery female discards for period-3: 2005-P
    selSCFF_slpA3 = wts::jitterParameter(selSCFF_slpA3,fac,rng);
    selSCFF_z50A3 = wts::jitterParameter(selSCFF_z50A3,fac,rng);
    
    // snow fishery male discards for period-1: 1989-1996
    selSCFM_slpA1   = wts::jitterParameter(selSCFM_slpA1,fac,rng);
    selSCFM_z50A1   = wts::jitterParameter(selSCFM_z50A1,fac,rng);
    selSCFM_slpD1   = wts::jitterParameter(selSCFM_slpD1,fac,rng);
    selSCFM_lnZ50D1 = wts::jitterParameter(selSCFM_lnZ50D1,fac,rng);
    
    // snow fishery male discards for period-2: 1997-2004
    selSCFM_slpA2   = wts::jitterParameter(selSCFM_slpA2,fac,rng);
    selSCFM_z50A2   = wts::jitterParameter(selSCFM_z50A2,fac,rng);
    selSCFM_slpD2   = wts::jitterParameter(selSCFM_slpD2,fac,rng);
    selSCFM_lnZ50D2 = wts::jitterParameter(selSCFM_lnZ50D2,fac,rng);
    
    // snow fishery male discards for period-3: 2005-P
    selSCFM_slpA3   = wts::jitterParameter(selSCFM_slpA3,fac,rng);
    selSCFM_z50A3   = wts::jitterParameter(selSCFM_z50A3,fac,rng);
    selSCFM_slpD3   = wts::jitterParameter(selSCFM_slpD3,fac,rng);
    selSCFM_lnZ50D3 = wts::jitterParameter(selSCFM_lnZ50D3,fac,rng);
    
    // red king fishery female discards

    selRKFF_slpA1 = wts::jitterParameter(selRKFF_slpA1,fac,rng);
    selRKFF_z50A1 = wts::jitterParameter(selRKFF_z50A1,fac,rng);
    selRKFF_slpA2 = wts::jitterParameter(selRKFF_slpA2,fac,rng);
    selRKFF_z50A2 = wts::jitterParameter(selRKFF_z50A2,fac,rng);
    selRKFF_slpA3 = wts::jitterParameter(selRKFF_slpA3,fac,rng);
    selRKFF_z50A3 = wts::jitterParameter(selRKFF_z50A3,fac,rng);
    
    // red king fishery male discards
    selRKFM_slpA1 = wts::jitterParameter(selRKFM_slpA1,fac,rng);
    selRKFM_z50A1 = wts::jitterParameter(selRKFM_z50A1,fac,rng);
    selRKFM_slpA2 = wts::jitterParameter(selRKFM_slpA2,fac,rng);
    selRKFM_z50A2 = wts::jitterParameter(selRKFM_z50A2,fac,rng);
    selRKFM_slpA3 = wts::jitterParameter(selRKFM_slpA3,fac,rng);
    selRKFM_z50A3 = wts::jitterParameter(selRKFM_z50A3,fac,rng);
    
    // Trawl fishery selectivity female, 1973-1987
    selGTFF_slpA1 = wts::jitterParameter(selGTFF_slpA1,fac,rng);
    selGTFF_z50A1 = wts::jitterParameter(selGTFF_z50A1,fac,rng);
    // Trawl fishery selectivity female, 1988-1996
    selGTFF_slpA2 = wts::jitterParameter(selGTFF_slpA2,fac,rng);
    selGTFF_z50A2 = wts::jitterParameter(selGTFF_z50A2,fac,rng);
    // Trawl fishery selectivity female, 1997-P
    selGTFF_slpA3 = wts::jitterParameter(selGTFF_slpA3,fac,rng);
    selGTFF_z50A3 = wts::jitterParameter(selGTFF_z50A3,fac,rng);
    // Trawl fishery selectivity male, 1973-1987
    selGTFM_slpA1 = wts::jitterParameter(selGTFM_slpA1,fac,rng);
    selGTFM_z50A1 = wts::jitterParameter(selGTFM_z50A1,fac,rng);
    // Trawl fishery selectivity male, 1988-1996
    selGTFM_slpA2 = wts::jitterParameter(selGTFM_slpA2,fac,rng);
    selGTFM_z50A2 = wts::jitterParameter(selGTFM_z50A2,fac,rng);
    // Trawl fishery selectivity male, 1997-P
    selGTFM_slpA3 = wts::jitterParameter(selGTFM_slpA3,fac,rng);
    selGTFM_z50A3 = wts::jitterParameter(selGTFM_z50A3,fac,rng);
    //1974 to 1981 
    srv2_q       = wts::jitterParameter(srv2_q,fac,rng);
    srv2_seldiff = wts::jitterParameter(srv2_seldiff,fac,rng);
    srv2_sel50   = wts::jitterParameter(srv2_sel50,fac,rng);
    //1982-86 net change; 1982 first year of 83-112; burn-in period
    srv2a_q       = wts::jitterParameter(srv2a_q,fac,rng);
    srv2a_seldiff = wts::jitterParameter(srv2a_seldiff,fac,rng);
    srv2a_sel50   = wts::jitterParameter(srv2a_sel50,fac,rng);
    //1987-P
    srv3_q       = wts::jitterParameter(srv3_q,fac,rng);
    srv3_seldiff = wts::jitterParameter(srv3_seldiff,fac,rng);
    srv3_sel50   = wts::jitterParameter(srv3_sel50,fac,rng);
    
    matestf = wts::jitterParameter(matestf,fac,rng);
    matestm = wts::jitterParameter(matestm,fac,rng);
    
    srv2_qFem      = wts::jitterParameter(srv2_qFem,fac,rng);
    srv2_seldiff_f = wts::jitterParameter(srv2_seldiff_f,fac,rng);    
    srv2_sel50_f = wts::jitterParameter(srv2_sel50_f,fac,rng);
    
    srv2a_qFem      = wts::jitterParameter(srv2a_qFem,fac,rng);
    srv2a_seldiff_f = wts::jitterParameter(srv2a_seldiff_f,fac,rng);
    srv2a_sel50_f   = wts::jitterParameter(srv2a_sel50_f,fac,rng);
    
    srv3_qFem      = wts::jitterParameter(srv3_qFem,fac,rng);
    srv3_seldiff_f = wts::jitterParameter(srv3_seldiff_f,fac,rng);
    srv3_sel50_f   = wts::jitterParameter(srv3_sel50_f,fac,rng);
    
    
    prNewShellRecruits = wts::jitterParameter(prNewShellRecruits,fac,rng);
    
    pAvgLnF_TCFF = wts::jitterParameter(pAvgLnF_TCFF,fac,rng);
    pAvgLnF_SCFF = wts::jitterParameter(pAvgLnF_SCFF,fac,rng);
    pAvgLnF_RKFF = wts::jitterParameter(pAvgLnF_RKFF,fac,rng);
    pAvgLnF_GTFF = wts::jitterParameter(pAvgLnF_GTFF,fac,rng);
    
    cout<<"finished jittering"<<endl;
//    cout<<"enter 1 to continue >> ";
//    int dummy = 0;
//    cin>>dummy;
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION openMCMCFile                                     //wts: new
    mcmc.open("TCSAM_WTS.MCMC.R", ofstream::out|ofstream::trunc);
    mcmc<<"mcmc<-list("<<endl;
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION closeMCMCFile                                     //wts: new
    mcmc<<"dummy=0)"<<endl;
    mcmc.close();
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION WriteMCMC                                     //wts: checked
    post<<
    // srv1_slope <<","<<
    // srv1_sel50 <<","<<
    fish_slope_mn <<","<<
    fish_sel50_mn <<","<<
    // fish_fit_slope_mn <<","<<
    // fish_fit_sel50_mn <<","<<
    selTCFF_slp <<","<<
    selTCFF_z50 <<","<<
    //selGTFF_slpA <<","<<
    //selGTFF_z50A <<","<<
    endl;

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
FUNCTION get_maturity                                  //wts: revised
    if(active(matestm)){
        maturity_est(FEMALE)       = 1.0;
        maturity_est(FEMALE)(1,16) = mfexp(matestf);//females> length_bins(16) assumed mature
        maturity_est(MALE)         = mfexp(matestm);
    } else{    
        maturity_est(FEMALE) = obsAvgMatNS_xz(FEMALE);
        maturity_est(MALE)   = obsAvgMatNS_xz(MALE);
    }
//     CheckFile<<"maturity_est"<<endl;
//     CheckFile<<maturity_est(1)<<endl;
//     CheckFile<<maturity_est(2)<<endl;
    
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// NOTE: function renamed from get_growth() and growth_switch flag check removed
FUNCTION get_growth1                                   //wts: revised
    
    mean_length(FEMALE)= mfexp(af1)* pow(length_bins,bf1);
    mean_length(MALE)  = mfexp(am1)* pow(length_bins,bm1);
    
    //  CheckFile << af<<" " <<bf<<" "<<am<<" "<<bm<<endl;
//     CheckFile<<"mean_length"<<endl;
//     CheckFile<<mean_length(1)<<endl;
//     CheckFile<<mean_length(2)<<endl;
    
    // using Gamma function for transition matrix
    // devia is the bounds of growth bins to evaluate
    // the gamma function (x) in prop = integral(i1 to i2) g(x|alpha,beta) dx
    // alpha and growth_beta are parameters 
    // alpha is the mean growth increment per molt for some premolt length class
    // alpha = mean growth increment per molt divided by beta
    // beta is the shape parameter - larger beta - more variance 

    double devia;
    dvariable alpha;    
    len_len.initialize();
    for (int sex=FEMALE;sex<=MALE;sex++) {
        for (int ilen=1;ilen<=nZBs;ilen++){    
            // subract the 2.5 from the midpoint of the length bin to get the lower bound
            alpha = (mean_length(sex,ilen)-(length_bins(ilen)-2.5))/growth_beta(sex);
            //    cout<<"alpha = "<<alpha<<endl;
            //    cout<<"growth_beta = "<<growth_beta<<endl;
            //    cout<<"mean_length = "<<mean_length(sex,ilen)<<endl;
             //    truncate growth transition to max=10 bins
            for (int il2=ilen;il2<=ilen+min(10,nZBs-ilen);il2++) {
                devia = length_bins(il2)+2.5-length_bins(ilen);
                len_len(sex,ilen,il2) = pow(devia,(alpha-1.))*exp(-devia/growth_beta(sex));
            }  
            //standardize so each row sums to 1.0
            len_len(sex,ilen) /= sum(len_len(sex,ilen));
        }
    }
//     CheckFile<<"len_len"<<endl;
//     CheckFile<<len_len(1)<<endl;
//     CheckFile<<len_len(2)<<endl;
    
    // Fraction recruiting
    dvariable alpha_rec = alpha1_rec/beta_rec;
    for (int ilen=1;ilen<=nZBs;ilen++) {
        devia = length_bins(ilen)+2.5-length_bins(1);     //wts: think adding 2.5 here is wrong! just want offset from 1st bin (unless these are cutputs)
        rec_len(ilen) = pow(devia,alpha_rec-1.)*exp(-devia/beta_rec);//mfexp(...) doesn't work here?
    }
    rec_len /= sum(rec_len); //standardize so each row sums to 1.0
    
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
FUNCTION get_moltingp                         //wts: revised
    //assuming a declining logistic function
    moltp(FEMALE)=1.0-(1.0/(1.0+mfexp(-1.*moltp_af*(length_bins-moltp_bf))));
    moltp(MALE)  =1.0-(1.0/(1.0+mfexp(-1.*moltp_am*(length_bins-moltp_bm))));
    
    // set molting prob for mature females at 0.0
    moltp_mat(FEMALE)=0.0;
    
    // molting probability for mature males can be zero (or estimated)
    if(phase_moltingp > 0){
        moltp_mat(MALE) = 1.0-(1.0/(1.0+mfexp(-1.*moltp_ammat*(length_bins-moltp_bmmat))));
    } else {
        moltp_mat(MALE) = 0.0;
    }

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
FUNCTION get_selectivity                  //wts: revised
//    cout<<"get_selectivity"<<endl;
    int ii = 1;
    
    selTCFM.initialize();
    retFcn.initialize();
    // logistic selectivity curves
    if(active(log_sel50_dev_mn)) {
        fish_sel50_mn.initialize();
        for(int iy=styr;iy<endyr;iy++){      //used to be iy<=endyr
            // length at 50% selectivity
            if( (iy<=1979) || (1985<=iy && iy<=1987) || (1997<=iy && iy<=2004) || (2010<=iy) )
                fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn);
            else {
                fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn+log_sel50_dev_mn(ii));
                ii=ii+1;
            }
            
            //logistic selectivity curve
            selTCFM(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*fish_slope_mn*(length_bins-fish_sel50_mn(iy))));
            if(iy<=1992)    
                retFcn(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*fish_fit_slope_mn1*(length_bins-fish_fit_sel50_mn1)));
            else
                retFcn(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*fish_fit_slope_mn2*(length_bins-fish_fit_sel50_mn2)));
                
            if(phase_fishsel > 0) {
                //for dome shaped add this part
                selTCFM(NEW_SHELL,iy)=elem_prod(selTCFM(NEW_SHELL,iy),1./(1.+mfexp(fish_slope_mn2*(length_bins-fish_sel50_mn2))));
            }               
            // set new and old selTCFM same
            selTCFM(OLD_SHELL,iy) = selTCFM(NEW_SHELL,iy);
            retFcn(OLD_SHELL,iy)= retFcn(NEW_SHELL,iy);
        }//year loop
    } //<-if(active(log_sel50_dev_mn))
//    cout<<"get_sel: 1"<<endl;
    
    dvariable tmpSel50 = mean(exp(log_avg_sel50_3+log_sel50_dev_3(1,6)));
    for(int iy=styr;iy<=1990;iy++){ 
        selTCFM(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*fish_slope_1*(length_bins-tmpSel50)));    
        retFcn(NEW_SHELL, iy) = 1./(1.+mfexp(-1.*fish_fit_slope_mn1*(length_bins-fish_fit_sel50_mn1)));
    }
//    cout<<"get_sel: 1a"<<endl;
    int ctr = 1;
//    cout<<"max index of log_sel50_dev_3: "<<log_sel50_dev_3.indexmax()<<endl;
    for(int iy=1991;iy<=1996;iy++){ 
        if (hasDirectedFishery(iy)) {
//            cout<<"yr = "<<iy<<".  ctr = "<<ctr<<endl;
            selTCFM(NEW_SHELL,iy)=1./(1.+mfexp(-1.*fish_slope_1*(length_bins-exp(log_avg_sel50_3+log_sel50_dev_3(ctr++)))));//ctr was iy-1990
        } else {
//            cout<<"yr = "<<iy<<".  no fishery."<<endl;
        }
    }
//    cout<<"get_sel: 1b"<<endl;
    //no directed fishery, set 50% selectivity to mean
    for(int iy=1997;iy<endyr;iy++){ 
        if (hasDirectedFishery(iy)) {
//            cout<<"yr = "<<iy<<".  ctr = "<<ctr<<endl;
            selTCFM(NEW_SHELL,iy)=1./(1.+mfexp(-1.*fish_slope_yr_3*(length_bins-exp(log_avg_sel50_3+log_sel50_dev_3(ctr++)))));//ctr was iy-1998
        } else {
//            cout<<"yr = "<<iy<<".  no fishery."<<endl;
            selTCFM(NEW_SHELL,iy)=1./(1.+mfexp(-1.*fish_slope_yr_3*(length_bins-exp(log_avg_sel50_3))));
        }
    }
//    cout<<"get_sel: 1d"<<endl;
    
    for(int iy=styr;iy<=1990;iy++) retFcn(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*fish_fit_slope_mn1*(length_bins-fish_fit_sel50_mn1)));
//    cout<<"get_sel: 1f"<<endl;
    for(int iy=1991;iy<endyr;iy++) retFcn(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*fish_fit_slope_mn2*(length_bins-fish_fit_sel50_mn2)));
//    cout<<"get_sel: 1g"<<endl;
    
    for(int iy=styr;iy<endyr;iy++){      //used to be iy<=endyr            
        //for dome shaped add this part
        if(phase_fishsel > 0) selTCFM(NEW_SHELL,iy) = elem_div(selTCFM(NEW_SHELL,iy),1.0+mfexp(fish_slope_mn2*(length_bins-fish_sel50_mn2)));
        
        // set new and old selTCFM same
        selTCFM(OLD_SHELL,iy) = selTCFM(NEW_SHELL,iy);
        retFcn(OLD_SHELL,iy)  = retFcn(NEW_SHELL,iy);
    }//year loop
//    cout<<"get_sel: 2"<<endl;
    
    // female discards ascending logistic curve 
    selTCFF=1./(1.+mfexp(-1.*selTCFF_slp*(length_bins-selTCFF_z50)));
    
    //  snow fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selSCF(1,FEMALE)=1./(1.+mfexp(-1.*selSCFF_slpA1*(length_bins-selSCFF_z50A1))); 
    selSCF(2,FEMALE)=1./(1.+mfexp(-1.*selSCFF_slpA2*(length_bins-selSCFF_z50A2))); 
    selSCF(3,FEMALE)=1./(1.+mfexp(-1.*selSCFF_slpA3*(length_bins-selSCFF_z50A3))); 
//    cout<<"get_sel: 2a"<<endl;
        
    //  snow fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selSCF(1,MALE)=elem_prod(1./(1.+mfexp(-1.*selSCFM_slpA1*(length_bins-selSCFM_z50A1))),
                                 1./(1.+mfexp(selSCFM_slpD1*(length_bins-(selSCFM_z50A1+mfexp(selSCFM_lnZ50D1))))));
    selSCF(2,MALE)=elem_prod(1./(1.+mfexp(-1.*selSCFM_slpA2*(length_bins-selSCFM_z50A2))),
                                 1./(1.+mfexp(selSCFM_slpD2*(length_bins-(selSCFM_z50A2+mfexp(selSCFM_lnZ50D2))))));
    selSCF(3,MALE)=elem_prod(1./(1.+mfexp(-1.*selSCFM_slpA3*(length_bins-selSCFM_z50A3))),
                                 1./(1.+mfexp(selSCFM_slpD3*(length_bins-(selSCFM_z50A3+mfexp(selSCFM_lnZ50D3))))));
//    cout<<"get_sel: 2b"<<endl;
    
    //  red fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selRKF(1,FEMALE)=1./(1.+mfexp(-1.*selRKFF_slpA1*(length_bins-selRKFF_z50A1))); 
    selRKF(2,FEMALE)=1./(1.+mfexp(-1.*selRKFF_slpA2*(length_bins-selRKFF_z50A2))); 
    selRKF(3,FEMALE)=1./(1.+mfexp(-1.*selRKFF_slpA3*(length_bins-selRKFF_z50A3))); 
    
    //  red fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selRKF(1,MALE)=1./(1.+mfexp(-1.*selRKFM_slpA1*(length_bins-selRKFM_z50A1))); 
    selRKF(2,MALE)=1./(1.+mfexp(-1.*selRKFM_slpA2*(length_bins-selRKFM_z50A2))); 
    selRKF(3,MALE)=1./(1.+mfexp(-1.*selRKFM_slpA3*(length_bins-selRKFM_z50A3))); 
//    cout<<"get_sel: 2c"<<endl;
    
    //  trawl fishery selectivity for 3 time periods, #1 (1973-1987), #2 (1988-1996) and #3 (1997-P)
    selGTF(1,FEMALE)=1./(1.+mfexp(-1.*selGTFF_slpA1*(length_bins-selGTFF_z50A1)));
    selGTF(2,FEMALE)=1./(1.+mfexp(-1.*selGTFF_slpA2*(length_bins-selGTFF_z50A2)));
    selGTF(3,FEMALE)=1./(1.+mfexp(-1.*selGTFF_slpA3*(length_bins-selGTFF_z50A3)));
    
    selGTF(1,MALE)=1./(1.+mfexp(-1.*selGTFM_slpA1*(length_bins-selGTFM_z50A1)));    
    selGTF(2,MALE)=1./(1.+mfexp(-1.*selGTFM_slpA2*(length_bins-selGTFM_z50A2)));
    selGTF(3,MALE)=1./(1.+mfexp(-1.*selGTFM_slpA3*(length_bins-selGTFM_z50A3)));
//    cout<<"get_sel: 2d"<<endl;
        
    selSrv2.initialize();
    selSrv2a.initialize();
    selSrv3.initialize();
//    // somerton and otto curve for survey selectivities
//    if (survsel_phase<0)
//        selSrv3(MALE) = sel_som(1)/(1.+sel_som(2)*mfexp(-1.*sel_som(3)*length_bins));
//    else
//        selSrv3(MALE) = srv3_q*1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50)/(srv3_seldiff)));
//    // this sets time periods 1 and 2 survey selectivities to somerton otto as well
//    if (survsel1_phase < 0)
//        selSrv2(MALE) = selSrv3(MALE);
//    else { 
//        selSrv2(MALE)  = srv2_q*1./(1.+mfexp(-1.*log(19.)*(length_bins-srv2_sel50)/(srv2_seldiff)));
//        selSrv2a(MALE) = srv3_q*1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50)/(srv3_seldiff)));
//    }
//        
//    //set male and female equal unless estimating qFem
//    selSrv2(FEMALE)  = srv2_qFem*1./(1.+mfexp(-1.*log(19.)*(length_bins-srv2_sel50_f)/(srv2_seldiff_f)));
//    selSrv2a(FEMALE) = srv3_qFem*1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50_f)/(srv3_seldiff_f)));
//    selSrv3(FEMALE)  = srv3_qFem*1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50_f)/(srv3_seldiff_f)));
//    cout<<"get_sel: 3"<<endl;
//    
    dvariable maxsel;
    selTCFR.initialize();
    for(int iy=styr;iy<endyr;iy++){          //used to be iy<=endyr
        maxsel = max(selTCFM(NEW_SHELL,iy));
        if(maxsel<max(selTCFM(OLD_SHELL,iy))) maxsel = max(selTCFM(OLD_SHELL,iy));  //wts: is this differentiable??
        for (int shell=NEW_SHELL;shell<=OLD_SHELL;shell++){
            selTCFM(shell,iy) = selTCFM(shell,iy)/maxsel;
            selTCFR(shell,iy) = elem_prod(retFcn(shell,iy),selTCFM(shell,iy));
        }
    }
//    cout<<"get_sel: 4"<<endl;
//    cout<<"done"<<endl;
    
    //new 20150901-->
    //calculate survey selectivities
    selSrv2( MALE) = 1./(1.+mfexp(-1.*log(19.)*(length_bins-srv2_sel50)/(srv2_seldiff)));
    selSrv2a(MALE) = 1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50)/(srv3_seldiff)));
    selSrv3( MALE) = 1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50)/(srv3_seldiff)));
        
    selSrv2( FEMALE) = 1./(1.+mfexp(-1.*log(19.)*(length_bins-srv2_sel50_f)/(srv2_seldiff_f)));
    selSrv2a(FEMALE) = 1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50_f)/(srv3_seldiff_f)));
    selSrv3( FEMALE) = 1./(1.+mfexp(-1.*log(19.)*(length_bins-srv3_sel50_f)/(srv3_seldiff_f)));
//    cout<<"get_sel: 3"<<endl;
    
     if (optFshSel==1){//set logistic selectivity = 1 in largest size bin
        //TCFM and retFcn
        for(int iy=styr;iy<endyr;iy++){
            selTCFM(NEW_SHELL,iy) /= selTCFM(NEW_SHELL,iy,nZBs);
            retFcn( NEW_SHELL,iy) /= retFcn( NEW_SHELL,iy,nZBs);
            selTCFM(OLD_SHELL,iy) /= selTCFM(OLD_SHELL,iy,nZBs);
            retFcn( OLD_SHELL,iy) /= retFcn( OLD_SHELL,iy,nZBs);
        }
        //TCFF
        selTCFF /= selTCFF(nZBs);
        //SCF females (only)
        selSCF(1,FEMALE) /= selSCF(1,FEMALE,nZBs);
        selSCF(2,FEMALE) /= selSCF(2,FEMALE,nZBs);
        selSCF(3,FEMALE) /= selSCF(3,FEMALE,nZBs);
        //RKF
        selRKF(1,  MALE) /= selRKF(1,  MALE,nZBs);
        selRKF(2,  MALE) /= selRKF(2,  MALE,nZBs);
        selRKF(3,  MALE) /= selRKF(3,  MALE,nZBs);
        selRKF(1,FEMALE) /= selRKF(1,FEMALE,nZBs);
        selRKF(2,FEMALE) /= selRKF(2,FEMALE,nZBs);
        selRKF(3,FEMALE) /= selRKF(3,FEMALE,nZBs);
        //GTF
        selGTF(1,  MALE) /= selGTF(1,  MALE,nZBs);
        selGTF(2,  MALE) /= selGTF(2,  MALE,nZBs);
        selGTF(3,  MALE) /= selGTF(3,  MALE,nZBs);
        selGTF(1,FEMALE) /= selGTF(1,FEMALE,nZBs);
        selGTF(2,FEMALE) /= selGTF(2,FEMALE,nZBs);
        selGTF(3,FEMALE) /= selGTF(3,FEMALE,nZBs);
     }
    
    //survey selectivities
     if (optSrvSel==1){//set logistic selectivity = 1 in largest size bin
        selSrv2(   MALE) /= selSrv2(   MALE,nZBs);
        selSrv2a(  MALE) /= selSrv2a(  MALE,nZBs);
        selSrv3(   MALE) /= selSrv3(   MALE,nZBs);
        selSrv2( FEMALE) /= selSrv2( FEMALE,nZBs);
        selSrv2a(FEMALE) /= selSrv2a(FEMALE,nZBs);
        selSrv3( FEMALE) /= selSrv3( FEMALE,nZBs);
    }
    
    if (survsel_phase<0)
        // use somerton and otto curve for survey selectivities
        selSrv3(MALE) = sel_som(1)/(1.+sel_som(2)*mfexp(-1.*sel_som(3)*length_bins));
    else
        //scale survey selectivity by survey catchability
        selSrv3(MALE) *= srv3_q;
    
    if (survsel1_phase < 0)
        // this sets time periods 1 and 2 survey selectivities to somerton otto as well
        selSrv2(MALE) = selSrv3(MALE);
    else { 
        //scale survey selectivity by survey catchability
        selSrv2(MALE)  *= srv2_q;
        selSrv2a(MALE) *= srv3_q;
    }
        
    //scale survey selectivity by survey catchability
    selSrv2(FEMALE)  *= srv2_qFem;
    selSrv2a(FEMALE) *= srv3_qFem;
    selSrv3(FEMALE)  *= srv3_qFem;
//    cout<<"get_sel: 3"<<endl;
    //<--
    
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
FUNCTION get_mortality
    int debug = 0;
    if (debug) cout<<"get_mortality"<<endl;
    int ii;
    int inc;
    
    M_msx(IMMATURE,NEW_SHELL,FEMALE) = baseM_msx(IMMATURE,NEW_SHELL,FEMALE)*Mmult_imat;
    M_msx(IMMATURE,NEW_SHELL,  MALE) = baseM_msx(IMMATURE,NEW_SHELL,  MALE)*Mmult_imat;
    M_msx(IMMATURE,OLD_SHELL,FEMALE) = baseM_msx(IMMATURE,OLD_SHELL,FEMALE)*Mmult_imat;
    M_msx(IMMATURE,OLD_SHELL,  MALE) = baseM_msx(IMMATURE,OLD_SHELL,  MALE)*Mmult_imat;
    M_msx(  MATURE,NEW_SHELL,FEMALE) = baseM_msx(  MATURE,NEW_SHELL,FEMALE)*Mmultf;
    M_msx(  MATURE,NEW_SHELL,  MALE) = baseM_msx(  MATURE,NEW_SHELL,  MALE)*Mmultm;
    M_msx(  MATURE,OLD_SHELL,FEMALE) = baseM_msx(  MATURE,OLD_SHELL,FEMALE)*Mmultf;
    M_msx(  MATURE,OLD_SHELL,  MALE) = baseM_msx(  MATURE,OLD_SHELL,  MALE)*Mmultm;
    if (debug) cout<<"0"<<endl;
    
    //first year retained catch 1965(1966 fishery) no fishery 1985, 1986 or 1997-2004 or 2010-2012
    fTCF_xy.initialize();
    fTCF_xy(MALE)(styr,1964) = 0.05;//was 1965!!
//    cout<<"0a"<<endl;
    int idx = 1;
    for(int iy =1965;iy<endyr;iy++){
        if(hasDirectedFishery(iy)) fTCF_xy(MALE,iy) = mfexp(pAvgLnF_TCF+pF_DevsTCF(idx++));
    }
    fTCF_xy(FEMALE) = fTCF_xy(MALE)*mfexp(pAvgLnF_TCFF);
    if (debug) cout<<"1"<<endl;
    
    // fmortdf=mfexp(log_avg_fmortdf+fmortdf_dev); using overall fmTCFM_syz for females as well as males in directed fishery
    //Fs in snow and BBRKC fishery are scalars need to multiply in projections by retained snow crab/average snow catch * fmTCFM_syz to get fmTCFM_syz.
    //20150601: ratio is now either mortality rate/effort OR fishing capture rate/effort
    dvar_vector f_SCF1 = mfexp(pAvgLnF_SCF+pF_DevsSCF);
    brSCF = mean(f_SCF1)/(mean(effSCF_y(1992,endyr-1)));
    if (debug) cout<<" brSCF = "<<brSCF<<endl; 
    if (debug) cout<<"2"<<endl;
    
    fSCF_xy.initialize();
    fSCF_xy(MALE)(styr,1977)= 0.01;
    for(int iy=1978;iy<=1991;iy++) fSCF_xy(MALE)(iy) = brSCF*effSCF_y(iy);
    fSCF_xy(MALE)(1992,endyr-1) = f_SCF1;
    fSCF_xy(FEMALE) = fSCF_xy(MALE)*mfexp(pAvgLnF_SCFF);
    if (debug) cout<<"3"<<endl;
    
    // need to have the devs 1992 to present
    //20150601: ratio is now either mortality rate/effort OR fishing capture rate/effort
    dvar_vector f_RKF1(1,nObsDscRKF);   //was nObsDscRKF-1
    f_RKF1 = mfexp(pAvgLnF_RKF+pF_DevsRKF);
    brRKF = mean(1-exp(-f_RKF1))/(mean(effRKF_y(yrsObsDscRKF)));
//    cout<<"4a"<<endl;
    fRKF_xy.initialize();
    fRKF_xy(MALE)(styr,1952)= 0.02; //brRKF*mean(rkccatch(1969,1973));
//    cout<<"4a"<<endl;
    for (int iy=1953;iy<=1965;iy++) fRKF_xy(MALE)(iy) = -log(1-brRKF*effRKF_y(iy));//WTS: used to be rkceffortjap(iy)
//    cout<<"4a"<<endl;
    for (int iy=1966;iy<=1972;iy++) fRKF_xy(MALE)(iy) = -log(1-brRKF*effRKF_y(iy));//WTS: used to be effRKF_y(iy)+rkceffortjap(iy)
//    cout<<"4b"<<endl;
    for (int iy=1973;iy<=1991;iy++) fRKF_xy(MALE)(iy) = -log(1-brRKF*effRKF_y(iy));
//    cout<<"4c"<<endl;
    for (int iy=1953;iy<=1991;iy++) if(fRKF_xy(MALE)(iy)< 0.01) fRKF_xy(MALE)(iy) = 0.01;       
//    cout<<"4d"<<endl;
    for (int iy=1984;iy<=1985;iy++) fRKF_xy(MALE)(iy) = 0.0;
//    cout<<"4e"<<endl;
    for (int iy=1994;iy<=1995;iy++) fRKF_xy(MALE)(iy) = 0.0;
//    cout<<"4f"<<endl;
    for (int iy=1;iy<=nObsDscRKF;iy++) {
        int y = yrsObsDscRKF(iy);
        fRKF_xy(MALE)(y) = mfexp(pAvgLnF_RKF+pF_DevsRKF(iy));
    }
    fRKF_xy(FEMALE) = fRKF_xy(MALE)*mfexp(pAvgLnF_RKFF);
    if (debug) cout<<"5"<<endl;
    
    fGTF_xy.initialize();
    for (int iy=styr;iy<=1972;iy++) fGTF_xy(MALE)(iy) = mean(mfexp(pAvgLnF_GTF+pF_DevsGTF));    
    for (int iy=1973;iy<endyr;iy++) fGTF_xy(MALE)(iy) = mfexp(pAvgLnF_GTF+pF_DevsGTF(iy));
    fGTF_xy(FEMALE) = fGTF_xy(MALE)*mfexp(pAvgLnF_GTFF);
    if (debug) cout<<"4"<<endl;
    
    //  cout<<"fTCF_xy "<<fTCF_xy<<endl;
    //  cout<<"fSCF_xy "<<fSCF_xy<<endl;
    //  cout<<"fRKF_xy "<<fRKF_xy<<endl;
    //  cout<<"rkcatch = "<<rkccatch<<endl;
    
    //initialize capture rates
    fcSCF_xyz.initialize();   
    fcRKF_xyz.initialize();   
    fcGTF_xyz.initialize();   
    fcTCFF_yz.initialize();
    fcTCFM_syz.initialize();   
            
    //initialize discard rates
    fdSCF_xyz.initialize();   
    fdRKF_xyz.initialize();   
    fdGTF_xyz.initialize();   
    fdTCFF_yz.initialize();
    fdTCFM_syz.initialize();   
            
    //initialize fishing mortality rates
    fmSCF_xyz.initialize();   
    fmRKF_xyz.initialize();   
    fmGTF_xyz.initialize();   
    fmTCFF_yz.initialize();
    fmTCFM_syz.initialize();   
    fmTCFR_syz.initialize();   
    fmTCFD_syz.initialize();       
    fmTOT_xsyz.initialize();
    
    //initialize fishery survival rates
    S_xsyz.initialize();
    
    //calculate fully-expanded fishery capture and mortality rates
    dvar_matrix tmpSelSCF(1,nSXs,1,nZBs);                                               
    dvar_matrix tmpSelRKF(1,nSXs,1,nZBs);                                               
    dvar_matrix tmpSelGTF(1,nSXs,1,nZBs);                                          
    for (int iy=styr;iy<endyr;iy++) {
        tmpSelGTF.initialize();
        tmpSelSCF.initialize();
        tmpSelRKF.initialize();
        
        //need to set fmTCFM_syz in directed fishery to 0.0 when was closed 1985-1986 and 1997-2004 and?
        //set fmTCFM_syz in red king to 0 when closed 84-85 and 94-95
        
        // test on year for 3 snow selectivity periods
        if (iy<=1996) {
            tmpSelSCF(FEMALE)=selSCF(1,FEMALE);
            tmpSelSCF(MALE)  =selSCF(1,MALE);
        }
        if (1997<=iy && iy<=2004) {
            tmpSelSCF(FEMALE)=selSCF(2,FEMALE);
            tmpSelSCF(MALE)  =selSCF(2,MALE);
        }
        if (2005<=iy) {
            tmpSelSCF(FEMALE)=selSCF(3,FEMALE);
            tmpSelSCF(MALE)  =selSCF(3,MALE);
        }
        
        // test on year for 3 red selectivity periods
        if (iy<=1996) {
            tmpSelRKF(FEMALE)=selRKF(1,FEMALE);
            tmpSelRKF(MALE)  =selRKF(1,MALE);
        }
        if (1997<=iy && iy<=2004) {
            tmpSelRKF(FEMALE)=selRKF(2,FEMALE);
            tmpSelRKF(MALE)  =selRKF(2,MALE);
        }
        if (2005<=iy) {
            tmpSelRKF(FEMALE)=selRKF(3,FEMALE);
            tmpSelRKF(MALE)  =selRKF(3,MALE);
        }
        
        // test on year for 3 trawl selectivity periods
        if (iy<=1986) {
            tmpSelGTF(FEMALE)=selGTF(1,FEMALE);
            tmpSelGTF(MALE)  =selGTF(1,MALE);
        }
        if (1987<=iy && iy<=1996) {
            tmpSelGTF(FEMALE)=selGTF(2,FEMALE);
            tmpSelGTF(MALE)  =selGTF(2,MALE);
        }
        if (1997<=iy) {
            tmpSelGTF(FEMALE)=selGTF(3,FEMALE);
            tmpSelGTF(MALE)  =selGTF(3,MALE);
        }
        
        //20150601: Added option to use gmacs model
        if (optFM==0){//original fishing mortality formulation
            //fishing mortality
            fmSCF_xyz(FEMALE,iy) = tmpSelSCF(FEMALE)*fSCF_xy(FEMALE,iy);   
            fmSCF_xyz(  MALE,iy) = tmpSelSCF(  MALE)*fSCF_xy(  MALE,iy);   
            fmRKF_xyz(FEMALE,iy) = tmpSelRKF(FEMALE)*fRKF_xy(FEMALE,iy);   
            fmRKF_xyz(  MALE,iy) = tmpSelRKF(  MALE)*fRKF_xy(  MALE,iy);   
            fmGTF_xyz(FEMALE,iy) = tmpSelGTF(FEMALE)*fGTF_xy(FEMALE,iy);   
            fmGTF_xyz(  MALE,iy) = tmpSelGTF(  MALE)*fGTF_xy(  MALE,iy);   
            fmTCFF_yz(iy)= selTCFF*fTCF_xy(FEMALE,iy);
            for(int s=NEW_SHELL;s<=OLD_SHELL;s++) {
                fmTCFM_syz(s,iy)     = selTCFM(s,iy)*fTCF_xy(MALE,iy);//total fishing mortality on males in directed fishery       
                fmTCFR_syz(s,iy)     = selTCFR(s,iy)*fTCF_xy(MALE,iy);//retained fishing mortality on males in directed fishery
                fmTCFD_syz(s,iy)     = fmTCFM_syz(s,iy)-fmTCFR_syz(s,iy);//discard fishing mortality on males in directed fishery
                fmTOT_xsyz(  MALE,s,iy) = fmTCFM_syz(s,iy)+fmSCF_xyz(  MALE,iy)+fmRKF_xyz(  MALE,iy)+fmGTF_xyz(  MALE,iy);            
                fmTOT_xsyz(FEMALE,s,iy) = fmTCFF_yz(iy)   +fmSCF_xyz(FEMALE,iy)+fmRKF_xyz(FEMALE,iy)+fmGTF_xyz(FEMALE,iy);//wts: does not depend on s
                S_xsyz(  MALE,s,iy) = mfexp(-1.0*fmTOT_xsyz(  MALE,s,iy));
                S_xsyz(FEMALE,s,iy) = mfexp(-1.0*fmTOT_xsyz(FEMALE,s,iy));//wts: does not depend on s
           }//shell category
            if (debug) cout<<"6"<<endl;
        } else 
        if (optFM==1){//gmacs-style fishing capture/mortality
            //capture rates
            fcSCF_xyz(FEMALE,iy) = tmpSelSCF(FEMALE)*fSCF_xy(FEMALE,iy);   
            fcSCF_xyz(  MALE,iy) = tmpSelSCF(  MALE)*fSCF_xy(  MALE,iy);   
            fcRKF_xyz(FEMALE,iy) = tmpSelRKF(FEMALE)*fRKF_xy(FEMALE,iy);   
            fcRKF_xyz(  MALE,iy) = tmpSelRKF(  MALE)*fRKF_xy(  MALE,iy);   
            fcGTF_xyz(FEMALE,iy) = tmpSelGTF(FEMALE)*fGTF_xy(FEMALE,iy);   
            fcGTF_xyz(  MALE,iy) = tmpSelGTF(  MALE)*fGTF_xy(  MALE,iy);   
            fcTCFF_yz(iy)= selTCFF*fTCF_xy(FEMALE,iy);
            for(int s=NEW_SHELL;s<=OLD_SHELL;s++) fcTCFM_syz(s,iy) = selTCFM(s,iy)*fTCF_xy(MALE,iy);       
        
            //mortality rates
            fmSCF_xyz(FEMALE,iy) = hm_pot*fcSCF_xyz(FEMALE,iy);   
            fmSCF_xyz(  MALE,iy) = hm_pot*fcSCF_xyz(  MALE,iy);   
            fmRKF_xyz(FEMALE,iy) = hm_pot*fcRKF_xyz(FEMALE,iy);   
            fmRKF_xyz(  MALE,iy) = hm_pot*fcRKF_xyz(  MALE,iy);  
            fmGTF_xyz(FEMALE,iy) = hm_trawl*fcGTF_xyz(FEMALE,iy);  
            fmGTF_xyz(  MALE,iy) = hm_trawl*fcGTF_xyz(  MALE,iy);   
            fmTCFF_yz(iy)        = hm_pot*fcTCFF_yz(iy);
            for (int s=NEW_SHELL;s<=OLD_SHELL;s++) { //over new (shell=1) and old (shell=2) shell...
                fmTCFR_syz(s,iy) = elem_prod(retFcn(s,iy),             fcTCFM_syz(s,iy)); //retention rate
                fmTCFD_syz(s,iy) = elem_prod(hm_pot*(1.0-retFcn(s,iy)),fcTCFM_syz(s,iy));  //discard mortality rate
                fmTCFM_syz(s,iy) = fmTCFR_syz(s,iy)+fmTCFD_syz(s,iy);  //total fishing mortality rate
                fmTOT_xsyz(  MALE,s,iy) = fmTCFM_syz(s,iy)+fmSCF_xyz(  MALE,iy)+fmRKF_xyz(  MALE,iy)+fmGTF_xyz(  MALE,iy);            
                fmTOT_xsyz(FEMALE,s,iy) = fmTCFF_yz(iy)   +fmSCF_xyz(FEMALE,iy)+fmRKF_xyz(FEMALE,iy)+fmGTF_xyz(FEMALE,iy);//wts: does not depend on s
                S_xsyz(  MALE,s,iy) = mfexp(-1.0*fmTOT_xsyz(  MALE,s,iy));
                S_xsyz(FEMALE,s,iy) = mfexp(-1.0*fmTOT_xsyz(FEMALE,s,iy));//wts: does not depend on s
            }//shell category
            if (debug) cout<<"7"<<endl;
        } else {
            cout<<"Error!! Fishing mortality model option '"<<optFM<<"' not recognized!!"<<endl;
            cout<<"Please fix in model options file."<<endl;
            exit(-1);
        }
    }//year
    
    if (debug) cout<<"done"<<endl;

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
FUNCTION get_numbers_at_len                                    //wts: revised
//    cout<<"get_numbers_at_len"<<endl;
    dvar_matrix tmpo(1,2,styr,endyr);
    dvariable tmpi,Surv1,Surv2,Surv3,Surv4,Surv5,Surv6;
    
    rec_y.initialize();
    natlength.initialize();
    natlength_inew.initialize();
    natlength_iold.initialize();
    natlength_mnew.initialize();
    natlength_mold.initialize();
    natlength_new.initialize();
    natlength_old.initialize();
    natlength_i.initialize();
    natlength_mat.initialize();
    
    //nAtZ_msxy.initialize();
    rec_y(styr,mnYrRecCurr-1)            = mfexp(pMnLnRecHist);
    rec_y(mnYrRecDevsHist,mnYrRecCurr-1) = mfexp(pMnLnRecHist+pRecDevsHist);
    rec_y(mnYrRecCurr,endyr)             = mfexp(pMnLnRec+pRecDevs);
//    cout<<"1"<<endl;
    
    //numbers at length from styr to endyr
//     cout<<"lyr_mort = "<<lyr_mort<<endl;
//     cout<<"uyr_mort = "<<uyr_mort<<endl;
//     cout<<"mort_switch = "<<mort_switch<<endl;
//     cout<<"mat_big = "<<mat_big<<endl;
    if (sd_phase()){
        for (int x=1;x<=nSXs;x++){
            for (int yr=sdrNatMortImm(x).indexmin();yr<=sdrNatMortImm(x).indexmax();yr++) {
                sdrNatMortImm(x,yr) = M_msx(IMMATURE,NEW_SHELL,x);
                if((lyr_mort<=yr) && (yr<=uyr_mort) && (mort_switch==1)) {
                    sdrNatMortNS(x,yr) = M_msx(MATURE,NEW_SHELL,x)*mat_big(x);
                    sdrNatMortOS(x,yr) = M_msx(MATURE,OLD_SHELL,x)*mat_big(x);
                } else {
                    sdrNatMortNS(x,yr) = M_msx(MATURE,NEW_SHELL,x);
                    sdrNatMortOS(x,yr) = M_msx(MATURE,OLD_SHELL,x);
                }
            }
        }
    }
    
    for (int sex=1;sex<=nSXs;sex++) {  
        //initialize natlength in styr with recruitment in new shell, immature
        natlength_inew(sex,styr) += rec_y(styr)*rec_len;  
        for (int yr=styr;yr<endyr;yr++) {
            
            Surv1 = mfexp(-M_msx(IMMATURE,NEW_SHELL,sex));
            Surv2 = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,sex));
            if((lyr_mort<=yr) && (yr<=uyr_mort) && (mort_switch==1)) {
//                 cout<<yr<<": applying big_mort"<<endl;
                Surv3 = mfexp(-M_msx(MATURE,NEW_SHELL,sex)*mat_big(sex));
                Surv4 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,sex)*mat_big(sex));
                Surv5 = mfexp(-M_msx(MATURE,OLD_SHELL, sex)*mat_big(sex));
                Surv6 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL, sex)*mat_big(sex));
            } else {
                Surv3 = mfexp(-M_msx(MATURE,NEW_SHELL,sex));
                Surv4 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,sex));
                Surv5 = mfexp(-M_msx(MATURE,OLD_SHELL,sex));
                Surv6 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,sex));
            }
//             if (1978<=yr && yr<=1990) {
//                 CheckFile<<"nAtZ("<<sex<<";IMM;NS;"<<yr<<") = "<<natlength_inew(sex,yr)<<endl;
//                 CheckFile<<"nAtZ("<<sex<<";MAT;NS;"<<yr<<") = "<<natlength_mnew(sex,yr)<<endl;
//                 CheckFile<<"nAtZ("<<sex<<";MAT;OS;"<<yr<<") = "<<natlength_mold(sex,yr)<<endl;
//                 CheckFile<<"Surv1 = "<<Surv1<<endl;
//                 CheckFile<<"Surv2 = "<<Surv2<<endl;
//                 CheckFile<<"Surv3 = "<<Surv3<<endl;
//                 CheckFile<<"Surv4 = "<<Surv4<<endl;
//                 CheckFile<<"Surv5 = "<<Surv5<<endl;
//                 CheckFile<<"Surv6 = "<<Surv6<<endl;
//                 CheckFile<<"S_xsyz("<<sex<<";NEW_SHELL;"<<yr<<") = "<<S_xsyz(sex,NEW_SHELL,yr)<<endl;
//                 CheckFile<<"S_xsyz("<<sex<<";OLD_SHELL;"<<yr<<") = "<<S_xsyz(sex,OLD_SHELL,yr)<<endl;
//             }
            //     cout<<" to 1 "<<endl;
            // Numbers advancing to new shell...
            dvar_vector tmp = Surv1*elem_prod(moltp(sex),elem_prod(S_xsyz(sex,NEW_SHELL,yr),natlength_inew(sex,yr)));
            natlength_new(sex,yr+1) =  tmp * len_len(sex);
            
            dvar_vector tmpo = Surv1*elem_prod(moltp(sex),elem_prod(S_xsyz(sex,OLD_SHELL,yr),natlength_iold(sex,yr)));
            natlength_new(sex,yr+1) +=  tmpo * len_len(sex);
            
            natlength_iold(sex,yr+1) = Surv1*(elem_prod(S_xsyz(sex,NEW_SHELL,yr),natlength_inew(sex,yr))+
                                              elem_prod(S_xsyz(sex,OLD_SHELL,yr),natlength_iold(sex,yr))) - tmp-tmpo;
            
            dvar_vector tmpm = Surv3*elem_prod(moltp_mat(sex),elem_prod(S_xsyz(sex,NEW_SHELL,yr),natlength_mnew(sex,yr)));
            natlength_mnew(sex,yr+1) = tmpm * len_len(sex);
            
            dvar_vector tmpmo = Surv5*elem_prod(moltp_mat(sex),elem_prod(S_xsyz(sex,OLD_SHELL,yr),natlength_mold(sex,yr)));
            natlength_mnew(sex,yr+1) +=   tmpmo * len_len(sex);
            
            natlength_mold(sex,yr+1) = Surv3 * elem_prod(S_xsyz(sex,NEW_SHELL,yr),natlength_mnew(sex,yr)) + 
                                       Surv5 * elem_prod(S_xsyz(sex,OLD_SHELL,yr),natlength_mold(sex,yr)) - tmpm-tmpmo;
            
            // this is for estimating the fraction of new shell that move to old shell to fit
            // the survey data that is split by immature and mature
            natlength_mnew(sex,yr+1) += elem_prod(    maturity_est(sex),natlength_new(sex,yr+1));//add new shell that become mature
            natlength_inew(sex,yr+1)  = elem_prod(1.0-maturity_est(sex),natlength_new(sex,yr+1));//new shell that stay immature
            //     cout<<" to 2 "<<endl;
            // add in recruits for next year
            // put all recruits in new shell immature
            natlength_inew(sex,yr+1) += rec_y(yr+1)*rec_len;
            natlength_new(sex,yr+1)   = natlength_inew(sex,yr+1) + natlength_mnew(sex,yr+1);
            natlength_old(sex,yr+1)   = natlength_mold(sex,yr+1) + natlength_iold(sex,yr+1);
            natlength_mat(sex,yr+1)   = natlength_mnew(sex,yr+1) + natlength_mold(sex,yr+1);
            natlength_i(sex,yr+1)     = natlength_inew(sex,yr+1) + natlength_iold(sex,yr+1);
            natlength(sex,yr+1)       = natlength_mat(sex,yr+1)  + natlength_i(sex,yr+1);
            
        }
    }
//    cout<<"2"<<endl;
    
    //n-at-length just before fishing occurs
    for (int yr=styr;yr<endyr;yr++){
        for (int sex=1;sex<=nSXs;sex++) {
            if (mort_switch2) {//correct way of applying this; wts: new 2013-08-28
                natl_inew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,sex))*natlength_inew(sex,yr);
                natl_iold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,OLD_SHELL,sex))*natlength_iold(sex,yr);
                if(lyr_mort<=yr && yr<=uyr_mort && mort_switch==1) {
                    natl_mnew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,sex)*mat_big(sex))*natlength_mnew(sex,yr);
                    natl_mold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,sex)*mat_big(sex))*natlength_mold(sex,yr);
                } else {  
                    natl_mnew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,sex))*natlength_mnew(sex,yr);
                    natl_mold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,sex))*natlength_mold(sex,yr);
                }
                natl_new_fishtime(sex,yr)  = natl_inew_fishtime(sex,yr)+natl_mnew_fishtime(sex,yr); 
                natl_old_fishtime(sex,yr)  = natl_iold_fishtime(sex,yr)+natl_mold_fishtime(sex,yr);
            } else { //2012 way of doing it               
                if(lyr_mort<=yr && yr<=uyr_mort && mort_switch==1) {
                    natl_inew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,sex))*natlength_inew(sex,yr);
                    natl_iold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,OLD_SHELL,sex))*natlength_iold(sex,yr);
                    natl_mnew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,sex)*mat_big(sex))*natlength_mnew(sex,yr);
                    natl_mold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,sex)*mat_big(sex))*natlength_mold(sex,yr);
                    natl_new_fishtime(sex,yr)  = natl_inew_fishtime(sex,yr)+natl_mnew_fishtime(sex,yr); 
                    natl_old_fishtime(sex,yr)  = natl_iold_fishtime(sex,yr)+natl_mold_fishtime(sex,yr);
                }
                natl_inew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,sex))*natlength_inew(sex,yr);
                natl_iold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,OLD_SHELL,sex))*natlength_iold(sex,yr);
                natl_mnew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,sex))*natlength_mnew(sex,yr);
                natl_mold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,sex))*natlength_mold(sex,yr);
                natl_new_fishtime(sex,yr)  = natl_inew_fishtime(sex,yr)+natl_mnew_fishtime(sex,yr); 
                natl_old_fishtime(sex,yr)  = natl_iold_fishtime(sex,yr)+natl_mold_fishtime(sex,yr);
            }
        }
    }
//    cout<<"3"<<endl;
    //assume mdptFshs_y(endyr) = mdptFshs_y(endyr-1)
    {   int yr = endyr;    //don't worry about mort_switch here
        for (int sex=1;sex<=nSXs;sex++) {
            natl_inew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr-1)*M_msx(IMMATURE,NEW_SHELL,sex))*natlength_inew(sex,yr);
            natl_iold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr-1)*M_msx(IMMATURE,OLD_SHELL,sex))*natlength_iold(sex,yr);
            natl_mnew_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr-1)*M_msx(MATURE,NEW_SHELL,sex))*natlength_mnew(sex,yr);
            natl_mold_fishtime(sex,yr) = mfexp(-mdptFshs_y(yr-1)*M_msx(MATURE,OLD_SHELL,sex))*natlength_mold(sex,yr);
            natl_new_fishtime(sex,yr)  = natl_inew_fishtime(sex,yr)+natl_mnew_fishtime(sex,yr); 
            natl_old_fishtime(sex,yr)  = natl_iold_fishtime(sex,yr)+natl_mold_fishtime(sex,yr);
        }
    }
//    cout<<"3a"<<endl;
    
    // predicted survey values 
    for (int yr=styr;yr<=endyr;yr++){
        for(int sex=FEMALE;sex<=MALE;sex++) {
            if (yr<1982)             {modSrvNum_xy(sex,yr) = (natlength(sex,yr)*selSrv2(sex));}  else
            if (1982<=yr && yr<1988) {modSrvNum_xy(sex,yr) = (natlength(sex,yr)*selSrv2a(sex));} else
            if (1988<=yr)            {modSrvNum_xy(sex,yr) = (natlength(sex,yr)*selSrv3(sex));}
        }
    }
//    cout<<"4"<<endl;
    
    dvariable totSrvNum;
    dvar_matrix useSelSrv(1,nSXs,1,nZBs);
    pred_bio.initialize();
    fspbio.initialize();
    mspbio.initialize(); 
    for (int yr=styr;yr<=endyr;yr++) {
        fspbio(yr) = natlength_mat(FEMALE,yr)*wtf_mz(MATURE);//dot product sum
        mspbio(yr) = natlength_mat(  MALE,yr)*wtm_z;        //dot product sum
        
        // Selection pattern
        if (yr<1982)             {useSelSrv = selSrv2;}  else
        if (1982<=yr && yr<1988) {useSelSrv = selSrv2a;} else
        if (1988<=yr)            {useSelSrv = selSrv3;}
        
        modSrvMatBio_xy(FEMALE,yr) = multQ*natlength_mat(FEMALE,yr)*elem_prod(wtf_mz(MATURE),useSelSrv(FEMALE));//dot product sum
        modSrvMatBio_xy(  MALE,yr) = multQ*natlength_mat(  MALE,yr)*elem_prod(wtm_z,         useSelSrv(  MALE));//dot product sum
        
        // this is predicted survey in numbers not biomass-don't adjust by max selectivity 
        for(int sex=FEMALE;sex<=MALE;sex++) modSrvNum_xy(sex,yr) = (natlength(sex,yr)*useSelSrv(sex));
        totSrvNum = modSrvNum_xy(FEMALE,yr) + modSrvNum_xy(MALE,yr);
        if(totSrvNum<0.001) totSrvNum = 1.0;                     //this is non-differentiable, but PROBABLY just means no survey was done
        for(int sex=FEMALE;sex<=MALE;sex++) {
            modSrvPrNatZ_NS_mxyz(IMMATURE,sex,yr) = elem_prod(useSelSrv(sex),natlength_inew(sex,yr))/totSrvNum;
            modSrvPrNatZ_OS_mxyz(IMMATURE,sex,yr) = elem_prod(useSelSrv(sex),natlength_iold(sex,yr))/totSrvNum;
            modSrvPrNatZ_NS_mxyz(  MATURE,sex,yr) = elem_prod(useSelSrv(sex),natlength_mnew(sex,yr))/totSrvNum;
            modSrvPrNatZ_OS_mxyz(  MATURE,sex,yr) = elem_prod(useSelSrv(sex),natlength_mold(sex,yr))/totSrvNum;
        } 
        pred_bio(yr) += natlength_inew(FEMALE,yr)*wtf_mz(IMMATURE)+(natlength_mnew(FEMALE,yr)+natlength_mold(FEMALE,yr))*wtf_mz(MATURE)
                      +(natlength_inew(MALE,  yr)              + natlength_mnew(MALE,  yr)+natlength_mold(MALE,  yr))*wtm_z;
    }  
    //    cout<<" end srv 2"<<endl;
    sdrDepletion = pred_bio(endyr) / pred_bio(styr);
    sdrSpBioF   = fspbio(sdrSpBioF.indexmin(),sdrSpBioF.indexmax());
    sdrSpBioM   = mspbio(sdrSpBioM.indexmin(),sdrSpBioM.indexmax());
//    cout<<"5"<<endl;
    
    // Legal males
    numLegalMales_y.initialize();
    modSrvNumLegal_y.initialize();
    for (int yr=styr;yr<=endyr;yr++) {
        // Selection pattern//
        if (yr<1982)             {useSelSrv = selSrv2;}  else
        if (1982<=yr && yr<1988) {useSelSrv = selSrv2a;} else
        if (1988<=yr)            {useSelSrv = selSrv3;}        
        // legal is >=138mm take half the numbers in the 135-139 bin
        numLegalMales_y(yr)    = 0.5*natlength(MALE,yr,23);
        modSrvNumLegal_y(yr) = 0.5*natlength(MALE,yr,23)*useSelSrv(MALE,23);  //fixed indices; need vector of 0's, 0.5 and 1's to mult here
        for(int j=24;j<=nZBs;j++) {
            numLegalMales_y(yr)    += natlength(MALE,yr,j);
            modSrvNumLegal_y(yr) += natlength(MALE,yr,j)*useSelSrv(MALE,j);
        }
    }
    //  cout<<"6"<<endl;
    
    sdrLegalMales = numLegalMales_y(sdrLegalMales.indexmin(),sdrLegalMales.indexmax());                                       //fixed index
    sdrRecEarly   = rec_y(sdrRecEarly.indexmin(),sdrRecEarly.indexmax());
    sdrRecF       = rec_y(sdrRecF.indexmin(),sdrRecF.indexmax());//was "endyr-1"
    sdrRecM       = rec_y(sdrRecM.indexmin(),sdrRecM.indexmax());//was "endyr-1"
    //  cout<<" to end of number at len "<<endl;
    //  cout<<"done"<<endl;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
FUNCTION get_catch_at_len                  //wts: new version
//    cout<<"get_catch_at_len"<<endl;
    dvar_vector ratio1(1,nZBs);
    dvar_vector ratio2(1,nZBs);
    
    predRetBioMortTCFM_y.initialize();
    predTotBioMortTCFM_y.initialize();
    predDscBioMortTCFM_y.initialize();
    predDscBioMortTCFF_y.initialize();
    predDscBioMortSCF_xy.initialize();
    predDscBioMortRKF_xy.initialize();
    predDscBioMortGTF_y.initialize();
    //   cout<<" to get catch at length "<<endl;
    for (int yr=styr;yr<endyr;yr++){                //(IMPORTANT CHANGE: used to be "endyr")        
        //total male directed catch mortality
        ratio1 = elem_prod(elem_div(fmTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
        tmN_fyxmsz(iTCF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
        tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(MALE,yr));
        tmN_fyxmsz(iTCF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(MALE,yr));
        catch_lmale_new(yr) = tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL);
        catch_lmale_old(yr) = tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL);
        catch_lmale(yr)     = catch_lmale_new(yr)+catch_lmale_old(yr);
        predTotBioMortTCFM_y(yr) = catch_lmale(yr)*wtm_z;//note dot product sum over size bins here
        //         cout<<ratio1<<endl;
        //         cout<<ratio2<<endl;
        
        //retained male directed catch mortality
        ratio1 = elem_prod(elem_div(fmTCFR_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmTCFR_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        rmN_ymsz(yr,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
        rmN_ymsz(yr,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
        rmN_ymsz(yr,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(MALE,yr));
        rmN_ymsz(yr,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(MALE,yr));
        catch_male_ret_new(yr) = rmN_ymsz(yr,IMMATURE,NEW_SHELL) + rmN_ymsz(yr,MATURE,NEW_SHELL); 
        catch_male_ret_old(yr) = rmN_ymsz(yr,IMMATURE,OLD_SHELL) + rmN_ymsz(yr,MATURE,OLD_SHELL);
        catch_male_ret(yr)     = catch_male_ret_new(yr)+catch_male_ret_old(yr);
        predRetBioMortTCFM_y(yr) = catch_male_ret(yr)*wtm_z;
        
        predDscBioMortTCFM_y(yr) = predTotBioMortTCFM_y(yr) - predRetBioMortTCFM_y(yr);
        
        //snow crab discard catch male mortality
        ratio1 = elem_prod(elem_div(fmSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
        tmN_fyxmsz(iSCF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
        tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(MALE,yr));
        tmN_fyxmsz(iSCF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(MALE,yr));
        catch_male_snowd_new(yr) = tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iSCF,yr,MALE,MATURE,NEW_SHELL); 
        catch_male_snowd_old(yr) = tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iSCF,yr,MALE,MATURE,OLD_SHELL);
        catch_male_snowd(yr)     = catch_male_snowd_new(yr)+catch_male_snowd_old(yr);
        predDscBioMortSCF_xy(MALE,yr) = catch_male_snowd(yr)*wtm_z;
        
        //red king crab discard catch male     
        ratio1 = elem_prod(elem_div(fmRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
        tmN_fyxmsz(iRKF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
        tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(MALE,yr));
        tmN_fyxmsz(iRKF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(MALE,yr));
        catch_male_rkd_new(yr) = tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iRKF,yr,MALE,MATURE,NEW_SHELL); 
        catch_male_rkd_old(yr) = tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iRKF,yr,MALE,MATURE,OLD_SHELL);
        catch_male_rkd(yr)     = catch_male_rkd_new(yr)+catch_male_rkd_old(yr);
        predDscBioMortRKF_xy(MALE,yr) = catch_male_rkd(yr)*wtm_z;
        
        //trawl bycatch male (20150601: is this correct?!->rates are indep of shell condition)
        ratio1 = elem_prod(elem_div(fmGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
        tmN_fyxmsz(iGTF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
        tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(MALE,yr));
        tmN_fyxmsz(iGTF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(MALE,yr));
        catch_trawl_male(yr)  = tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iGTF,yr,MALE,MATURE,NEW_SHELL) 
                              + tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iGTF,yr,MALE,MATURE,OLD_SHELL); 
        predDscBioMortGTF_y(yr) += catch_trawl_male(yr)*wtm_z;
        
        //directed tanner discard catch female
        ratio1 = elem_prod(elem_div(fmTCFF_yz(yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmTCFF_yz(yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iTCF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(FEMALE,yr));
        tmN_fyxmsz(iTCF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(FEMALE,yr));
        catch_female_d_new(yr)  = tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iTCF,yr,FEMALE,MATURE,NEW_SHELL); 
        catch_female_d_old(yr)  = tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iTCF,yr,FEMALE,MATURE,OLD_SHELL);
        catch_female_d(yr)      = catch_female_d_new(yr)+catch_female_d_old(yr);
        predDscBioMortTCFF_y(yr)  = catch_female_d(yr)*wtf_mz(MATURE);            
        
        //snow crab discard catch female
        ratio1 = elem_prod(elem_div(fmSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iSCF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(FEMALE,yr));
        tmN_fyxmsz(iSCF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(FEMALE,yr));
        catch_female_snowd_new(yr)  = tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iSCF,yr,FEMALE,MATURE,NEW_SHELL); 
        catch_female_snowd_old(yr)  = tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iSCF,yr,FEMALE,MATURE,OLD_SHELL);
        catch_female_snowd(yr)      = catch_female_snowd_new(yr)+catch_female_snowd_old(yr);
        predDscBioMortSCF_xy(FEMALE,yr) = catch_female_snowd(yr)*wtf_mz(MATURE);
        
        //red king crab discard catch female
        ratio1 = elem_prod(elem_div(fmRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iRKF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(FEMALE,yr));
        tmN_fyxmsz(iRKF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(FEMALE,yr));
        catch_female_rkd_new(yr)   = tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iRKF,yr,FEMALE,MATURE,NEW_SHELL); 
        catch_female_rkd_old(yr)   = tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iRKF,yr,FEMALE,MATURE,OLD_SHELL);
        catch_female_rkd(yr)       = catch_female_rkd_new(yr)+catch_female_rkd_old(yr);
        predDscBioMortRKF_xy(FEMALE,yr) = catch_female_rkd(yr)*wtf_mz(MATURE);
        
        //trawl bycatch female (20150601: explicitly calculating ratio 2 rather than assuming ratio2=ratio1)
        ratio1 = elem_prod(elem_div(fmGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iGTF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
        tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,natl_iold_fishtime(FEMALE,yr));
        tmN_fyxmsz(iGTF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,natl_mold_fishtime(FEMALE,yr));
        catch_trawl_female(yr) = tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iGTF,yr,FEMALE,MATURE,NEW_SHELL) 
                               + tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iGTF,yr,FEMALE,MATURE,OLD_SHELL); 
        predDscBioMortGTF_y(yr) += catch_trawl_female(yr)*wtf_mz(MATURE);//add in females to males for trawl catches
    }    //end of year loop
    
    modPrNatZ_TCFR_syz.initialize();
    modPrNatZ_TCFM_syz.initialize();
    modPrNatZ_TCFF_yz.initialize();
    modPrNatZ_SCF_xyz.initialize();
    modPrNatZ_RKF_xyz.initialize();
    modPrNatZ_GTF_xyz.initialize();
    for (int yr=styr;yr<endyr;yr++) {          //(IMPORTANT CHANGE: used to be "endyr")        
        // Retained catch (males)
        if(sum(catch_male_ret(yr))>0.000001){                         //non-differentiable--does it matter
            dvariable tot = sum(catch_male_ret(yr));
            modPrNatZ_TCFR_syz(NEW_SHELL,yr) = catch_male_ret_new(yr)/tot;
            modPrNatZ_TCFR_syz(OLD_SHELL,yr) = catch_male_ret_old(yr)/tot;
        }
        
        // Total catch (males)
        if(sum(catch_lmale(yr))>0.0000001){                         //non-differentiable--does it matter
            dvariable tot = sum(catch_lmale(yr));
            modPrNatZ_TCFM_syz(NEW_SHELL,yr) = catch_lmale_new(yr)/tot;
            modPrNatZ_TCFM_syz(OLD_SHELL,yr) = catch_lmale_old(yr)/tot;
        }
         
        
        // female discards
        if(sum(catch_female_d(yr))>0.0000001){                         //non-differentiable--does it matter
            modPrNatZ_TCFF_yz(yr) = catch_female_d(yr)/sum(catch_female_d(yr));
        }
        
        // snow crab discards female male
        if(sum(catch_female_snowd(yr))>0.00000001){                    //non-differentiable--does it matter
            modPrNatZ_SCF_xyz(FEMALE,yr) = catch_female_snowd(yr)/sum(catch_female_snowd(yr));
        }
        if(sum(catch_male_snowd(yr))>0.00000001){                       //non-differentiable--does it matter
            modPrNatZ_SCF_xyz(MALE,yr) = catch_male_snowd(yr)/sum(catch_male_snowd(yr));
        }
        
        // red king crab discards female male
        if(sum(catch_male_rkd(yr))>0.00000001){                         //non-differentiable--does it matter
            modPrNatZ_RKF_xyz(MALE,yr) = catch_male_rkd(yr)/sum(catch_male_rkd(yr));
        }
        if(sum(catch_female_rkd(yr))>0.00000001){                       //non-differentiable--does it matter
            modPrNatZ_RKF_xyz(FEMALE,yr) = catch_female_rkd(yr)/sum(catch_female_rkd(yr));
        }
        
        // total trawl selected numbers
        {
            dvariable tot = sum(catch_trawl_male(yr))+sum(catch_trawl_female(yr));
            // Trawl proportions (adds to 1 over sex, shell and length)
            modPrNatZ_GTF_xyz(FEMALE,yr) = catch_trawl_female(yr)/tot;
            modPrNatZ_GTF_xyz(  MALE,yr) = catch_trawl_male(yr)/tot;
        }
    }//yr
    //  cout<<" done catch at len "<<endl;
//    cout<<"done"<<endl;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
FUNCTION evaluate_the_objective_function    //wts: revising
//    cout<<"evaluate the objective funtion"<<endl;
    int yr;
    dvariable nextf;
    dvar_matrix cv_srv1(1,nSXs,styr,endyr);
    dvariable multi;
    
    f.initialize();
    objfOut.initialize();
    likeOut.initialize();
    wgtsOut.initialize();
    //cout<<" to obj func "<<endl;
    // PENALTIES
    // =========
    
    // Constraints on recruitment
    penal_rec.initialize();
    if (active(pRecDevs)) {        
        //recruitment likelihood - norm2 is sum of square values   
        penal_rec = 1.0*like_wght_recf*norm2(pRecDevs); //+ like_wght_rec*norm2(rec_devm);
        //   first difference on recruitment in period-1     
        penal_rec += 1.0*norm2(first_difference(pRecDevsHist));
        
        f += penal_rec; objfOut(1) = penal_rec; likeOut(1) = penal_rec; wgtsOut(1) = 1;
        
//         // Deviations on difference between make and female recruites (not used)
//         penal_sexr.initialize();
//         if (active(rec_devm)) {
// //            for(int i=styr;i<endyr;i++) penal_sexr += like_wght_sexr*square((mnLnRec_x(FEMALE)+pRecDevs(i))-(mnLnRec_x(MALE)+rec_devm(i)));
//             penal_sexr += like_wght_sexr*norm2((mnLnRec_x(FEMALE)+pRecDevs)-(mnLnRec_x(MALE)+rec_devm));
//         }
//         objfOut(2) = penal_sexr; f += penal_sexr; likeOut(2) = penal_sexr; wgtsOut(2) = 1;
    } 
    
    //nat Mort. penalty
    if(active(Mmult_imat)) {
        nat_penalty = 0.5 * square((Mmult_imat - 1.0) / 0.05);                 //hard-wired
        f += nat_penalty; objfOut(3) = nat_penalty; likeOut(3) = nat_penalty; wgtsOut(3) = 1;
    }
    if(active(Mmultm)) {  
        nat_penalty = 0.5 * square((Mmultm - 1.0) / 0.05);                     //hard-wired
        f += nat_penalty; objfOut(4) = nat_penalty; likeOut(4) = nat_penalty; wgtsOut(4) = 1;
    }
    if(active(Mmultf)) {  
        nat_penalty = 0.5 * square((Mmultf - 1.0) / 0.05);                     //hard-wired
        f += nat_penalty; objfOut(5) = nat_penalty; likeOut(5) = nat_penalty; wgtsOut(5) = 1;
    }
    
    //penalty on survey Q
//    if(active(srv3_q) && srv3_qPriorWgt>0) {  
//        //max of underbag at 182.5 mm is 0.873   
//        srv3q_penalty = 0.5 * square((srv3_q - 0.88) / 0.05);                  //hard-wired
//        //    srv3q_penalty = 0.0 * square((srv3_q - 0.88) / 0.05);
//        f += srv3q_penalty; objfOut(6) = srv3q_penalty; likeOut(6) = srv3q_penalty; wgtsOut(6) = 1;
//    }
//    if(active(srv3_qFem) && srv3_qPriorWgt>0) {  
//        //peak of females is at about 80mm underbag is 0.75 at this size - less uncertainty  
//        srv3q_penalty = 0.5 * square((srv3_qFem - 0.88) / 0.05);                //hard-wired
//        //    srv3q_penalty = 0.0 * square((srv3_qFem - 0.88) / 0.05);
//        f += srv3q_penalty; objfOut(7) = srv3q_penalty; likeOut(7) = srv3q_penalty; wgtsOut(7) = 1;
//    }
    if(active(srv3_q) && srv3_qPriorWgt>=0) {  
        //max of underbag at 182.5 mm is 0.873   
        srv3q_penalty = 0.5 * square((srv3_q - srv3_qPriorMean) / srv3_qPriorStD);
        f += srv3_qPriorWgt*srv3q_penalty; objfOut(6) = srv3_qPriorWgt*srv3q_penalty; likeOut(6) = srv3q_penalty; wgtsOut(6) = srv3_qPriorWgt;
    }
    if(active(srv3_qFem) && srv3_qFemPriorWgt>=0) {  
        //peak of females is at about 80mm underbag is 0.75 at this size - less uncertainty  
        srv3q_penalty = 0.5 * square((srv3_qFem - srv3_qFemPriorMean) / srv3_qFemPriorStD);
        f += srv3_qFemPriorWgt*srv3q_penalty; objfOut(7) = srv3_qFemPriorWgt*srv3q_penalty; likeOut(7) = srv3q_penalty; wgtsOut(7) = srv3_qFemPriorWgt;
    }
    
    // bayesian part - likelihood on growth parameters af,am,bf,bm
    // not used in this case
    af_penal = 0; bf_penal = 0; am_penal = 0; bm_penal = 0;
    if(active(af1)) {  
        af_penal = 0.5 * square((af1 - 0.56560241)/0.1);                  //hard-wired
        f += af_penal; objfOut(8) = af_penal; likeOut(8) = af_penal; wgtsOut(8) = 1;
    }
    if(active(bf1)) {  
        bf_penal = 0.5 * square((bf1 - 0.9132661)/0.025);                    //hard-wired
        f += bf_penal; objfOut(9) = bf_penal; likeOut(9) = bf_penal; wgtsOut(9) = 1;
    }
    if(active(am1)) {
        am_penal   = 0.5 * square((am1 - 0.437941)/0.025);                     //hard-wired
        f += am_penal; objfOut(10) = am_penal; likeOut(10) = am_penal; wgtsOut(10) = 1;
    }
    if(active(bm1)) {
        bm_penal = 0.5 * square((bm1 - 0.9487)/0.1);                          //hard-wired
        f += bm_penal; objfOut(11) = bm_penal; likeOut(11) = bm_penal; wgtsOut(11) = 1;
    }
    
    if(active(matestf)) {
        like_mat = norm2(first_difference(first_difference(matestf)));
        f += 1.0*like_mat; objfOut(12)= 1.0*like_mat; likeOut(12) = like_mat; wgtsOut(12) = 1;
    }
        
    if(active(matestm)) {
        like_mat = norm2(first_difference(first_difference(matestm)));//wts: why different weights? 1.0 vs. 0.5?
        f += 0.5*like_mat; objfOut(13)= 0.5*like_mat; likeOut(13) = like_mat; wgtsOut(13) = 0.5;
    }
    
    //jim said take the log
    sel_50m_penal.initialize();
    if(active(log_sel50_dev_mn)) {
        sel_50m_penal = like_wght_sel50*norm2(first_difference(log_sel50_dev_mn));
        f += sel_50m_penal; objfOut(14) = sel_50m_penal; likeOut(14) = norm2(first_difference(log_sel50_dev_mn)); wgtsOut(14) = like_wght_sel50;
    }
    
    // various penalties
    // =================
    double llw = 0.0;
    double red = 0.01;
    int max_number_phases = 8;
    fpen.initialize();
    if (active(log_sel50_dev_3)) { 
        int phs = log_sel50_dev_3.get_phase_start();
        llw = llw_sel50_dev_3; 
//        if (doPenRed) llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(max_number_phases-phs)))*llw;
        nextf = norm2(log_sel50_dev_3);
        fpen += llw*nextf; objfOut(38) = llw*nextf; likeOut(38) = nextf; wgtsOut(38) = llw;   //wts: need to turn this off in last phase?        
    }
    if (active(pF_DevsTCF)) { 
        int phs = pF_DevsTCF.get_phase_start();
        llw = 1.0; 
        if (doPenRed) llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(max_number_phases-phs)))*llw;
        nextf = norm2(pF_DevsTCF);
        fpen += llw*nextf; objfOut(15) = llw*nextf; likeOut(15) = nextf; wgtsOut(15) = llw;   //wts: need to turn this off in last phase?        
    }
    if(active(pF_DevsSCF)) {
        int phs = pF_DevsSCF.get_phase_start();
        llw = 0.5; 
        if (doPenRed) llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(max_number_phases-phs)))*llw;
        nextf = norm2(pF_DevsSCF);
        fpen += llw*nextf; objfOut(16) = llw*nextf; likeOut(16) = nextf; wgtsOut(16) = llw; //wts: need to turn this off in last phase? note that relative weights are hard-wired
        
    }
    if(active(pF_DevsRKF)) {
        int phs = pF_DevsRKF.get_phase_start();
        llw = 3.0; 
        if (doPenRed) llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(max_number_phases-phs)))*llw;
        nextf = norm2(pF_DevsRKF);
        fpen += llw*nextf; objfOut(17) = llw*nextf; likeOut(17) = nextf; wgtsOut(17) = llw; //wts: need to turn this off in last phase?
        
    }
    if(active(pF_DevsGTF)) {
        int phs = pF_DevsGTF.get_phase_start();
        llw = 0.5; 
        if (doPenRed) llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(max_number_phases-phs)))*llw;
        nextf = norm2(pF_DevsGTF);
        fpen += llw*nextf; objfOut(18) = llw*nextf; likeOut(18) = nextf; wgtsOut(18) = llw; //wts: need to turn this off in last phase?        
    }
    
    f += fpen;
//    cout<<"1"<<endl;    
    
    
    // LIKELIHOODS
    // ===========
    
    lkZCs.initialize();
    
    // retained (male) catch in TCF length likelihood (old and new shell together)
    for (int n=1; n <= nObsRetZCsTCF; n++) {
        yr = yrsObsRetZCsTCF_n(n);
        lkZCs(1) -= ssRetZCsTCF_sn(NEW_SHELL,n)*((obsPrNatZ_TCFR_sn(NEW_SHELL,n)+obsPrNatZ_TCFR_sn(OLD_SHELL,n))
                                                    *log(modPrNatZ_TCFR_syz(NEW_SHELL,yr)+modPrNatZ_TCFR_syz(OLD_SHELL,yr)+p_const));
    }
    //  CheckFile << "Yrs fish "<<yrsObsRetZCsTCF_n<<endl;
    //  CheckFile << "obs ret length "<<obsPrNatZ_TCFR_sn<<endl;
    //  CheckFile << "pred ret length "<<modPrNatZ_TCFR_syz<<endl;
//    cout<<"2"<<endl;    
    
    //  cout<<" ret length "<<endl;
    // total male catch in TCF length likelihood (old and new shell together) AEP???
    for (int n=1; n <= nObsZCsTCFM; n++) {
        yr = yrsObsZCsTCFM_n(n);
        lkZCs(2) -= ssTotZCsTCFM_sn(NEW_SHELL,n)*((obsPrNatZ_TCFM_snz(NEW_SHELL,n)+obsPrNatZ_TCFM_snz(OLD_SHELL,n))
                                                     *log(modPrNatZ_TCFM_syz(NEW_SHELL,yr)+modPrNatZ_TCFM_syz(OLD_SHELL,yr)+p_const));
    }
    //  cout<<" discm length "<<endl;
//    cout<<"3"<<endl;    
    
    // total female catch (discards) in TCF fishery length likelihood
    for (int n=1; n <= nObsZCsTCFF; n++) {
        yr=yrsObsZCsTCFF_n(n);
        lkZCs(3) -= ssZCsTCFF_n(n)*(obsPrNatZ_TCFF_nz(n)*log(modPrNatZ_TCFF_yz(yr)+p_const));
    }
    //      cout<<" disc f length "<<endl;
//    cout<<"4"<<endl;    
    
    // total male catch (discards) in SCF fishery length likelihood
    for (int n=1; n <= nObsZCsSCF; n++) {
        yr=yrsObsZCsSCF_n(n);
        lkZCs(4) -= ssZCsSCFM_sn(NEW_SHELL,n)*(obsPrNatZ_SCF_xnz(MALE,n)*log(modPrNatZ_SCF_xyz(MALE,yr)+p_const));
    }
//    cout<<" snow m length "<<endl;
//    CheckFile<<"+++++++++++++++++"<<endl;
//    CheckFile<<"snf males LL = "<<lkZCs(4)<<endl;
//    CheckFile<<modPrNatZ_SCF_xyz(MALE)<<endl;
//    CheckFile<<"+++++++++++++++++"<<endl;
    
    // total female catch (discards) in SCF fishery length likelihood
    for (int n=1; n <= nObsZCsSCF; n++) {
        yr=yrsObsZCsSCF_n(n);
        lkZCs(5) -= ssZCsSCFF_n(n)*(obsPrNatZ_SCF_xnz(FEMALE,n)*log(modPrNatZ_SCF_xyz(FEMALE,yr)+p_const));
    }
//    cout<<"1"<<endl;    
    
    // total male/female catch (discards) in RKF fishery length likelihood
    for (int n=1; n <= nObsZCsRKF; n++) {
        yr=yrsObsZCsRKF_n(n);
        lkZCs(6) -= ssZCsRKFM_sn(NEW_SHELL,n)*(obsPrNatZ_RKF_xnz(MALE,n)*log(modPrNatZ_RKF_xyz(MALE,yr)+p_const));
        lkZCs(7) -= ssZCsRKFF_n(n)*(obsPrNatZ_RKF_xnz(FEMALE,n)*log(modPrNatZ_RKF_xyz(FEMALE,yr)+p_const));
    }
//    cout<<" red f length "<<endl;
//    CheckFile<<"+++++++++++++++++"<<endl;
//    CheckFile<<"rkf males LL = "<<lkZCs(6)<<endl;
//    CheckFile<<modPrNatZ_RKF_xyz(MALE)<<endl;
//    CheckFile<<"+++++++++++++++++"<<endl;
    
    // total male/female catch (discards) in GTF fishery length likelihood
    for (int n=1; n <= nObsZCsGTF; n++) {
        yr=yrsObsZCsGTF(n);
        if (optPrNatZ_GTF==0){
            //old way: 
            //1. normalize observed extended size comp by total counts
            //2. weight sex-specific components of extended size comp by sex-specific ss
            for(int x=1;x<=nSXs;x++){
                lkZCs(8) -= ssZCsGTF_xn(x,n)*(obsPrNatZ_GTF_xnz(x,n)*log(modPrNatZ_GTF_xyz(x,yr)+p_const));
            }
        } else if (optPrNatZ_GTF==1){
            //new way 1: 
            //1. normalize observed extended size comp by total counts
            //2. weight extended size comp by ss summed over sexes
            for(int x=1;x<=nSXs;x++){
                lkZCs(8) -= ssObsZCsGTF_n(n)*(obsPrNatZ_GTF_xnz(x,n)*log(modPrNatZ_GTF_xyz(x,yr)+p_const));
            }
        } else if (optPrNatZ_GTF==2){
            //new way 2: 
            //1. normalize observed sex-specific size comp by count
            //2. create extended comp by weighting sex-specific parts by sex-specific ss 
            //3. weight extended size comp by ss summed over sexes
            for(int x=1;x<=nSXs;x++){
                lkZCs(8) -= ssObsZCsGTF_n(n)*(obsPrNatZ_GTF_xnz(x,n)*log(modPrNatZ_GTF_xyz(x,yr)+p_const));
            }
        } 
    }//n
//    cout<<" trawl length "<<endl;
    
    // survey likelihood
    {int x;
        for (int n=1; n <=nObsZCsSrv; n++) {
            yr=yrsObsZCsSrv_n(n);         

            x = MALE;   
            // obs(maturity, SC, x, year), pred(maturity,x, year)
            // immature new and old together
            lkZCs( 9) -= ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,x,n)*(
                             (obsPrNatZ_Srv_msxnz(IMMATURE,NEW_SHELL,x,n)+obsPrNatZ_Srv_msxnz(IMMATURE,OLD_SHELL,x,n))*
                              log(modSrvPrNatZ_NS_mxyz(IMMATURE,x,yr)+modSrvPrNatZ_OS_mxyz(IMMATURE,x,yr)+p_const)
                            );
            // this is for mature new and old shell together
            lkZCs(10) -= ssObsZCsSrv_msxn(MATURE,NEW_SHELL,x,n)*(
                            (obsPrNatZ_Srv_msxnz(MATURE,NEW_SHELL,x,n)+obsPrNatZ_Srv_msxnz(MATURE,OLD_SHELL,x,n))*
                             log(modSrvPrNatZ_NS_mxyz(MATURE,x,yr)+modSrvPrNatZ_OS_mxyz(MATURE,x,yr)+p_const)
                            );       

            x = FEMALE;   
            // obs(maturity, SC, x, year), pred(maturity,x, year)
            // immature new and old together
            lkZCs(11) -= ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,x,n)*(
                             (obsPrNatZ_Srv_msxnz(IMMATURE,NEW_SHELL,x,n)+obsPrNatZ_Srv_msxnz(IMMATURE,OLD_SHELL,x,n))*
                              log(modSrvPrNatZ_NS_mxyz(IMMATURE,x,yr)+modSrvPrNatZ_OS_mxyz(IMMATURE,x,yr)+p_const)
                            );
            // this is for mature new and old shell together
            lkZCs(12) -= ssObsZCsSrv_msxn(MATURE,NEW_SHELL,x,n)*(
                            (obsPrNatZ_Srv_msxnz(MATURE,NEW_SHELL,x,n)+obsPrNatZ_Srv_msxnz(MATURE,OLD_SHELL,x,n))*
                             log(modSrvPrNatZ_NS_mxyz(MATURE,x,yr)+modSrvPrNatZ_OS_mxyz(MATURE,x,yr)+p_const)
                            );        
        }// year loop
    }
//    cout<<"1"<<endl;    
    
    //add the offset to the likelihood   
    for (int n=1;n<=NUM_LEN_LIKE;n++) lkZCs(n) -= offset(n);
    
    // extra weight for start year length comp.
    if (current_phase() > 6) multi = 1.0; else multi = 1.0;//wts: does nothing! 
    
    nextf = like_wght( 1)*lkZCs( 1); objfOut(19) = nextf; f += nextf; likeOut(19) = lkZCs( 1); wgtsOut(19) = like_wght( 1);  // directed fishery: retained males     
    nextf = like_wght( 2)*lkZCs( 2); objfOut(20) = nextf; f += nextf; likeOut(20) = lkZCs( 2); wgtsOut(20) = like_wght( 2);  // directed fishery: total (ret+disc) males
    nextf = like_wght( 3)*lkZCs( 3); objfOut(21) = nextf; f += nextf; likeOut(21) = lkZCs( 3); wgtsOut(21) = like_wght( 3);  // directed fishery: females     
    nextf = like_wght( 4)*lkZCs( 4); objfOut(22) = nextf; f += nextf; likeOut(22) = lkZCs( 4); wgtsOut(22) = like_wght( 4);  // snow crab fishery: males
    nextf = like_wght( 5)*lkZCs( 5); objfOut(23) = nextf; f += nextf; likeOut(23) = lkZCs( 5); wgtsOut(23) = like_wght( 5);  // snow crab fishery: females
    nextf = like_wght( 6)*lkZCs( 6); objfOut(24) = nextf; f += nextf; likeOut(24) = lkZCs( 6); wgtsOut(24) = like_wght( 6);  // BBRKC fishery: males
    nextf = like_wght( 7)*lkZCs( 7); objfOut(25) = nextf; f += nextf; likeOut(25) = lkZCs( 7); wgtsOut(25) = like_wght( 7);  // BBRKC fishery: females
    nextf = like_wght( 8)*lkZCs( 8); objfOut(26) = nextf; f += nextf; likeOut(26) = lkZCs( 8); wgtsOut(26) = like_wght( 8);  // groundfish fishery: all
    nextf = like_wght( 9)*lkZCs( 9); objfOut(27) = nextf; f += nextf; likeOut(27) = lkZCs( 9); wgtsOut(27) = like_wght( 9);  // survey: immature males
    nextf = like_wght(10)*lkZCs(10); objfOut(28) = nextf; f += nextf; likeOut(28) = lkZCs(10); wgtsOut(28) = like_wght(10);  // survey: mature males
    nextf = like_wght(11)*lkZCs(11); objfOut(29) = nextf; f += nextf; likeOut(29) = lkZCs(11); wgtsOut(29) = like_wght(11);  // survey: immature females
    nextf = like_wght(12)*lkZCs(12); objfOut(30) = nextf; f += nextf; likeOut(30) = lkZCs(12); wgtsOut(30) = like_wght(12);  // survey: mature females
    
    // Fit to indices (lognormal) - AEP DIFFERENT Variance multipliers
    //weight each years estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass) 
    for(int i=1;i<=nObsSrvBio;i++) {
        cv_srv1(FEMALE,yrsObsSrvBio_n(i))  = like_wght_fbio*obsSrvCV_xn(FEMALE,i);  
        cv_srv1(  MALE,yrsObsSrvBio_n(i))  = like_wght_mbio*obsSrvCV_xn(  MALE,i);
    }
//    cout<<" survey biom "<<endl;
    //   CheckFile <<yrsObsSrvBio_n<<endl;
    //    CheckFile <<"obs surv = "<<obsSrvMatBio_xy<<endl;
    //    CheckFile <<"pred surv = "<<biom_tmp<<endl;
    // this fits mature biomass separate male and female
    zsSrvMatBio_xn.initialize();
    lkSrvMatBio_x.initialize();
    dvar_vector stdSrv(1,nObsSrvBio);
    stdSrv.initialize();
    for (int x=1;x<=nSXs;x++){
        stdSrv = sqrt(log(elem_prod(cv_srv1(x)(yrsObsSrvBio_n),cv_srv1(x)(yrsObsSrvBio_n))+1.0));
        zsSrvMatBio_xn(x) = elem_div(log(obsSrvMatBio_xy(x)(yrsObsSrvBio_n)+smlValSrv)-log(modSrvMatBio_xy(x)(yrsObsSrvBio_n)+smlValSrv),stdSrv);
        lkSrvMatBio_x(x)  = 0.5*norm2(zsSrvMatBio_xn(x));
        if (lkSrvMatBio_x(x)==0.0){
            cout<<"x = "<<x<<endl;
            cout<<"yrsObsSrvBio_n = "<<endl<<yrsObsSrvBio_n<<endl;
            cout<<"stdSrv = "<<sqrt(log(elem_prod(cv_srv1(x)(yrsObsSrvBio_n),cv_srv1(x)(yrsObsSrvBio_n))+1.0))<<endl;
            cout<<"log(obsSrvMatBio_xy(x)(yrsObsSrvBio_n) = "<<endl<<log(obsSrvMatBio_xy(x)(yrsObsSrvBio_n))<<endl;
            cout<<"log(modSrvMatBio_xy(x)(yrsObsSrvBio_n) = "<<endl<<log(modSrvMatBio_xy(x)(yrsObsSrvBio_n))<<endl;
            cout<<"stdSrv = "<<stdSrv<<endl;
            cout<<"zsSrvMatBio_xn = "<<endl<<zsSrvMatBio_xn<<endl;
            cout<<"lkSrvMatBio_x  = "<<lkSrvMatBio_x<<endl;
            exit(-1.0);
        }
    }
//    cout<<"2"<<endl;    
    
    f += sum(lkSrvMatBio_x); objfOut(31) = sum(lkSrvMatBio_x); likeOut(31) = sum(lkSrvMatBio_x); wgtsOut(31) = 1;
    
    //Likelihoods for bulk fishery catch data
    // Male retained catch in TCF
    zsRetMortBio_TCFR_y.initialize();
    lkRetMortBio_TCFR.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(2) [fixed]
        zsRetMortBio_TCFR_y(1965,endyr-1) = (obsRetCatchBio(1965,endyr-1)-predRetBioMortTCFM_y(1965,endyr-1))/sqrt(0.5);
        lkRetMortBio_TCFR                  = 0.5*norm2(zsRetMortBio_TCFR_y);
        nextf = like_wght_CatchBio*lkRetMortBio_TCFR; objfOut(32) = nextf; f += nextf; likeOut(32) = lkRetMortBio_TCFR; wgtsOut(32) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrTCFR)));
        zsRetMortBio_TCFR_y(1965,endyr-1) = (log(obsRetCatchBio(1965,endyr-1)+smlValFsh)-log(predRetBioMortTCFM_y(1965,endyr-1)+smlValFsh))/stdv;
        lkRetMortBio_TCFR = 0.5*norm2(zsRetMortBio_TCFR_y);
        nextf = lkRetMortBio_TCFR; objfOut(32) = nextf; f += nextf; likeOut(32) = lkRetMortBio_TCFR; wgtsOut(32) = 1.0;
    }
    
    // male catch MORTALITY in TCF
    zsTotMortBio_TCFM_n.initialize();
    lkTotMortBio_TCFM.initialize();
    zsDscMortBio_TCFM_n.initialize();
    lkDscMortBio_TCFM.initialize();
    if (optTCFMfit==0){
        //fit to TOTAL male mortality
        if (optFshNLLs==0){
            //normal error assumption, sd=sqrt(2) [fixed]
            zsTotMortBio_TCFM_n = (obsTotBioMortTCFM - predTotBioMortTCFM_y(yrsObsDscTCF_n))/sqrt(0.5);
            lkTotMortBio_TCFM = 0.5*norm2(zsTotMortBio_TCFM_n);
            nextf = like_wght_CatchBio*lkTotMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkTotMortBio_TCFM; wgtsOut(33) = like_wght_CatchBio;
        } else {
            //lognormal error assumption
            double stdv = sqrt(log(1.0+square(obsErrTCFD)));
            zsTotMortBio_TCFM_n = (log(obsTotBioMortTCFM+smlValFsh)-log(predTotBioMortTCFM_y(yrsObsDscTCF_n)+smlValFsh))/stdv;
            lkTotMortBio_TCFM = 0.5*norm2(zsTotMortBio_TCFM_n);
            nextf = lkTotMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkTotMortBio_TCFM; wgtsOut(33) = 1.0;
        }
    } else if (optTCFMfit==1) {
        //fit to discard male mortality
        if (optFshNLLs==0){
            //normal error assumption, sd=sqrt(2) [fixed]
            zsDscMortBio_TCFM_n = (obsDscBioMortTCF(MALE) - predDscBioMortTCFM_y(yrsObsDscTCF_n))/sqrt(0.5);
            lkDscMortBio_TCFM = 0.5*norm2(zsDscMortBio_TCFM_n);
            nextf = like_wght_CatchBio*lkDscMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkDscMortBio_TCFM; wgtsOut(33) = like_wght_CatchBio;
        } else {
            //lognormal error assumption
            double stdv = sqrt(log(1.0+square(obsErrTCFD)));
            zsDscMortBio_TCFM_n = (log(obsDscBioMortTCF(MALE)+smlValFsh)-log(predDscBioMortTCFM_y(yrsObsDscTCF_n)+smlValFsh))/stdv;
            lkDscMortBio_TCFM = 0.5*norm2(zsDscMortBio_TCFM_n);
            nextf = lkDscMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkDscMortBio_TCFM; wgtsOut(33) = 1.0;
        }
    }
    
    //  female catch in directed fishery
    zsDscMortBio_TCFF_n.initialize();
    lkDscMortBio_TCFF.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(2) [fixed]
        zsDscMortBio_TCFF_n = (obsDscBioMortTCF(FEMALE)-predDscBioMortTCFF_y(yrsObsDscTCF_n))/sqrt(0.5);
        lkDscMortBio_TCFF = 0.5*norm2(zsDscMortBio_TCFF_n);
        nextf = like_wght_CatchBio*lkDscMortBio_TCFF; objfOut(34) = nextf; f+= nextf; likeOut(34) = lkDscMortBio_TCFF; wgtsOut(34) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrTCFD)));
        zsDscMortBio_TCFF_n = (log(obsDscBioMortTCF(FEMALE)+smlValFsh)-log(predDscBioMortTCFF_y(yrsObsDscTCF_n)+smlValFsh))/stdv;
        lkDscMortBio_TCFF = 0.5*norm2(zsDscMortBio_TCFF_n);
        nextf = lkDscMortBio_TCFF; objfOut(34) = nextf; f+= nextf; likeOut(34) = lkDscMortBio_TCFF; wgtsOut(34) = 1.0;
    }
    
    //snow crab fishery
    zsDscMortBio_SCF_xn.initialize();
    lkDscMortBio_SCF_x.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(2) [fixed]
        for (int x=1;x<=nSXs;x++){
            zsDscMortBio_SCF_xn(x) = (obsDscBioMortSCF(x) - predDscBioMortSCF_xy(x)(yrsObsDscSCF))/sqrt(0.5);
            lkDscMortBio_SCF_x(x)  = 0.5*norm2(zsDscMortBio_SCF_xn(x));
        }
        nextf = like_wght_CatchBio*sum(lkDscMortBio_SCF_x); objfOut(35) = nextf; f += nextf; likeOut(35) = sum(lkDscMortBio_SCF_x); wgtsOut(35) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrSCF)));
        for (int x=1;x<=nSXs;x++){
            zsDscMortBio_SCF_xn(x) = (log(obsDscBioMortSCF(x)+smlValFsh)-log(predDscBioMortSCF_xy(x)(yrsObsDscSCF)+smlValFsh))/stdv;
            lkDscMortBio_SCF_x(x)  = 0.5*norm2(zsDscMortBio_SCF_xn(x));
        }
        nextf = sum(lkDscMortBio_SCF_x); objfOut(35) = nextf; f += nextf; likeOut(35) = sum(lkDscMortBio_SCF_x); wgtsOut(35) = 1.0;
    }
    
    //BBRKC fishery
    zsDscMortBio_RKF_xn.initialize();
    lkDscMortBio_RKF_x.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(2) [fixed]
        for (int x=1;x<=nSXs;x++){
            zsDscMortBio_RKF_xn(x) = (obsDscBioMortRKF(x) - predDscBioMortRKF_xy(x)(yrsObsDscRKF))/sqrt(0.5);
            lkDscMortBio_RKF_x(x)  = 0.5*norm2(zsDscMortBio_RKF_xn(x));
        }
        nextf = like_wght_CatchBio*sum(lkDscMortBio_RKF_x); objfOut(36) = nextf; f += nextf; likeOut(36) = sum(lkDscMortBio_RKF_x); wgtsOut(36) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrRKF)));
        for (int x=1;x<=nSXs;x++){
            zsDscMortBio_RKF_xn(x) = (log(obsDscBioMortRKF(x)+smlValFsh)-log(predDscBioMortRKF_xy(x)(yrsObsDscRKF)+smlValFsh))/stdv;
            lkDscMortBio_RKF_x(x)  = 0.5*norm2(zsDscMortBio_RKF_xn(x));
        }
        nextf = sum(lkDscMortBio_RKF_x); objfOut(36) = nextf; f += nextf; likeOut(36) = sum(lkDscMortBio_RKF_x); wgtsOut(36) = 1.0;
    }
    
    //groundfish trawl fishery
    zsDscMortBio_GTF_n.initialize();
    lkDscMortBio_GTF.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(1/2) [fixed]
        zsDscMortBio_GTF_n = (obsDscBioMortGTF - predDscBioMortGTF_y(yrsObsDscGTF))/sqrt(0.5);
        lkDscMortBio_GTF = 0.5*norm2(zsDscMortBio_GTF_n);
        nextf = like_wght_CatchBio*lkDscMortBio_GTF;  objfOut(37) = nextf; f += nextf; likeOut(37) = lkDscMortBio_GTF; wgtsOut(37) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrGTF)));
        zsDscMortBio_GTF_n = (log(obsDscBioMortGTF+smlValFsh) - (predDscBioMortGTF_y(yrsObsDscGTF)+smlValFsh))/stdv;
        lkDscMortBio_GTF = 0.5*norm2(zsDscMortBio_GTF_n);
        nextf = lkDscMortBio_GTF;  objfOut(37) = nextf; f += nextf; likeOut(37) = lkDscMortBio_GTF; wgtsOut(37) = 1.0;
    }
    // ------------------------------------------
    
    call_no += 1;
//     cout <<"Likes = "<< objfOut << endl;
//     cout <<"phase = "<< current_phase() << " call = " << call_no << " Total Like = " << f << endl;

//    cout<<"done"<<endl;
    
// ==========================================================================
// ==========================================================================
FUNCTION Misc_output
//    cout<<"Misc_output"<<endl;

    int yr;
    dvar_matrix useSelSrv(1,2,1,nZBs);
    dvar_matrix cv_srv1_nowt(1,2,styr,endyr);
    
    //   cout<<" to misc output "<<endl;
    
//    //legal size for tanner is 138mm  //vectorize? the following<--CALCULATIONS BELOW ARE WRONG!!
//    modTotLglBio_TCFM_y.initialize();
//    modTotLglNum_TCFM_y.initialize();
//    for (int yr=styr;yr<endyr;yr++) {      //(IMPORTANT CHANGE: used to be "endyr")
//        for (int j=23;j<=nZBs;j++) {
//            modTotLglNum_TCFM_y(yr) += (fmTCFM_syz(NEW_SHELL,yr,j)/(fmTCFM_syz(NEW_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_inew_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,NEW_SHELL,yr,j)) + 
//                                           (fmTCFM_syz(NEW_SHELL,yr,j)/(fmTCFM_syz(NEW_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_mnew_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,NEW_SHELL,yr,j))+ 
//                                           (fmTCFM_syz(OLD_SHELL,yr,j)/(fmTCFM_syz(OLD_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_iold_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,OLD_SHELL,yr,j))+
//                                           (fmTCFM_syz(OLD_SHELL,yr,j)/(fmTCFM_syz(OLD_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_mold_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,OLD_SHELL,yr,j));
//            modTotLglBio_TCFM_y(yr) += (fmTCFM_syz(NEW_SHELL,yr,j)/(fmTCFM_syz(NEW_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_inew_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,NEW_SHELL,yr,j)) +
//                                      + (fmTCFM_syz(NEW_SHELL,yr,j)/(fmTCFM_syz(NEW_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_mnew_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,NEW_SHELL,yr,j))+
//                                      + (fmTCFM_syz(OLD_SHELL,yr,j)/(fmTCFM_syz(OLD_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_iold_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,OLD_SHELL,yr,j))
//                                      + (fmTCFM_syz(OLD_SHELL,yr,j)/(fmTCFM_syz(OLD_SHELL,yr,j)+fmGTF_xyz(MALE,yr,j)))*natl_mold_fishtime(MALE,yr,j)*(1-S_xsyz(MALE,OLD_SHELL,yr,j))  
//                                        * wtm_z(j);
//            if (j==23) {  //note this only occurs 1st time through loop
//                // AEP???
//                modTotLglBio_TCFM_y(yr)   = modTotLglBio_TCFM_y(yr)*0.5;
//                modTotLglNum_TCFM_y(yr)= modTotLglNum_TCFM_y(yr)*0.5;
//            }
//        }//j loop
//    }//yr loop
    modTotLglBio_TCFM_y.initialize();
    modTotLglNum_TCFM_y.initialize();
    for (int yr=styr;yr<endyr;yr++) {      //(IMPORTANT CHANGE: used to be "endyr")
        for (int z=iZLegal;z<=nZBs;z++) {
            modTotLglNum_TCFM_y(yr) +=  fmTCFM_syz(NEW_SHELL,yr,z)/fmTOT_xsyz(MALE,NEW_SHELL,yr,z)*natl_inew_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,NEW_SHELL,yr,z)) + 
                                        fmTCFM_syz(NEW_SHELL,yr,z)/fmTOT_xsyz(MALE,NEW_SHELL,yr,z)*natl_mnew_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,NEW_SHELL,yr,z))+ 
                                        fmTCFM_syz(OLD_SHELL,yr,z)/fmTOT_xsyz(MALE,OLD_SHELL,yr,z)*natl_iold_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,OLD_SHELL,yr,z))+
                                        fmTCFM_syz(OLD_SHELL,yr,z)/fmTOT_xsyz(MALE,OLD_SHELL,yr,z)*natl_mold_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,OLD_SHELL,yr,z));
            modTotLglBio_TCFM_y(yr) += (fmTCFM_syz(NEW_SHELL,yr,z)/fmTOT_xsyz(MALE,NEW_SHELL,yr,z)*natl_inew_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,NEW_SHELL,yr,z)) + 
                                        fmTCFM_syz(NEW_SHELL,yr,z)/fmTOT_xsyz(MALE,NEW_SHELL,yr,z)*natl_mnew_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,NEW_SHELL,yr,z))+ 
                                        fmTCFM_syz(OLD_SHELL,yr,z)/fmTOT_xsyz(MALE,OLD_SHELL,yr,z)*natl_iold_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,OLD_SHELL,yr,z))+
                                        fmTCFM_syz(OLD_SHELL,yr,z)/fmTOT_xsyz(MALE,OLD_SHELL,yr,z)*natl_mold_fishtime(MALE,yr,z)*(1-S_xsyz(MALE,OLD_SHELL,yr,z)))  
                                        * wtm_z(z);
            if (z==iZLegal) {  //note this only occurs 1st time through loop
                // AEP???
                modTotLglBio_TCFM_y(yr) = modTotLglBio_TCFM_y(yr)*0.5;
                modTotLglNum_TCFM_y(yr) = modTotLglNum_TCFM_y(yr)*0.5;
            }
        }//z loop
    }//yr loop
//    cout<<" to large males "<<endl;
    bio_males_gt101.initialize();
    num_males_gt101.initialize();
    for (int yr=styr;yr<=endyr;yr++){
        for(int j=23;j<=nZBs;j++) {
            num_males_gt101(yr)+=  natl_inew_fishtime(MALE,yr,j) + natl_iold_fishtime(MALE,yr,j) + natl_mnew_fishtime(MALE,yr,j) + natl_mold_fishtime(MALE,yr,j);
            bio_males_gt101(yr)+= (natl_inew_fishtime(MALE,yr,j) + natl_iold_fishtime(MALE,yr,j) + natl_mnew_fishtime(MALE,yr,j) + natl_mold_fishtime(MALE,yr,j))*wtm_z(j);
            if (j==23) {  //note this only occurs 1st time through loop
                num_males_gt101(yr)=num_males_gt101(yr)*0.5;
                bio_males_gt101(yr)=bio_males_gt101(yr)*0.5;
            }
        }//j loop
    }//year loop
//    cout<<" to eff N "<<endl;
    // Effective N's                    //can these be vectorized?
    //retained males in directed fishery
    dvariable num;
    for (int i=1;i<=nObsRetZCsTCF;i++){
        yr=yrsObsRetZCsTCF_n(i);
        num = modPrNatZ_TCFR_syz(NEW_SHELL,yr)*(1-modPrNatZ_TCFR_syz(NEW_SHELL,yr))+
              modPrNatZ_TCFR_syz(OLD_SHELL,yr)*(1-modPrNatZ_TCFR_syz(OLD_SHELL,yr));
        for(int shell=1;shell<=nSCs;shell++)
            effn_fish_ret(shell,yr) = num/norm2(modPrNatZ_TCFR_syz(shell,yr)-obsPrNatZ_TCFR_sn(shell,i));
    }//i loop
    //total males in directed fishery
    for (int i=1;i<=nObsZCsTCFM;i++) {
        yr=yrsObsZCsTCFM_n(i);
        num = modPrNatZ_TCFM_syz(NEW_SHELL,yr)*(1-modPrNatZ_TCFM_syz(NEW_SHELL,yr))+
              modPrNatZ_TCFM_syz(OLD_SHELL,yr)*(1-modPrNatZ_TCFM_syz(OLD_SHELL,yr));
        for(int shell=1;shell<=nSCs;shell++)
            effn_fish_tot(shell,yr) = num/norm2(modPrNatZ_TCFM_syz(shell,yr)-obsPrNatZ_TCFM_snz(shell,i));
    }// i loop
    for (int i=1;i<=nObsZCsSrv;i++) {
        yr=yrsObsZCsSrv_n(i);
        num = modSrvPrNatZ_NS_mxyz(IMMATURE,FEMALE,yr)*(1-modSrvPrNatZ_NS_mxyz(IMMATURE,FEMALE,yr))+
              modSrvPrNatZ_NS_mxyz(IMMATURE,  MALE,yr)*(1-modSrvPrNatZ_NS_mxyz(IMMATURE,  MALE,yr))+
              modSrvPrNatZ_OS_mxyz(IMMATURE,FEMALE,yr)*(1-modSrvPrNatZ_OS_mxyz(IMMATURE,FEMALE,yr))+
              modSrvPrNatZ_OS_mxyz(IMMATURE,  MALE,yr)*(1-modSrvPrNatZ_OS_mxyz(IMMATURE,  MALE,yr));//note dot product sums
        for (int sex=1;sex<=nSXs;sex++)
            effn_srv1(IMMATURE,NEW_SHELL,sex,yr) = num/norm2(modSrvPrNatZ_NS_mxyz(IMMATURE,sex,yr)-obsPrNatZ_Srv_msxnz(IMMATURE,NEW_SHELL,sex,i));
            
        for (int sex=1;sex<=nSXs;sex++)
            effn_srv1(IMMATURE,OLD_SHELL,sex,yr) = 0.0;//no immature, old_shell animals

        num = modSrvPrNatZ_NS_mxyz(  MATURE,FEMALE,yr)*(1-modSrvPrNatZ_NS_mxyz(  MATURE,FEMALE,yr))+
              modSrvPrNatZ_NS_mxyz(  MATURE,  MALE,yr)*(1-modSrvPrNatZ_NS_mxyz(  MATURE,  MALE,yr))+
              modSrvPrNatZ_OS_mxyz(  MATURE,FEMALE,yr)*(1-modSrvPrNatZ_OS_mxyz(  MATURE,FEMALE,yr))+
              modSrvPrNatZ_OS_mxyz(  MATURE,  MALE,yr)*(1-modSrvPrNatZ_OS_mxyz(  MATURE,  MALE,yr));
        for (int sex=1;sex<=nSXs;sex++)
            effn_srv1(  MATURE,NEW_SHELL,sex,yr) = num/norm2(modSrvPrNatZ_NS_mxyz(  MATURE,sex,yr)-obsPrNatZ_Srv_msxnz(  MATURE,NEW_SHELL,sex,i));
            
        num = modSrvPrNatZ_NS_mxyz(  MATURE,FEMALE,yr)*(1-modSrvPrNatZ_NS_mxyz(  MATURE,FEMALE,yr))+
              modSrvPrNatZ_NS_mxyz(  MATURE,  MALE,yr)*(1-modSrvPrNatZ_NS_mxyz(  MATURE,  MALE,yr))+
              modSrvPrNatZ_OS_mxyz(  MATURE,FEMALE,yr)*(1-modSrvPrNatZ_OS_mxyz(  MATURE,FEMALE,yr))+
              modSrvPrNatZ_OS_mxyz(  MATURE,  MALE,yr)*(1-modSrvPrNatZ_OS_mxyz(  MATURE,  MALE,yr));
        for (int sex=1;sex<=nSXs;sex++)
            effn_srv1(  MATURE,OLD_SHELL,sex,yr) = num/norm2(modSrvPrNatZ_OS_mxyz(  MATURE,sex,yr)-obsPrNatZ_Srv_msxnz(  MATURE,OLD_SHELL,sex,i));
    } // i loop
    
    // spawning biomass and related outputs
    mspnum_new_matetime.initialize();
    mspnum_old_matetime.initialize();
    mspnum_matetime.initialize();
    mspbio_old_matetime.initialize();
    mspbio_matetime.initialize();

    fspnum_new_matetime.initialize();
    fspnum_old_matetime.initialize();
    fspnum_matetime.initialize();
    fspbio_new_matetime.initialize();
    fspbio_matetime.initialize();
//    cout<<" to sp bio matetime "<<endl;
    //spawning occurs AFTER fisheries, so cannot be evaluated in endyr
    dvar_vector tmpNS(1,nZBs);
    dvar_vector tmpOS(1,nZBs);
    for (int yr=styr;yr<endyr;yr++) {  //IMPORTANT CHANGE: was <=endyr
//        mspbio_matetime(yr) = (elem_prod(S_xsyz(  MALE,NEW_SHELL,yr)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,  MALE)),
//                                         mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,  MALE))*natlength_mnew(MALE,yr))+
//                               elem_prod(S_xsyz(  MALE,OLD_SHELL,yr)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE)),
//                                         mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE))*natlength_mold(MALE,yr))
//                              )*wtm_z;
//        fspbio_matetime(yr) = (elem_prod(S_xsyz(FEMALE,NEW_SHELL,yr)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE)),
//                                         mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE))*natlength_mnew(FEMALE,yr))+
//                               elem_prod(S_xsyz(FEMALE,OLD_SHELL,yr)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,FEMALE)),
//                                         mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,FEMALE))*natlength_mold(FEMALE,yr))
//                              )*wtf_mz(MATURE);
//        if(yr>=lyr_mort && yr<=uyr_mort && mort_switch==1) { //recalculating
//            mspbio_matetime(yr) = (elem_prod(S_xsyz(2,1,yr)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,  MALE)*mat_big(2)),mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,  MALE)*mat_big(2))*natlength_mnew(2,yr))+elem_prod(S_xsyz(2,2,yr)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2)),mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2))*natlength_mold(2,yr)))*wtm_z;
//            fspbio_matetime(yr) = (elem_prod(S_xsyz(1,1,yr)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1)),mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1))*natlength_mnew(1,yr))+elem_prod(S_xsyz(1,2,yr)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,FEMALE)*mat_big(1)),mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,FEMALE)*mat_big(1))*natlength_mold(1,yr)))*wtf_mz(MATURE);
//        }
        if (yr>=lyr_mort && yr<=uyr_mort && mort_switch==1){
//            for (int j=1;j<=nZBs;j++) {
//                mspnum_old_matetime(yr)  += S_xsyz(2,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2))*natlength_mold(2,yr,j);
//                mspnum_matetime(yr)      += S_xsyz(2,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,  MALE)*mat_big(2))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,  MALE)*mat_big(2))*natlength_mnew(2,yr,j) + 
//                                            S_xsyz(2,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2))*natlength_mold(2,yr,j);
//                mspbio_old_matetime(yr)  += (S_xsyz(2,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(2))*natlength_mold(2,yr,j))*wtm_z(j);
//                fspnum_new_matetime(yr)  += (S_xsyz(1,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1))*natlength_mnew(1,yr,j));
//                fspbio_new_matetime(yr)  += (S_xsyz(1,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1))*natlength_mnew(1,yr,j))*wtf_mz(MATURE,j);
//                fspnum_matetime(yr)     += (S_xsyz(1,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(1))*natlength_mnew(1,yr,j)+
//                                             S_xsyz(1,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,FEMALE)*mat_big(1))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,FEMALE)*mat_big(1))*natlength_mold(1,yr,j));
//            }//j loop
            tmpNS = elem_prod(S_xsyz(  MALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,  MALE)*mat_big(  MALE))*natlength_mnew(  MALE,yr));
            tmpOS = elem_prod(S_xsyz(  MALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,  MALE)*mat_big(  MALE))*natlength_mold(  MALE,yr));
            mspnum_new_matetime(yr) = sum(tmpNS);
            mspnum_old_matetime(yr) = sum(tmpOS);
            mspnum_matetime(yr)     = mspnum_new_matetime(yr)+mspnum_old_matetime(yr);
            mspbio_old_matetime(yr) = tmpOS*wtm_z;         //dot product
            mspbio_matetime(yr)     = (tmpNS+tmpOS)*wtm_z; //dot product
            
            tmpNS = elem_prod(S_xsyz(FEMALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,FEMALE)*mat_big(FEMALE))*natlength_mnew(FEMALE,yr));
            tmpOS = elem_prod(S_xsyz(FEMALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,FEMALE)*mat_big(FEMALE))*natlength_mold(FEMALE,yr));
            fspnum_new_matetime(yr) = sum(tmpNS);
            fspnum_old_matetime(yr) = sum(tmpOS);
            fspnum_matetime(yr)     = fspnum_new_matetime(yr)+fspnum_old_matetime(yr);
            fspbio_new_matetime(yr) = tmpNS*wtf_mz(MATURE);     //dot product
            fspbio_matetime(yr) = (tmpNS+tmpOS)*wtf_mz(MATURE); //dot product
        } else {
//            for(int j=1;j<=nZBs;j++) {
//                mspnum_old_matetime(yr) += S_xsyz(2,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE))*natlength_mold(2,yr,j);
//                mspnum_matetime(yr)      += S_xsyz(2,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,  MALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,  MALE))*natlength_mnew(2,yr,j) + 
//                                            S_xsyz(2,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE))*natlength_mold(2,yr,j);
//                mspbio_old_matetime(yr)  += (S_xsyz(2,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,  MALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,  MALE))*natlength_mold(2,yr,j))*wtm_z(j);
//                fspnum_new_matetime(yr)  += (S_xsyz(1,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE))*natlength_mnew(1,yr,j));
//                fspbio_new_matetime(yr)  += (S_xsyz(1,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE))*natlength_mnew(1,yr,j))*wtf_mz(MATURE,j);
//                fspnum_matetime(yr)     += (S_xsyz(1,1,yr,j)*mfexp(-(spmo)*M_msx(MATURE,NEW_SHELL,FEMALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,FEMALE))*natlength_mnew(1,yr,j)+
//                                             S_xsyz(1,2,yr,j)*mfexp(-(spmo)*M_msx(MATURE,OLD_SHELL,FEMALE))*mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,FEMALE))*natlength_mold(1,yr,j));
//            }//j loop
            tmpNS = elem_prod(S_xsyz(  MALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,  MALE))*natlength_mnew(  MALE,yr));
            tmpOS = elem_prod(S_xsyz(  MALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,  MALE))*natlength_mold(  MALE,yr));
            mspnum_new_matetime(yr) = sum(tmpNS);
            mspnum_old_matetime(yr) = sum(tmpOS);
            mspnum_matetime(yr)     = mspnum_new_matetime(yr)+mspnum_old_matetime(yr);
            mspbio_old_matetime(yr) = tmpOS*wtm_z;         //dot product
            mspbio_matetime(yr)     = (tmpNS+tmpOS)*wtm_z; //dot product
            
            tmpNS = elem_prod(S_xsyz(FEMALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,FEMALE))*natlength_mnew(FEMALE,yr));
            tmpOS = elem_prod(S_xsyz(FEMALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,FEMALE))*natlength_mold(FEMALE,yr));
            fspnum_new_matetime(yr) = sum(tmpNS);
            fspnum_old_matetime(yr) = sum(tmpOS);
            fspnum_matetime(yr)     = fspnum_new_matetime(yr)+fspnum_old_matetime(yr);
            fspbio_new_matetime(yr) = tmpNS*wtf_mz(MATURE);        //dot product
            fspbio_matetime(yr)     = (tmpNS+tmpOS)*wtf_mz(MATURE);//dot product
        }
        // effective sp numbers
        emspbio_matetime(yr) = mspbio_old_matetime(yr);
        
        // for male old shell mating only (AEP ERROR IN OLD CODE HAS >=)
        efspbio_matetime(yr) = fspbio_matetime(yr);
        if (mspnum_old_matetime(yr) < (fspnum_matetime(yr)/mate_ratio)) efspbio_matetime(yr) = fspbio_matetime(yr)*((mspnum_old_matetime(yr)*mate_ratio)/fspnum_matetime(yr));
        
        // effective sp numbers for new shell females
        efspbio_new_matetime(yr) = fspbio_new_matetime(yr);
        if (mspnum_old_matetime(yr) < fspnum_new_matetime(yr)/mate_ratio) efspbio_new_matetime(yr) = fspbio_new_matetime(yr)*((mspnum_old_matetime(yr)*mate_ratio)/fspnum_new_matetime(yr));
    }//year loop
    
    //spawning biomass prior to fisheries CAN be evaluated in endyr
    for (int yr=styr;yr<= endyr;yr++) {
        mspbio_fishtime(yr) = (natl_mnew_fishtime(  MALE,yr)+natl_mold_fishtime(  MALE,yr))*wtm_z;
        fspbio_fishtime(yr) = (natl_mnew_fishtime(FEMALE,yr)+natl_mold_fishtime(FEMALE,yr))*wtf_mz(MATURE);
    }//year loop
//    cout<<" sex ratio "<<endl;
    // Sex ratio
    for (int yr=styr;yr<=endyr;yr++){
        if ((sum(natlength(FEMALE,yr))+sum(natlength(MALE,yr)))<0.01) { 
            predpop_sexr(yr)=0.0;
        } else {
            predpop_sexr(yr)=sum(natlength(FEMALE,yr))/(sum(natlength(FEMALE,yr))+sum(natlength(  MALE,yr)));
        }
    }
//    cout<<" to age - struct "<<endl;
    
    // Age-structure
    natlength_mold_age.initialize();
    
    // initialize
    dvariable tmpi = 1.0;
    for(int j=1;j<=(nages-3);j++) tmpi += mfexp(-j*M_msx(MATURE,OLD_SHELL,FEMALE));
    natlength_mold_age(1,styr,1) = natlength_mold(1,styr)/(tmpi+(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,FEMALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,FEMALE)))));
    
    for(int j=1;j<=(nages-2);j++) natlength_mold_age(1,styr,j+1) = natlength_mold_age(1,styr,1)*mfexp(-j*M_msx(MATURE,OLD_SHELL,FEMALE));
    natlength_mold_age(1,styr,nages) = natlength_mold_age(1,styr,1)*(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,FEMALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,FEMALE))));
    
    tmpi = 1.0;
    for(int j=1;j<=(nages-3);j++) tmpi += mfexp(-j*M_msx(MATURE,OLD_SHELL,  MALE));
    natlength_mold_age(2,styr,1) = natlength_mold(2,styr)/(tmpi+(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,  MALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,  MALE)))));
    
    for(int j=1;j<=(nages-2);j++) natlength_mold_age(2,styr,j+1) = natlength_mold_age(2,styr,1)*mfexp(-j*M_msx(MATURE,OLD_SHELL,  MALE));
    natlength_mold_age(2,styr,nages) = natlength_mold_age(2,styr,1)*(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,  MALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,  MALE))));
    
    //numbers at length from styr to endyr
    for (int sex=1;sex<=2;sex++) {
        for (int i=styr;i< endyr;i++) {
            // for numbers by length and age assumes no molting after maturity
            natlength_mold_age(sex,i+1,1) = mfexp(-(1-mdptFshs_y(i))*M_msx(MATURE,NEW_SHELL, sex)) * 
                                       elem_prod(S_xsyz(sex,1,i),mfexp(-mdptFshs_y(i)*M_msx(MATURE,NEW_SHELL, sex))*natlength_mnew(sex,i));
            for(int j=1;j<=(nages-1);j++){
                natlength_mold_age(sex,i+1,j+1) = (mfexp(-(1-mdptFshs_y(i))*M_msx(MATURE,OLD_SHELL, sex)) * 
                                       elem_prod(S_xsyz(sex,2,i),mfexp(-mdptFshs_y(i)*M_msx(MATURE,OLD_SHELL, sex))*natlength_mold_age(sex,i,j)));
            }
            natlength_mold_age(sex,i+1,nages) += (mfexp(-(1-mdptFshs_y(i))*M_msx(MATURE,OLD_SHELL, sex)) * 
                                   elem_prod(S_xsyz(sex,2,i),mfexp(-mdptFshs_y(i)*M_msx(MATURE,OLD_SHELL, sex))*natlength_mold_age(sex,i,nages)));
        }
    }
//    cout<<" to legal males "<<endl;
    // Legal males
    modPopNum_y.initialize();
    bioLegalMales_y.initialize();
    modSrvNumLegal_sy.initialize();
    modSrvNum_xyz.initialize();
    modSrvBio_xy.initialize();
    for (int yr=styr;yr<=endyr;yr++) {
        // Selection pattern
        if (yr<1982) useSelSrv = selSrv2;
        if (yr>1981 && yr<1988) useSelSrv = selSrv2a;
        if (yr>1987) useSelSrv = selSrv3;
        
//        // legal is >102mm take half the numbers in the 100-105 bin
//        bioLegalMales_y(yr) = numLegalMales_y(yr)*wtm_z(23);
//        modSrvNumLegal_sy(yr) = 0.5*natlength_new(2,yr,23)*useSelSrv(2,23);
//        legal_srv_males_o(yr) = 0.5*natlength_old(2,yr,23)*useSelSrv(2,23);
//        modSrvBioLegal_y(yr) = modSrvNumLegal_y(yr)*wtm_z(23);
//        for(int j=24;j<=nZBs;j++) {
//            bioLegalMales_y(yr) += natlength(2,yr,j)*wtm_z(j);
//            modSrvNumLegal_sy(yr) += natlength_new(2,yr,j)*useSelSrv(2,j);
//            legal_srv_males_o(yr) += natlength_old(2,yr,j)*useSelSrv(2,j);
//            modSrvBioLegal_y(yr) += natlength(2,yr,j)*useSelSrv(2,j)*wtm_z(j);
//        }
        // legal male size is based on zLegal
        bioLegalMales_y(yr)             = natlength(MALE,yr)*wtm_z;
        modSrvNumLegal_sy(NEW_SHELL,yr) = natlength_new(MALE,yr)*useSelSrv(MALE);
        modSrvNumLegal_sy(OLD_SHELL,yr) = natlength_old(MALE,yr)*useSelSrv(MALE);
        modSrvBioLegal_y(yr)            = elem_prod(natlength(MALE,yr),useSelSrv(MALE))*wtm_z;
        
         for (int x=1;x<=nSXs;x++){
            // survey numbers
            modSrvMatNum_xsy(x,NEW_SHELL,yr) = multQ*natlength_mnew(x,yr)*useSelSrv(x);
            modSrvMatNum_xsy(x,OLD_SHELL,yr) = multQ*natlength_mold(x,yr)*useSelSrv(x);
        
            // total survey summaries
            if(x==FEMALE) {
                modSrvBio_xy(x,yr) = multQ*((natlength_inew(x,yr)                     )*elem_prod(useSelSrv(x),wtf_mz(IMMATURE))+
                                            (natlength_mnew(x,yr)+natlength_mold(x,yr))*elem_prod(useSelSrv(x),wtf_mz(  MATURE)));
            } else {
                modSrvBio_xy(x,yr) = multQ*(natlength(x,yr)*elem_prod(useSelSrv(x),wtm_z));
            }
            modSrvNum_xyz(x,yr) = multQ*elem_prod(natlength(x,yr),useSelSrv(x));
            modPopNum_y(yr)   += sum(natlength(x,yr));
        }//x
    }//yr
//    cout<<" to surv like "<<endl;
    // Survey likelihood (by year)
    lkSrvZCs_msx.initialize();
    for(int x=1;x<=2;x++) {
        for (int i=1; i <=nObsZCsSrv; i++) {
            yr=yrsObsZCsSrv_n(i);            
            for (int z=1; z<=nZBs; z++) {
                // immature new and old together in likelihood indices are (mat,shell,x,year,length)
                lkSrvZCs_msx(1,1,x) -= ssObsZCsSrv_msxn(1,1,x,i)*(1e-9+obsPrNatZ_Srv_msxnz(1,1,x,i,z)+obsPrNatZ_Srv_msxnz(1,2,x,i,z))*log(modSrvPrNatZ_NS_mxyz(1,x,yr,z)+modSrvPrNatZ_OS_mxyz(1,x,yr,z)+1e-9);
                lkSrvZCs_msx(1,2,x) = 0.0;
                // mature
                lkSrvZCs_msx(2,1,x) -= ssObsZCsSrv_msxn(2,1,x,i)*(1e-9+obsPrNatZ_Srv_msxnz(2,1,x,i,z))*log(modSrvPrNatZ_NS_mxyz(2,x,yr,z)+1e-9);
                lkSrvZCs_msx(2,2,x) -= ssObsZCsSrv_msxn(2,2,x,i)*(1e-9+obsPrNatZ_Srv_msxnz(2,2,x,i,z))*log(modSrvPrNatZ_OS_mxyz(2,x,yr,z)+1e-9);            
            }  //z loop     
        } // year loop
    } //x loop
    
//    //weight each years estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass) 
//    for(int i=1;i<=nObsSrvBio;i++) {
//        cv_srv1_nowt(1,yrsObsSrvBio_n(i)) = obsSrvCV_xn(1,i);
//        cv_srv1_nowt(2,yrsObsSrvBio_n(i)) = obsSrvCV_xn(2,i);
//    }
//    
//    // Combined likelihood  --  different small value used in actual likelihood calculation
//    surv_like_nowt.initialize();
//    for(int sex=1;sex<=2;sex++) {
//        surv_like_nowt += norm2(elem_div( log(obsSrvBio_xy(sex)(yrsObsSrvBio_n)+.01)-log(biom_tmp(sex)(yrsObsSrvBio_n)+.01),
//                                          sqrt(2)*sqrt(log(elem_prod(cv_srv1_nowt(sex)(yrsObsSrvBio_n),cv_srv1_nowt(sex)(yrsObsSrvBio_n))+1.0))
//                                         ));
//    }

//    cout<<"before sd"<<endl;
    if(sd_phase()){
        sdrDepletion = pred_bio(endyr) / pred_bio(styr);
        sdrSpBioF    = fspbio(sdrSpBioF.indexmin(),sdrSpBioF.indexmax());
        sdrSpBioM    = mspbio(sdrSpBioM.indexmin(),sdrSpBioM.indexmax());
        //cout<<"0"<<endl;
        sdrLegalMales = numLegalMales_y(sdrLegalMales.indexmin(),sdrLegalMales.indexmax());
        //cout<<"1"<<endl;
        sdrRecEarly = rec_y(styr,mnYrRecCurr-1);
        sdrRecF     = rec_y(mnYrRecCurr,endyr);//was endyr-1
        sdrRecM     = rec_y(mnYrRecCurr,endyr);//was endyr-1
        //cout<<"2"<<endl;
        sdrMMB      = mspbio_matetime(styr+1,endyr-1);
        //cout<<"3"<<endl;
        for (int i=(styr+1);i<=(endyr-reclag);i++) {
            sdrRec(i)      = rec_y(i+reclag-1);
            sdrLnRec(i)    = log(rec_y(i+reclag-1));
            sdrLnRecMMB(i) = log(rec_y(i+reclag-1))-log(sdrMMB(i));
        }
    }
//    cout<<"done"<<endl;

  
// ==========================================================================
FUNCTION void writeMyReportToR(ostream& report)
    cout<<"Starting myReportToR"<<endl;
    


// ==========================================================================
FUNCTION void writeReport(ostream& report)
    CheckFile<<"start writeReport"<<endl;
    int ii,i,k,j;
    dvar_vector preds_sexr(styr,endyr);
    dvar_matrix tmpo(1,2,styr,endyr);
    dvar_matrix tmpp(1,2,styr,endyr);
    dvar_vector obs_tmp(styr,endyr);
    dvariable ghl,ghl_number;
    dvariable hrate;
    dvar_vector totcatch(styr, endyr-1);  
    Misc_output();
    tmpp1=0.0;
    tmpp2=0.0;
    tmpp3=0.0;
    tmpp4=0.0;
  
    writeLikelihoodComponents(report,0);
      
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"Natural mortality"<<endl;
  report<<"----------------------------------------------------------------------------"<<endl;
  dvar_vector mrt(styr,endyr-1);
  mrt = M_msx(IMMATURE,NEW_SHELL,FEMALE);
  report<<"Natural mortality, immature females"<<endl<<mrt<<endl;
  mrt = M_msx(IMMATURE,NEW_SHELL,  MALE);
  report<<"Natural mortality, immature males"<<endl<<mrt<<endl;
  mrt = M_msx(MATURE,NEW_SHELL,FEMALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= mat_big(FEMALE);
  report<<"Natural mortality, mature, new shell females"<<endl<<mrt<<endl;
  mrt = M_msx(MATURE,OLD_SHELL,FEMALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= mat_big(FEMALE);
  report<<"Natural mortality, mature, old shell females"<<endl<<mrt<<endl;
  mrt = M_msx(MATURE,NEW_SHELL,  MALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= mat_big(  MALE);
  report<<"Natural mortality, mature, new shell males"<<endl<<mrt<<endl;
  mrt = M_msx(MATURE,OLD_SHELL,  MALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= mat_big(  MALE);
  report<<"Natural mortality, mature, old shell males"<<endl<<mrt<<endl;
      
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"Population numbers"<<endl;
  report<<"----------------------------------------------------------------------------"<<endl;
    for (i=styr;i<=endyr;i++)
    if((modSrvNum_xy(1,i)+modSrvNum_xy(2,i))<0.01) {
        preds_sexr(i)=0.0;
    } else {
        preds_sexr(i)=modSrvNum_xy(1,i)/(modSrvNum_xy(1,i)+modSrvNum_xy(2,i));
    }
    report << "Estimated numbers of immature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report <<  i<<" "<<natlength_inew(1,i) << endl;
    report << "Estimated numbers of immature old shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report <<  i<<" "<<natlength_iold(1,i) << endl;
    report << "Estimated numbers of mature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report <<  i<<" "<<natlength_mnew(1,i) << endl;
    report << "Estimated numbers of mature old shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++)report <<  i<<" "<<natlength_mold(1,i) << endl;
    
    report << "Estimated numbers of immature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report << i<<" "<<natlength_inew(2,i) << endl;
    report << "Estimated numbers of immature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report << i<<" "<<natlength_iold(2,i) << endl;
    report << "Estimated numbers of mature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report << i<<" "<<natlength_mnew(2,i) << endl;
    report << "Estimated numbers of mature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report << i<<" "<<natlength_mold(2,i) << endl;
    
    report << "Observed numbers of immature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(1,1,1,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
    report << "Observed numbers of mature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,1,1,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
    report << "Observed numbers of mature old shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,2,1,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
    report << "Observed numbers of immature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(1,1,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
    report << "Observed numbers of immature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(1,2,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
    report << "Observed numbers of mature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,1,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
    report << "Observed numbers of mature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,2,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
    report << "Observed Survey Numbers by length females:  'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" " << obsSrvNum_xyz(1,yrsObsZCsSrv_n(i)) << endl;
    report << "Observed Survey Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" " << obsSrvNum_xyz(2,yrsObsZCsSrv_n(i))<< endl;
    
    report << "Predicted Survey Numbers by length females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "  << modSrvNum_xyz(1,yrsObsZCsSrv_n(i)) << endl;
    report << "Predicted Survey Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for (i=1; i <= nObsZCsSrv; i++) report<<yrsObsZCsSrv_n(i)<<" "  << modSrvNum_xyz(2,yrsObsZCsSrv_n(i)) << endl;
    report << "Predicted pop Numbers by length females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report<<i<<" "<< natlength(1,i)<< endl;
    report << "Predicted pop Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
    for(i=styr;i<=endyr;i++) report<<i<<" "<< natlength(2,i)<< endl;
    
    //actual years for obs survey male are 1969,1970,1972-2009
    report<<"observed number of males greater than 101 mm: seq(1965,"<<endyr<<")"<<endl;
    report<<obsSrvLegalNum_n<<endl;
    report<<"observed biomass of males greater than 101 mm: seq(1965,"<<endyr<<")"<<endl;
    report<<obsSrvLegalBio_n<<endl;
    report<<"pop estimate numbers of males >101: seq("<<styr<<","<<endyr<<")"<<endl;
    report<<numLegalMales_y<<endl;
    report<<"estimated population biomass of males > 101: seq("<<styr<<","<<endyr<<") "<<endl;
    report<<bioLegalMales_y<<endl;
    report<<"estimated survey numbers of males > 101: seq("<<styr<<","<<endyr<<") "<<endl;
    report<<modSrvNumLegal_y<<endl;
    report<<"estimated survey biomass of males > 101: seq("<<styr<<","<<endyr<<") "<<endl;
    report<<modSrvBioLegal_y<<endl;
    report << "Observed survey biomass: seq(1974,"<<endyr<<")"<<endl;
    report << obsSrvBio_y(1974,endyr)<<endl;
    report << "predicted survey biomass: seq("<<styr<<","<<endyr<<")"<<endl;
    report << modSrvBio_xy(1)+modSrvBio_xy(2)<<endl;
  
    //survey numbers
    for(k=1;k<=2;k++){
        for(i=styr;i<=endyr;i++)  {
            tmpo(k,i)=sum(obsSrvNum_xyz(k,i));
            tmpp(k,i)=sum(modSrvNum_xyz(k,i));
        }
    }
    report << "Observed survey numbers female: seq("<<styr<<","<<endyr<<")"<<endl;
    report << tmpo(1)<<endl;
    report << "Observed survey numbers male: seq("<<styr<<","<<endyr<<")"<<endl;
    report << tmpo(2)<<endl;
    report << "predicted survey numbers female: seq("<<styr<<","<<endyr<<")"<<endl;
    report << tmpp(1)<<endl;
    report << "predicted survey numbers male: seq("<<styr<<","<<endyr<<")"<<endl;
    report << tmpp(2)<<endl;
    report << "Observed survey female spawning biomass: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvMatBio_xy(1)<<endl;
    report << "Observed survey male spawning biomass: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvMatBio_xy(2)<<endl;
    report << "Observed survey female new spawning numbers: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvMatNum_sxy(1,1)<<endl;
    report << "Observed survey female old spawning numbers: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvMatNum_sxy(2,1)<<endl;
    report << "Observed survey male new spawning numbers: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvMatNum_sxy(1,2)<<endl;
    report << "Observed survey male old spawning numbers: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvMatNum_sxy(2,2)<<endl;
    report << "Observed survey female biomass: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvBio_xy(1)<<endl;
    report << "Observed survey male biomass: seq("<<styr<<","<<endyr<<")"<<endl;
    report << obsSrvBio_xy(2)<<endl;
    report << "natural mortality immature females, males: 'FemM','MaleM'" << endl;
    report << M_msx(IMMATURE,NEW_SHELL) << endl;
    report << "natural mortality mature new shell females, males: 'FemMm','MaleMm'" << endl;
    report << M_msx(MATURE,NEW_SHELL) << endl;
    report << "natural mortality mature old shell females, males: 'FemMmo','MaleMmo'" << endl;
    report << M_msx(MATURE,OLD_SHELL) << endl;
    report << "Predicted Biomass: seq("<<styr<<","<<endyr<<")" << endl;
    report << pred_bio << endl;
    report << "Predicted total population numbers: seq("<<styr<<","<<endyr<<") "<<endl;
    report <<modPopNum_y<<endl;
    report << "Female Spawning Biomass: seq("<<styr<<","<<endyr<<") " << endl;
    report << fspbio << endl;
    report << "Male Spawning Biomass (1000's t): seq("<<styr<<","<<endyr<<") " << endl;
    report << mspbio << endl;
    report << "Total Spawning Biomass: seq("<<styr<<","<<endyr<<") " << endl;
    report << fspbio+mspbio << endl;
    report << "Female Spawning Biomass at fish time: seq("<<styr<<","<<endyr<<") " << endl;
    report << fspbio_fishtime << endl;
    report << "Male Spawning Biomass at fish time: seq("<<styr<<","<<endyr<<") " << endl;
    report << mspbio_fishtime << endl;
    report << "Total Spawning Biomass at fish time: seq("<<styr<<","<<endyr<<") " << endl;
    report << fspbio_fishtime+mspbio_fishtime << endl;
    report << "Mating time Female Spawning Biomass: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << fspbio_matetime << endl;
    report << "Mating time Male Spawning Biomass: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << mspbio_matetime << endl;
    report << "Mating time Male old shell Spawning Biomass: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << mspbio_old_matetime << endl;
    report << "Mating time Female new shell Spawning Biomass: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << fspbio_new_matetime << endl;
    report << "Mating time Total Spawning Biomass : seq("<<styr<<","<<endyr-1<<") " << endl;
    report << fspbio_matetime+mspbio_matetime << endl;
    report << "Mating time effective Female Spawning Biomass: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << efspbio_matetime << endl;
    report << "Mating time effective Male Spawning Biomass(old shell only): seq("<<styr<<","<<endyr-1<<") " << endl;
    report << emspbio_matetime << endl;
    report << "Mating time Total effective Spawning Biomass: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << efspbio_matetime+emspbio_matetime << endl;
    report << "Mating time male Spawning numbers: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << mspnum_matetime << endl;
    report << "Mating time Female Spawning numbers: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << fspnum_matetime << endl;
    report << "Mating time Male Spawning numbers(old shell only): seq("<<styr<<","<<endyr-1<<") " << endl;
    report << mspnum_old_matetime << endl;
    //  report << "ratio Mating time Female Spawning numbers to male old shell mature numbers : seq("<<styr<<","<<endyr<<") " << endl;
    //  report << elem_div(fspnum_matetime,mspnum_old_matetime) << endl;
    report << "Mating time effective Female new shell Spawning biomass: seq("<<styr<<","<<endyr-1<<") " << endl;
    report <<efspbio_new_matetime << endl;
    report << "Mating time Female new shell Spawning numbers: seq("<<styr<<","<<endyr-1<<") " << endl;
    report << fspnum_new_matetime << endl;
    //  report << "ratio Mating time Female new shell Spawning numbers to male old shell mature numbers : seq("<<styr<<","<<endyr<<") " << endl;
    //            for(i=styr;i<=endyr;i++){if(mspnum_old_matetime(i)<0.001) mspnum_old_matetime(i)=1.0; }
    //  report << elem_div(fspnum_new_matetime,mspnum_old_matetime) << endl;
    report << "Predicted Female survey Biomass: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvBio_xy(1) << endl;
    report << "Predicted Male survey Biomass: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvBio_xy(2)<< endl;
    report << "Predicted Female survey mature Biomass: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvMatBio_xy(FEMALE) << endl;
    report << "Predicted Male survey mature Biomass: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvMatBio_xy(  MALE)<< endl;
    report << "Predicted total survey mature Biomass: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvMatBio_xy(FEMALE)+modSrvMatBio_xy(  MALE)<< endl;
    report << "Predicted Female survey new mature numbers: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvMatNum_xsy(FEMALE,NEW_SHELL) << endl;
    report << "Predicted Female survey old mature numbers: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvMatNum_xsy(FEMALE,OLD_SHELL) << endl;
    report << "Predicted Male survey new mature numbers: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvMatNum_xsy(  MALE,NEW_SHELL)<< endl;
    report << "Predicted Male survey old mature numbers: seq("<<styr<<","<<endyr<<") " << endl;
    report << modSrvMatNum_xsy(  MALE,OLD_SHELL)<< endl;

  report << "Observed Prop fishery ret new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for (i=1; i<=nObsRetZCsTCF; i++) report << yrsObsRetZCsTCF_n(i) << " " << obsPrNatZ_TCFR_sn(1,i)<< endl;
  report << "Predicted length prop fishery ret new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsRetZCsTCF; i++) 
   {
    ii=yrsObsRetZCsTCF_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_TCFR_syz(1,ii)  << endl;
   }
  report << "Observed Prop fishery ret old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for (i=1; i<=nObsRetZCsTCF; i++) report << yrsObsRetZCsTCF_n(i) << " " << obsPrNatZ_TCFR_sn(2,i)<< endl;
  report << "Predicted length prop fishery ret old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsRetZCsTCF; i++)
   {
    ii=yrsObsRetZCsTCF_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_TCFR_syz(2,ii)  << endl;
   }

  report << "Observed Prop fishery total new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for (i=1; i<=nObsZCsTCFM; i++) report << yrsObsZCsTCFM_n(i) << " " << obsPrNatZ_TCFM_snz(1,i) << endl;
  report << "Predicted length prop fishery total new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsTCFM; i++)
   {
    ii=yrsObsZCsTCFM_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_TCFM_syz(1,ii)  << endl;
   }
  report << "Observed Prop fishery total old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for (i=1; i<=nObsZCsTCFM; i++) report << yrsObsZCsTCFM_n(i) << " " << obsPrNatZ_TCFM_snz(2,i) << endl;
  report << "Predicted length prop fishery total old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsTCFM; i++)
   {
    ii=yrsObsZCsTCFM_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_TCFM_syz(2,ii)  << endl;
   }
  report << "Observed Prop fishery discard new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for (i=1; i<=nObsZCsTCFM; i++) report << yrsObsZCsTCFM_n(i) << " " << obs_p_fish_discm(NEW_SHELL,i) << endl;
  report << "Observed Prop fishery discard old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for (i=1; i<=nObsZCsTCFM; i++) report << yrsObsZCsTCFM_n(i) << " " << obs_p_fish_discm(OLD_SHELL,i)<< endl;

  report << "Observed length prop fishery discard all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsTCFF; i++) report <<  yrsObsZCsTCFF_n(i)  <<  " "  <<  obsPrNatZ_TCFF_nz(i)  << endl;
  report << "Predicted length prop fishery discard all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsTCFF; i++)
   {
    ii=yrsObsZCsTCFF_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_TCFF_yz(ii)  << endl;
   }
  report << "Observed length prop snow fishery females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSCF; i++)
   {
    report <<  yrsObsZCsSCF_n(i)  <<  " "  <<  obsPrNatZ_SCF_xnz(1,i)  << endl;
   }
  report << "Predicted length prop snow fishery females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSCF; i++)
   {
    ii=yrsObsZCsSCF_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_SCF_xyz(1,ii)  << endl;
   }
  report << "Observed length prop snow fishery males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSCF; i++)
   {
    report <<  yrsObsZCsSCF_n(i)  <<  " "  <<  obsPrNatZ_SCF_xnz(2,i)  << endl;
   }
  report << "Predicted length prop snow fishery males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSCF; i++)
   {
    ii=yrsObsZCsSCF_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_SCF_xyz(2,ii)  << endl;
   }
  report << "Observed length prop redk fishery females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsRKF; i++)
   {
    report <<  yrsObsZCsRKF_n(i)  <<  " "  <<  obsPrNatZ_RKF_xnz(1,i)  << endl;
   }
  report << "Predicted length prop redk fishery females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsRKF; i++)
   {
    ii=yrsObsZCsRKF_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_RKF_xyz(1,ii)  << endl;
   }
  report << "Observed length prop redk fishery males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsRKF; i++)
   {
    report <<  yrsObsZCsRKF_n(i)  <<  " "  <<  obsPrNatZ_RKF_xnz(2,i)  << endl;
   }
  report << "Predicted length prop redk fishery males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsRKF; i++)
   {
    ii=yrsObsZCsRKF_n(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_RKF_xyz(2,ii)  << endl;
   }

  report << "Predicted length prop trawl females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsGTF; i++)
   {
    ii=yrsObsZCsGTF(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_GTF_xyz(1,ii)  << endl;
   }
  report << "Observed length prop trawl females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsGTF; i++) report <<  yrsObsZCsGTF(i)  <<  " "  <<  obsPrNatZ_GTF_xnz(1,i)  << endl;
  report << "Predicted length prop trawl males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsGTF; i++)
   {
    ii=yrsObsZCsGTF(i);  
    report <<  ii  <<  " "  <<  modPrNatZ_GTF_xyz(2,ii)  << endl;
   }
  report << "Observed length prop trawl males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsGTF; i++) report <<  yrsObsZCsGTF(i)  <<  " "  <<  obsPrNatZ_GTF_xnz(2,i)  << endl;

  report << "Observed Length Prop survey immature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);
    report << ii <<" " <<obsPrNatZ_Srv_msxnz(1,1,1,i) << endl;
   }
  report << "Predicted length prop survey immature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_NS_mxyz(1,1,ii) << endl;
   }
  report << "Observed Length Prop survey immature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);
    report << ii <<" " <<obsPrNatZ_Srv_msxnz(1,2,1,i) << endl;
   }
  report << "Predicted length prop survey immature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_OS_mxyz(1,1,ii) << endl;
   }
 
  report << "Observed Length Prop survey immature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++) report << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(1,1,2,i) << endl;
  report << "Predicted length prop survey immature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_NS_mxyz(1,2,ii) << endl;
   }
  report << "Observed Length Prop survey immature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++) report << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(1,2,2,i) << endl;
  report << "Predicted length prop survey immature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
  {
   ii=yrsObsZCsSrv_n(i);  
   report << ii << " " << modSrvPrNatZ_OS_mxyz(1,2,ii) << endl;
  }
  report << "Observed Length Prop survey mature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);
    report << ii <<" " <<obsPrNatZ_Srv_msxnz(2,1,1,i) << endl;
   }
  report << "Predicted length prop survey mature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_NS_mxyz(2,1,ii) << endl;
   }
  report << "Observed Length Prop survey mature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);
    report << ii <<" " <<obsPrNatZ_Srv_msxnz(2,2,1,i) << endl;
   }
  report << "Predicted length prop survey mature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_OS_mxyz(2,1,ii) << endl;
   }
 
  report << "Observed Length Prop survey mature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++) report << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(2,1,2,i) << endl;
  report << "Predicted length prop survey mature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++)
   {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_NS_mxyz(2,2,ii) << endl;
   }
  report << "Observed Length Prop survey mature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++) report << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(2,2,2,i) << endl;
  report << "Predicted length prop survey mature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++){
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_OS_mxyz(2,2,ii) << endl;
   }
//for females don't have length data in first four years first year is 1974
     report << "Observed Length Prop survey all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++){
    ii=yrsObsZCsSrv_n(i);
    report << ii <<" " <<obsPrNatZ_Srv_msxnz(1,1,1,i)+obsPrNatZ_Srv_msxnz(2,1,1,i)+obsPrNatZ_Srv_msxnz(2,2,1,i)<< endl;
              tmpp4+=obsPrNatZ_Srv_msxnz(1,1,1,i)+obsPrNatZ_Srv_msxnz(2,1,1,i)+obsPrNatZ_Srv_msxnz(2,2,1,i);
   }
  report << "Predicted length prop survey all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++) {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_NS_mxyz(1,1,ii)+modSrvPrNatZ_NS_mxyz(2,1,ii)+modSrvPrNatZ_OS_mxyz(2,1,ii) << endl;
    tmpp1+=modSrvPrNatZ_NS_mxyz(1,1,ii)+modSrvPrNatZ_NS_mxyz(2,1,ii)+modSrvPrNatZ_OS_mxyz(2,1,ii);
   }
  report << "Observed Length Prop survey all males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++) {
    ii=yrsObsZCsSrv_n(i);
    report << ii <<" " <<obsPrNatZ_Srv_msxnz(1,1,2,i)+obsPrNatZ_Srv_msxnz(1,2,2,i)+obsPrNatZ_Srv_msxnz(2,1,2,i)+obsPrNatZ_Srv_msxnz(2,2,2,i)<< endl;
         tmpp2+=obsPrNatZ_Srv_msxnz(1,1,2,i)+obsPrNatZ_Srv_msxnz(1,2,2,i)+obsPrNatZ_Srv_msxnz(2,1,2,i)+obsPrNatZ_Srv_msxnz(2,2,2,i);
   }
  report << "Predicted length prop survey all males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
  for (i=1; i<=nObsZCsSrv; i++) {
    ii=yrsObsZCsSrv_n(i);  
    report << ii << " " << modSrvPrNatZ_NS_mxyz(1,2,ii)+modSrvPrNatZ_NS_mxyz(2,2,ii)+modSrvPrNatZ_OS_mxyz(2,2,ii) << endl;
    tmpp3+=modSrvPrNatZ_NS_mxyz(1,2,ii)+modSrvPrNatZ_NS_mxyz(2,2,ii)+modSrvPrNatZ_OS_mxyz(2,2,ii);
   }
  report << "Sum of predicted prop survey all females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
            report <<tmpp1<<endl;
  report << "Sum of predicted prop survey all males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
            report <<tmpp3<<endl;
  report << "Sum of Observed prop survey all females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
            report <<tmpp4<<endl;
  report << "Sum of Observed prop survey all males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'" << endl;
            report <<tmpp2<<endl;

  report << "Predicted mean postmolt length females:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << mean_length(1) << endl;
  report << "Predicted mean postmolt length males:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << mean_length(2)<<endl; 
  report << "af1: 'females'" << endl;
  report << af1 << endl;
  report << "am1: 'males'" << endl;
  report << am1 << endl;
  report << "bf1: 'females'" << endl;
  report << bf1 << endl;
  report << "bm1: 'males'" << endl;
  report << bm1 << endl;
  report<<"Predicted probability of maturing females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<<endl;
  report<<maturity_est(1)<<endl;
  report<<"Predicted probability of maturing males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<<endl;
  report<<maturity_est(2)<<endl;
  report<<"molting probs female: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<<endl;
  report<<moltp(1)<<endl;
  report<<"molting probs male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report<<moltp(2)<<endl;
  report <<"Molting probability mature males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<moltp_mat(2)<<endl;
  
  //fishery info
  report << "observed TCF retained catch biomass: seq(1965,"<<endyr-1<<")" << endl;
  report << obsRetCatchBio(1965,endyr-1) << endl;
  report << "predicted TCF retained catch biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << predRetBioMortTCFM_y(styr,endyr-1)<<endl;
  report << "predicted TCF retained new catch biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << (catch_male_ret_new*wtm_z)(styr,endyr-1)<<endl;
  report << "predicted TCF retained old catch biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << (catch_male_ret_old*wtm_z)(styr,endyr-1)<<endl;
  
  report << "observed TCF total male catch mortality biomass: ( "<<yrsObsDscTCF_n<<" )" << endl;
  report << obsTotBioMortTCFM << endl;
  report << "predicted TCF total male catch mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << predTotBioMortTCFM_y(styr,endyr-1) << endl;
  report << "predicted TCF total new male catch mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << (catch_lmale_new*wtm_z)(styr,endyr-1) << endl;
  report << "predicted TCF total old male catch mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << (catch_lmale_old*wtm_z)(styr,endyr-1) << endl;
  
  report << "observed TCF discard male mortality biomass: ( "<<yrsObsDscTCF_n<<" )"<<endl;
  report << (obsTotBioMortTCFM-obsRetCatchBio(yrsObsDscTCF_n)) <<endl;
  report << "observed TCF male discard mortality biomass: ("<<yrsObsDscTCF_n<<")" << endl;
  report << obsDscBioMortTCF(MALE) << endl;
  report << "predicted TCF male discard mortality biomass1: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << predTotBioMortTCFM_y - predRetBioMortTCFM_y<< endl;
  report << "observed TCF female discard mortality biomass: ( "<<yrsObsDscTCF_n<<" )" << endl;
  report << obsDscBioMortTCF(FEMALE) << endl;
  report << "predicted TCF female discard mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;  
  report << predDscBioMortTCFF_y(styr,endyr-1) << endl;
  
  report << "observed SCF female discard mortality biomass: ( "<<yrsObsDscSCF<<" )" << endl;
  report << obsDscBioSCF_xn(FEMALE)<<endl;
  report << "predicted SCF female discard mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << predDscBioMortSCF_xy(FEMALE) << endl;
  report << "observed SCF male discard mortality biomass: ( "<<yrsObsDscSCF<<" )" << endl;
  report << obsDscBioSCF_xn(MALE)<<endl;
  report << "predicted SCF male discard mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << predDscBioMortSCF_xy(MALE) << endl;  
  
  report << "observed RKC female discard mortality biomass: ( "<<yrsObsDscRKF<<" )" << endl;
  report << obsDscBioRKF_xn(FEMALE) << endl;
  report << "predicted RKC female discard mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << predDscBioMortRKF_xy(FEMALE) << endl;
  report << "observed RKC male discard mortality biomass: ( "<<yrsObsDscRKF<<" )" << endl;
  report << obsDscBioRKF_xn(MALE) << endl;
  report << "predicted RKC male discard mortality biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << predDscBioMortRKF_xy(MALE) << endl;
  
  report << "observed GTF discard mortality biomass: ( "<<yrsObsDscGTF<<" )"<<endl;
  report << obsDscBioMortGTF<<endl;
  report << "predicted GTF discard mortality biomass: seq("<<styr<<","<<endyr-1<<")"<<endl;
  report <<predDscBioMortGTF_y<<endl;
  
  report << "predicted total male catch biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report<<predTotBioMortTCFM_y(styr,endyr-1)+predDscBioMortRKF_xy(  MALE)(styr,endyr-1)+predDscBioMortSCF_xy(  MALE)(styr,endyr-1)+predDscBioMortGTF_y(styr,endyr-1)/2.0<<endl;
  report << "predicted total female catch biomass: seq("<<styr<<","<<endyr-1<<")" << endl;
  report<<predDscBioMortTCFF_y(styr,endyr-1)+predDscBioMortRKF_xy(FEMALE)(styr,endyr-1)+predDscBioMortSCF_xy(FEMALE)(styr,endyr-1)+predDscBioMortGTF_y(styr,endyr-1)/2.0<<endl;
  
  //fishing mortality rates
  report << "estimated annual total male directed fishing mortality: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << fTCF_xy(MALE)(styr,endyr-1) << endl;
  report << "estimated annual male snow fishing mortality: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << fSCF_xy(MALE)(styr,endyr-1) << endl;
  report << "estimated annual male red king fishing mortality: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << fRKF_xy(MALE)(styr,endyr-1) << endl;
  report << "estimated annual male fishing mortality trawl bycatch: seq("<<styr<<","<<endyr-1<<")" << endl;
  report << fGTF_xy(MALE)(styr,endyr-1) <<endl;
  report << "estimated annual total fishing mortality: seq("<<styr<<","<<endyr-1<<")" << endl;
  for(i=styr;i<endyr;i++) report << fmTCFM_syz(1,i)(nZBs)+ fmGTF_xyz(2,i)(nZBs)+fmSCF_xyz(MALE,i)(nZBs)+fmRKF_xyz(2,i)(nZBs) <<" "; report<< endl;
  report <<"retained fmTCFM_syz: seq("<<styr<<","<<endyr-1<<")" << endl;
  for(i=styr;i<endyr;i++) report <<fmTCFR_syz(1,i)(nZBs)<<" "; report<<endl;
//   report <<"ghl: seq(1979,"<<endyr<<")" << endl;
//   report <<catch_ghl/2.2<<endl;
  report << "estimated annual fishing mortality females pot: seq("<<styr<<","<<endyr-1<<")" << endl;
  for(i=styr;i<endyr;i++) report << fmTCFF_yz(i)(nZBs) <<" "; report<<endl;
//recruits in the model are 1978 to 2004, the 1978 recruits are those that enter the population
//in spring of 1979, before the 1979 survey - since using the survey as the start of the year
// in the model spring 1979 is stil 1978.  the last recruits are 2003 that come in spring 2004
  report << "estimated number of recruits female (1000's): seq("<<styr<<","<<endyr<<")" << endl;      //IMPORTANT CHANGE: was 'seq("<<styr+1<<","<<endyr<<")"'
  report << rec_y <<endl;//IMPORTANT CHANGE: was i<=(endyr-1)
  report <<endl<< "estimated number of recruits male (1000's): seq("<<styr<<","<<endyr<<")" << endl;//IMPORTANT CHANGE: was 'seq("<<styr+1<<","<<endyr<<")"'
  report << rec_y <<endl;//IMPORTANT CHANGE: was i<=(endyr-1)
  
  report<<"distribution of recruits to length bins: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<<endl;
  report<<rec_len<<endl;
//  report<<"fishery total selectivity new shell 50% parameter: seq("<<1981+1<<","<<endyr<<")"<<endl;
//  report <<fish_sel50_mn(1981,endyr-1)<<endl;
//  report <<"fishery total selectivity old shell 50% parameter: seq("<<styr+1<<","<<endyr<<")"<<endl;
//  report <<mfexp(log_avg_sel50_mo+log_sel50_dev_mo)(styr,endyr-1)<<endl;
  report << "selectivity fishery total new males styr to 1991: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selTCFM(1,1990) << endl;
//  report << "selectivity fishery total new males 1992 to 1996: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
//  report << selTCFM(1,1991) << endl;
  report << "selectivity fishery total new males 1991 to 1996: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
   for(i=1991;i<=1996; i++) { report << selTCFM(1,i) << endl;}
  report << "selectivity fishery total new males 2005 to present: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
   for(i=2005;i<endyr; i++) { report << selTCFM(1,i) << endl;}
//  report << "selectivity fishery total old males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
//  report << selTCFM(2) << endl;
  //wts 2014-04-25 changed from selTCFR to retFcn
  report << "retention function new males styr to 1991: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;  
  report << retFcn(1,1990) << endl;
//  report << "selectivity fishery ret new males 1992 to 1996: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;  
//  report << selTCFR(1,1991) << endl;
  //wts 2014-04-25 changed from selTCFR to retFcn
  report << "retention function new males 1991 to 1996: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;  
  for(i=1991;i<=1996; i++) { 
      for(j=1; j<=nZBs; j++){report <<retFcn(1,i,j)<<" ";}
      report<<endl;
   }
  //wts 2014-04-25 changed from selTCFR to retFcn
  report << "retention function new males 2005 to present: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;  
  for(i=2005;i<endyr; i++) {        //IMPORTANT CHANGE: was <=endyr 
    for(j=1; j<=nZBs; j++){report <<retFcn(1,i,j)<<" ";}
    report<<endl;
   }
  report <<"FINAL YEAR retention curve males new: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<retFcn(NEW_SHELL,endyr-1)<<endl;
  report <<"FINAL YEAR retention curve males old: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<retFcn(OLD_SHELL,endyr-1)<<endl;
  report << "selectivity discard females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<selTCFF<<endl;
  report << "selectivity trawl females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<selGTF(1,FEMALE)<<endl;
  report <<selGTF(2,FEMALE)<<endl;
  report <<selGTF(3,FEMALE)<<endl;
  report << "selectivity trawl males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<selGTF(1,MALE)<<endl;
  report <<selGTF(2,MALE)<<endl;
  report <<selGTF(3,MALE)<<endl;
  report << "selectivity snow females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<selSCF(1,FEMALE)<<endl;
  report <<selSCF(2,FEMALE)<<endl;
  report <<selSCF(3,FEMALE)<<endl;
  report << "selectivity snow males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<selSCF(1,MALE)<<endl;
  report <<selSCF(2,MALE)<<endl;
  report <<selSCF(3,MALE)<<endl;
  report << "selectivity redk females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<selRKF(1,FEMALE)<<endl;
  report <<selRKF(2,FEMALE)<<endl;
  report <<selRKF(3,FEMALE)<<endl;
  report << "selectivity redk males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report <<selRKF(1,MALE)<<endl;
  report <<selRKF(2,MALE)<<endl;
  report <<selRKF(3,MALE)<<endl;
  report << "selectivity survey females 1974 to 1981: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selSrv2(FEMALE) << endl;
  report << "selectivity survey males 1974 to 1981: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selSrv2(MALE) << endl;
  report << "selectivity survey females 1982 to 1987: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selSrv2a(FEMALE) << endl;
  report << "selectivity survey males 1982 to 1987: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selSrv2a(MALE) << endl;
  report << "selectivity survey females 1988 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selSrv3(FEMALE) << endl;
  report << "selectivity survey males 1988 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selSrv3(MALE) << endl;
  report << "numbers of mature females by age and length: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for(i=styr;i<=endyr;i++)
   { report << natlength_mnew(FEMALE,i)<<endl; report << natlength_mold_age(FEMALE,i)<<endl; }
  
  report << "numbers of mature males by age and length: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  for(i=styr;i<=endyr;i++) 
   { report << natlength_mnew(MALE,i)<<endl; report << natlength_mold_age(MALE,i)<<endl; }
  report << "pred_sexr population: seq("<<styr<<","<<endyr<<")" << endl;
  report << predpop_sexr <<endl;
  report << "pred_sexr survey: seq("<<styr<<","<<endyr<<")" << endl;
  report << preds_sexr <<endl;
  
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"Size transition matrices"<<endl;
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"length - length transition matrix Females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report<<len_len(1)<<endl;  
  report<<"length - length transition matrix Males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report<<len_len(2)<<endl;  
  
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"Effective N's"<<endl;
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"effective N survey lengths immature new shell female "<<endl;
  report<<effn_srv1(IMMATURE,NEW_SHELL,FEMALE)<<endl;
  report<<" effective N survey lengths mature new shell female "<<endl;
  report<<effn_srv1(MATURE,NEW_SHELL,FEMALE)<<endl;
  report<<" effective N survey lengths mature old shell female "<<endl;
  report<<effn_srv1(MATURE,OLD_SHELL,FEMALE)<<endl;
  report<<" effective N survey lengths immature new shell male "<<endl;
  report<<effn_srv1(IMMATURE,NEW_SHELL,MALE)<<endl;
  report<<" effective N survey lengths immature old shell male "<<endl;
  report<<effn_srv1(IMMATURE,OLD_SHELL,MALE)<<endl;
  report<<" effective N survey lengths mature new shell male "<<endl;
  report<<effn_srv1(MATURE,NEW_SHELL,MALE)<<endl;
  report<<" effective N survey lengths mature old shell male "<<endl;
  report<<effn_srv1(MATURE,OLD_SHELL,MALE)<<endl;
  report<<" effective N retained lengths new, old shell "<<endl;
  report<<effn_fish_ret<<endl;
  report<<" effective N total lengths new, old shell"<<endl;
  report<<effn_fish_tot<<endl;
  
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"Exploitation rates"<<endl;
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"male new shell total pot fishery exploitation rates"<<endl;
  report<<1-mfexp(-1.0*fmTCFM_syz(NEW_SHELL))<<endl;
  report<<"male old shell total pot fishery exploitation rates"<<endl;
  report<<1-mfexp(-1.0*fmTCFM_syz(OLD_SHELL))<<endl;
  
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"Numbers at fish time and catch"<<endl;
  report<<"----------------------------------------------------------------------------"<<endl;
  report<<"numbers new shell males at time of pop fishery"<<endl;
  report<<natl_new_fishtime(MALE)<<endl;
  report<<"numbers old shell males at time of pop fishery"<<endl;
  report<<natl_old_fishtime(MALE)<<endl;
  report<<"total catch in numbers new shell males"<<endl;
  report<<catch_lmale_new<<endl;
  report<<"total catch in numbers old shell males"<<endl;
  report<<catch_lmale_old<<endl;
  report<<"retained catch in numbers new shell males"<<endl;
  report<<catch_male_ret_new<<endl;
  report<<"retained catch in numbers old shell males"<<endl;
  report<<catch_male_ret_old<<endl;
  report<<"observed retained catch new shell males"<<endl;
  for (i=1; i<=nObsRetZCsTCF; i++) 
   report<<yrsObsRetZCsTCF_n(i)<<" "<<obsPrNatZ_TCFR_sn(NEW_SHELL,i)*obsRetCatchNum(yrsObsRetZCsTCF_n(i))<<endl;
  report<<"observed retained catch old shell males"<<endl;
  for (i=1; i<=nObsRetZCsTCF; i++) 
   report<<yrsObsRetZCsTCF_n(i)<<" "<<obsPrNatZ_TCFR_sn(OLD_SHELL,i)*obsRetCatchNum(yrsObsRetZCsTCF_n(i))<<endl;
  
  // stuff for input to projection model
  report<<"#--------------------------------------------------------------------"<<endl;
  report<<"#FOR PROJECTION MODEL------------------------------------------------"<<endl;
  report<<"#--------------------------------------------------------------------"<<endl;
  report<<"#number of length bins"<<endl;
  report<<nZBs<<endl;
  report<<"#Nat mort immature female/male"<<endl;
  report<<M_msx(IMMATURE,NEW_SHELL)<<endl;
  report<<"#nat mort mature new shell female/male"<<endl;
  report<<M_msx(MATURE,NEW_SHELL)<<endl;
  report<<"#nat mort mature old shell female/male"<<endl;
  report<<M_msx(MATURE,OLD_SHELL)<<endl;
  report<<"#constant recruitment"<<endl;
  report<<"1000000"<<endl;
  report<<"#average of last 4 years selTCFM total male new old shell"<<endl;
  report<<(selTCFM(1,endyr-4)+selTCFM(1,endyr-3)+selTCFM(1,endyr-2)+selTCFM(1,endyr-1))/4.0<<endl;
  report<<(selTCFM(1,endyr-4)+selTCFM(2,endyr-3)+selTCFM(2,endyr-2)+selTCFM(2,endyr-1))/4.0<<endl;
  report<<"#average of last 4 years selTCFM retained curve male new old shell"<<endl;
  report<<(selTCFR(1,endyr-4)+selTCFR(1,endyr-3)+selTCFR(1,endyr-2)+selTCFR(1,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
  report<<(selTCFR(2,endyr-4)+selTCFR(2,endyr-3)+selTCFR(2,endyr-2)+selTCFR(2,endyr-1))/4.0<<endl;
  report<<"#trawl selectivity female male"<<endl;
  report<<selGTF(3)<<endl;
  report<<"#female pot discard selectivity"<<endl;
  report<<selTCFF<<endl;
  report <<"#selectivity snow females"<< endl;
  report <<selSCF(3,1)<<endl;
  report << "#selectivity snow males"<< endl;
  report <<selSCF(3,2)<<endl;
  report <<"#selectivity redk females"<< endl;
  report <<selRKF(3,1)<<endl;
  report <<"#selectivity redk males"<< endl;
  report <<selRKF(3,2)<<endl; 
  report<<"#maturity curve new shell female male"<<endl;
  report<<maturity_est(1)<<endl;
  report<<maturity_est(2)<<endl;
  report<<"#maturity curve old shell female male"<<endl;
  report<<obsAvgMatOS_xz<<endl;
  report<<"#molting probability immature female male"<<endl;
  report<<moltp<<endl;
  report<<"#molting probability mature female male"<<endl;
  report<<moltp_mat<<endl;
  report<<"#prop recruits to new shell"<<endl;
  report<<prNewShellRecruits<<endl;
  report<<"#distribution of recruits to length bins"<<endl;
  report<<rec_len<<endl;
  report<<"#time of catch in fraction of year from survey - 7 months"<<endl;
  report<<mdptFshs_y(endyr-1)<<endl;//IMPORTANT CHANGE; was endyr
  report<<"#number at length new shell females males at time of fishery endyr from model"<<endl;
  report<<natl_new_fishtime(1,endyr)<<endl;
  report<<natl_new_fishtime(2,endyr)<<endl;
  report<<"#number at length old shell females males at time of fishery endyr from model"<<endl;
  report<<natl_old_fishtime(1,endyr)<<endl;
  report<<natl_old_fishtime(2,endyr)<<endl;
  report<<"#last year male spawning biomass (1000's t)"<<endl;
  report<<mspbio(endyr)<<endl;
  report<<"#last year female spawning biomass"<<endl;
  report<<fspbio(endyr)<<endl;
  report<<"#last year male spawning biomass at matingtime (endyr-1)"<<endl;
  report<<mspbio_matetime(endyr-1)<<endl;//IMPORTANT CHANGE; was endyr
  report<<"#last year female spawning biomass at matingtime (endyr-1)"<<endl;
  report<<fspbio_matetime(endyr-1)<<endl;//IMPORTANT CHANGE; was endyr
  report<<"#numbers at length immature new shell female male last year"<<endl;
  report<<natlength_inew(1,endyr)<<endl;
  report<<natlength_inew(2,endyr)<<endl;
  report<<"#numbers at length immature old shell female male last year"<<endl;
  report<<natlength_iold(1,endyr)<<endl;
  report<<natlength_iold(2,endyr)<<endl;
  report<<"#numbers at length mature new shell female male last year"<<endl;
  report<<natlength_mnew(1,endyr)<<endl;
  report<<natlength_mnew(2,endyr)<<endl;
  report<<"#numbers at length mature old shell female male last year"<<endl;
  report<<natlength_mold(1,endyr)<<endl;
  report<<natlength_mold(2,endyr)<<endl;
  report<<"#weight at length female juvenile (t)"<<endl;
  report<<wtf_mz(IMMATURE)<<endl;
  report<<"#weight at length female mature (t)"<<endl;
  report<<wtf_mz(MATURE)<<endl;
  report<<"#weight at length male"<<endl;
  report<<wtm_z<<endl;
  report<<"#length-length transition matrix"<<endl;
  report<<len_len<<endl;
  report<<"#female discard pot fishing fmTCFM_syz average last 5 yrs (removed 20150601)"<<endl;
//  report<<mean(fmortdf(endyr-5,endyr-1))<<endl;
  report<<"#fSCF_xy male average last 5 yrs"<<endl;
  report<<mean(fSCF_xy(MALE)(endyr-5,endyr-1))<<endl;
  report<<"#fRKF_xy male average last 5 yrs"<<endl;
  report<<mean(fRKF_xy(MALE)(endyr-5,endyr-1))<<endl;
  report<<"#fGTF_xy(MALE) average last 5 yrs"<<endl;
  report<<mean(fGTF_xy(MALE)(endyr-5,endyr-1))<<endl;
  report<<"#number of recruits from the model styr to endyr"<<endl;//was endyr-1
  report<<endyr-styr+1<<endl;
  report<<"#number of recruits for avg to estimate B35%"<<endl;
  report<<endyr-1982+1<<endl;
  report <<"#recruitments female, male start year to endyr from model (1000's)" << endl;//was endyr-1
  report << rec_y <<endl;//was endyr-1
  report << "#recruitments male, male start 1960 to endyr from model (1000's)" << endl;//was endyr-1
  report << rec_y <<endl;//was endyr-1
  
  report<<"#male spawning biomass at matetime for endyr-5 to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
  report<<mspbio_matetime(endyr-5,endyr-1)<<endl;
  report<<"#male spawning biomass at matetime for str year to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
  report<<mspbio_matetime(styr,endyr-1)<<endl;
  report <<"#selectivity survey males 1989 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5','137.5','142.5','147.5','152.5','157.5','162.5','167.5','172.5','177.5','182.5'"<< endl;
  report << selSrv3(2) << endl;
  
  CheckFile<<"end writeReport"<<endl;
  
// ==========================================================================
FUNCTION void writeMyProjectionFile(ofstream& os)
//    cout<<"starting writeToR"<<endl;
  // stuff for input to projection model
      os<<"#--------------------------------------------------------------------"<<endl;
      os<<"#FOR WTS PROJECTION MODEL--------------------------------------------"<<endl;
      os<<"#--------------------------------------------------------------------"<<endl;
      os<<styr <<tb<<tb<<"#start year for assessment model"<<endl;
      os<<endyr<<tb<<tb<<"#end year for assessment model/start year for projections"<<endl;
      os<<10   <<tb<<tb<<"#number of years for projections"<<endl;
      os<<1000 <<tb<<tb<<"#number of projections to make"<<endl;
      os<<nZBs<<tb<<tb<<"#number of size bins in model"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#reference point calculations"<<endl;
      os<<"#---------------------------"<<endl;
      os<<0.35<<tb<<tb<<"#target SBPR reduction ratio (e.g. 0.35)"<<endl;
      os<<2*mean(rec_y(1982,endyr))<<tb<<tb<<"#total average recruitment for BXX/Bmsy calculation (1000's of recruits)"<<endl;
      os<<mspbio_matetime(endyr-1)<<tb<<tb<<"#'current' spawning biomass (MMB, 1000's t) "<<endl;
      os<<"???"<<tb<<tb<<"#cv of 'current' spawning biomass"<<endl;
      os<<1    <<tb<<tb<<"#harvest strategy"<<endl;
      os<<0.1  <<tb<<tb<<"#alpha for control rule"<<endl;
      os<<0.25 <<tb<<tb<<"#beta  for control rule"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#SRR info"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"1"   <<tb<<tb<<"#srType  : flag indicating type of SRR"<<endl;
      os<<"0"   <<tb<<tb<<"#recGamma: recruitment autocorrelation parameter"<<endl;
      os<<-1    <<tb<<tb<<"#recDep  : recruitment depensation flag"<<endl;
      os<<-1    <<tb<<tb<<"#inpRecH : input value for steepness"<<endl;
      os<<-1    <<tb<<tb<<"#inpRecR0: input value for virgin recruitment"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#Recruitment and spawning biomass time series"<<endl;
      os<<"#---------------------------"<<endl;
      os<<1982  <<tb<<tb<<"#mMnYrForRecAvg: min assessment model year for averaging recruitment to estimate BXX, Bmsy"<<endl;
      os<<endyr <<tb<<tb<<"#mMxYrForRecAvg: max assessment model year for averaging recruitment to estimate BXX, Bmsy"<<endl;
      os <<"#asmtModRec(nSXs,mMnYr,mMxYr): unlagged recruitments female, male start year to endyr from model (1000's)" << endl;//was endyr-1
      os << rec_y <<endl;//females; was endyr-1
      os << rec_y <<endl;//  males; was endyr-1     
      os<<"#asmtModSpB(mMnYr,mMxYr-1): male spawning biomass at matetime (1000's t) for str year to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
      os<<mspbio_matetime(styr,endyr-1)<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#Pop info in last year of assessment model"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#numbers at size immature new shell female, male in final year (1000's)"<<endl;
      os<<natlength_inew(FEMALE,endyr)<<endl;
      os<<natlength_inew(  MALE,endyr)<<endl;
      os<<"#numbers at length immature old shell female male last year (1000's)"<<endl;
      os<<natlength_iold(FEMALE,endyr)<<endl;
      os<<natlength_iold(  MALE,endyr)<<endl;
      os<<"#numbers at length mature new shell female male last year (1000's)"<<endl;
      os<<natlength_mnew(FEMALE,endyr)<<endl;
      os<<natlength_mnew(  MALE,endyr)<<endl;
      os<<"#numbers at length mature old shell female male last year (1000's) "<<endl;
      os<<natlength_mold(FEMALE,endyr)<<endl;
      os<<natlength_mold(  MALE,endyr)<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#Fisheries info"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#time of catch in fraction of year from survey - 7 months"<<endl;
      os<<mdptFshs_y(endyr-1)<<endl;//IMPORTANT CHANGE; was endyr
      
      os<<0<<tb<<tb<<"#inpFmTCF: input F for directed Tanner crab fishing mortality"<<endl;
      os<<mean(fSCF_xy(MALE)(endyr-5,endyr-1))<<tb<<tb<<"#inpFmSCF: input male F for snow crab fishing mortality"<<endl;
      os<<mean(fRKF_xy(MALE)(endyr-5,endyr-1))<<tb<<tb<<"#inpFmRKF: input male F for BBRKC  fishing mortality"<<endl;
      os<<mean(fGTF_xy(MALE)(endyr-5,endyr-1))<<tb<<tb<<"#inpFmGTF: input male F for groundfish fishery fishing mortality"<<endl;
      
      os<<"#selTCF_TotMale(nSCs,nSXs): average of last 4 years selTCFM total male new old shell"<<endl;
      os<<(selTCFM(NEW_SHELL,endyr-4)+selTCFM(NEW_SHELL,endyr-3)+selTCFM(NEW_SHELL,endyr-2)+selTCFM(NEW_SHELL,endyr-1))/4.0<<endl;
      os<<(selTCFM(OLD_SHELL,endyr-4)+selTCFM(OLD_SHELL,endyr-3)+selTCFM(OLD_SHELL,endyr-2)+selTCFM(OLD_SHELL,endyr-1))/4.0<<endl;
      os<<"#selTCF_RetMale(nSCs,nSXs): average of last 4 years selTCFM retained curve male new old shell"<<endl;
      os<<(selTCFR(NEW_SHELL,endyr-4)+selTCFR(NEW_SHELL,endyr-3)+selTCFR(NEW_SHELL,endyr-2)+selTCFR(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
      os<<(selTCFR(OLD_SHELL,endyr-4)+selTCFR(OLD_SHELL,endyr-3)+selTCFR(OLD_SHELL,endyr-2)+selTCFR(OLD_SHELL,endyr-1))/4.0<<endl;
      os<<"#selTCF_TotMaleEast(nSCs,nSXs): set same as average total"<<endl;
      os<<(selTCFM(NEW_SHELL,endyr-4)+selTCFM(NEW_SHELL,endyr-3)+selTCFM(NEW_SHELL,endyr-2)+selTCFM(NEW_SHELL,endyr-1))/4.0<<endl;
      os<<(selTCFM(OLD_SHELL,endyr-4)+selTCFM(OLD_SHELL,endyr-3)+selTCFM(OLD_SHELL,endyr-2)+selTCFM(OLD_SHELL,endyr-1))/4.0<<endl;
      os<<"#selTCF_RetMaleEast(nSCs,nSXs): set same as avg retained"<<endl;
      os<<(selTCFR(NEW_SHELL,endyr-4)+selTCFR(NEW_SHELL,endyr-3)+selTCFR(NEW_SHELL,endyr-2)+selTCFR(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
      os<<(selTCFR(OLD_SHELL,endyr-4)+selTCFR(OLD_SHELL,endyr-3)+selTCFR(OLD_SHELL,endyr-2)+selTCFR(OLD_SHELL,endyr-1))/4.0<<endl;
      os<<"#selTCF_TotMaleWest(nSCs,nSXs): set same as average total"<<endl;
      os<<(selTCFM(NEW_SHELL,endyr-4)+selTCFM(NEW_SHELL,endyr-3)+selTCFM(NEW_SHELL,endyr-2)+selTCFM(NEW_SHELL,endyr-1))/4.0<<endl;
      os<<(selTCFM(OLD_SHELL,endyr-4)+selTCFM(OLD_SHELL,endyr-3)+selTCFM(OLD_SHELL,endyr-2)+selTCFM(OLD_SHELL,endyr-1))/4.0<<endl;
      os<<"#selTCF_RetMaleWest(nSCs,nSXs): SET SAME AS AVG RETAINED, BUT SHIFTED TO LOWER END BY 10 mm"<<endl;
      os<<(selTCFR(NEW_SHELL,endyr-4)+selTCFR(NEW_SHELL,endyr-3)+selTCFR(NEW_SHELL,endyr-2)+selTCFR(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
      os<<(selTCFR(OLD_SHELL,endyr-4)+selTCFR(OLD_SHELL,endyr-3)+selTCFR(OLD_SHELL,endyr-2)+selTCFR(OLD_SHELL,endyr-1))/4.0<<endl;
      os<<"#selTCF_Female(nZs): selectivity for females in directed fishery"<<endl;
      os<<selTCFF<<endl;
      os<<"#selSCF(nSXs,nZs): selectivity in snow crab fishery"<<endl;
      os<<selSCF(3,FEMALE)<<endl;
      os<<selSCF(3,  MALE)<<endl;
      os<<"#selRKF(nSXs,nZs): selectivity in BBRKC fishery"<<endl;
      os<<selRKF(3,FEMALE)<<endl;
      os<<selRKF(3,  MALE)<<endl;      
      os<<"#selGTF(nSXs,nZs): selectivity in groundfish fishery"<<endl;
      os<<selGTF(3)<<endl;
      
      os<<"#---------------------------"<<endl;
      os<<"#Biological info"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#M_f(nSCs,nMSs): natural mortality for females"<<endl;
      os<<M_msx(IMMATURE,NEW_SHELL,FEMALE)<<tb<<M_msx(MATURE,NEW_SHELL,FEMALE)<<endl;
      os<<M_msx(IMMATURE,OLD_SHELL,FEMALE)<<tb<<M_msx(MATURE,OLD_SHELL,FEMALE)<<endl;
      os<<"#M_m(nSCs,nMSs): natural mortality for males"<<endl;
      os<<M_msx(IMMATURE,NEW_SHELL,  MALE)<<tb<<M_msx(MATURE,NEW_SHELL,  MALE)<<endl;
      os<<M_msx(IMMATURE,OLD_SHELL,  MALE)<<tb<<M_msx(MATURE,OLD_SHELL,  MALE)<<endl;
      os<<"#weight at length female juvenile (t)"<<endl;
      os<<wtf_mz(IMMATURE)<<endl;
      os<<"#weight at length female mature (t)"<<endl;
      os<<wtf_mz(MATURE)<<endl;
      os<<"#weight at length male (t)"<<endl;
      os<<wtm_z<<endl;
      os<<"#tmZtoZ_xzz: size transition matrix"<<endl;
      os<<len_len<<endl;      
      os<<"#prMatNS(nSXs,nZs): maturity curve new shell female male"<<endl;
      os<<maturity_est(FEMALE)<<endl;
      os<<maturity_est(MALE)<<endl;
      os<<"#prMoltImm(nSXs,nZs): molting probability immature female male"<<endl;
      os<<moltp<<endl;
      os<<"#prMoltMat(nSXs,nZs): molting probability mature female male"<<endl;
      os<<moltp_mat<<endl;
      os<<0.5     <<tb<<tb<<"#recPropAsMale: proportion recruiting as males"<<endl;
      os<<prNewShellRecruits<<tb<<tb<<"#recPropAsNewShell: prop recruits to new shell"<<endl;
      os<<"#recPropAtZ(nZs): distribution of recruits to length bins"<<endl;
      os<<rec_len<<endl;
      os<<"#propEast(nZs): proportion of population at size east of 166W"<<endl;
      os<<"???????????????????????????????"<<endl;
  
// ==========================================================================
//"Old style" R file
FUNCTION void writeToR(ofstream& R_out)
//    cout<<"starting writeToR"<<endl;

        int ii;
        dvar_vector preds_sexr(styr,endyr);
        dvar_matrix tmpo(1,2,styr,endyr);
        dvar_matrix tmpp(1,2,styr,endyr);
        dvar_vector obs_tmp(styr,endyr);
        dvariable ghl,ghl_number;
        dvariable hrate;
        dvar_vector totcatch(styr, endyr-1);  
        tmpp1=0.0;
        tmpp2=0.0;
        tmpp3=0.0;
        tmpp4=0.0;
        
        Misc_output();
        
        R_out<<"$styr" <<endl<<styr<<endl;
        R_out<<"$endyr"<<endl<<endyr<<endl;
        R_out<<"$obsyr"<<endl<<ptrMDS->pTSD->yrsAbund[1]<<endl;
        R_out<<"$pltyr"<<endl<<1969<<endl;
        R_out<<"$length.bins"<<endl<<length_bins<<endl;
        
        R_out<<"$optFM"<<endl<<optFM<<endl;
        
        R_out<<"$years.survey.abundance"                   <<endl<<ptrMDS->pTSD->yrsAbund<<endl;
        R_out<<"$years.obs.retained.catch.directed.fishery"<<endl<<ptrMDS->pTCFR->yrsCatch<<endl;
        R_out<<"$years.obs.total.catch.directed.fishery"   <<endl<<ptrMDS->pTCFD->yrsCatch<<endl;
        R_out<<"$years.obs.bycatch.snow.fishery"           <<endl<<ptrMDS->pSCF->yrsCatch<<endl;
        R_out<<"$years.obs.bycatch.redk.fishery"           <<endl<<ptrMDS->pRKF->yrsCatch<<endl;
        R_out<<"$years.obs.bycatch.trawl.fishery"          <<endl<<ptrMDS->pGTF->yrsCatch<<endl;
        
        R_out << "$Estimated.numbers.of.immature.new.shell.female.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++) R_out <<  i<<" "<<natlength_inew(1,i) << endl;
        R_out << "$Estimated.numbers.of.immature.old.shell.female.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++) R_out <<  i<<" "<<natlength_iold(1,i) << endl;
        R_out << "$Estimated.numbers.of.mature.new.shell.female.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++) R_out <<  i<<" "<<natlength_mnew(1,i) << endl;
        R_out << "$Estimated.numbers.of.mature.old.shell.female.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++)R_out <<  i<<" "<<natlength_mold(1,i) << endl;
        
        R_out << "$Estimated.numbers.of.immature.new.shell.male.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++) R_out << i<<" "<<natlength_inew(2,i) << endl;
        R_out << "$Estimated.numbers.of.immature.old.shell.male.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++) R_out << i<<" "<<natlength_iold(2,i) << endl;
        R_out << "$Estimated.numbers.of.mature.new.shell.male.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++) R_out << i<<" "<<natlength_mnew(2,i) << endl;
        R_out << "$Estimated.numbers.of.mature.old.shell.male.crab.by.length"<< endl;
        for (int i=styr;i<=endyr;i++) R_out << i<<" "<<natlength_mold(2,i) << endl;
        
        R_out << "$Observed.numbers.of.immature.new.shell.female.crab.by.length"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(1,1,1,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
        R_out << "$Observed.numbers.of.mature.new.shell.female.crab.by.length"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,1,1,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
        R_out << "$Observed.numbers.of.mature.old.shell.female.crab.by.length"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,2,1,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
        R_out << "$Observed.numbers.of.immature.new.shell.male.crab.by.length"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(1,1,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
        R_out << "$Observed.numbers.of.immature.old.shell.male.crab.by.length"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(1,2,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
        R_out << "$Observed.numbers.of.mature.new.shell.male.crab.by.length"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,1,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
        R_out << "$Observed.numbers.of.mature.old.shell.male.crab.by.length"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "<<obsPrNatZ_Srv_msxnz(2,2,2,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))<<endl;
        R_out << "$Observed.Survey.Numbers.by.length.females"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" " << obsSrvNum_xyz(1,yrsObsZCsSrv_n(i)) << endl;
        R_out << "$Observed.Survey.Numbers.by.length.males"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" " << obsSrvNum_xyz(2,yrsObsZCsSrv_n(i))<< endl;
        
        R_out << "$Predicted.Survey.Numbers.by.length.females"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "  << modSrvNum_xyz(1,yrsObsZCsSrv_n(i)) << endl;
        R_out << "$Predicted.Survey.Numbers.by.length.males"<< endl;
        for (int i=1; i <= nObsZCsSrv; i++) R_out<<yrsObsZCsSrv_n(i)<<" "  << modSrvNum_xyz(2,yrsObsZCsSrv_n(i)) << endl;
        R_out << "$Predicted.pop.Numbers.by.length.females"<< endl;
        for (int i=styr;i<=endyr;i++) R_out<<i<<" "<< natlength(1,i)<< endl;
        R_out << "$Predicted.pop.Numbers.by.length.males"<< endl;
        for (int i=styr;i<=endyr;i++) R_out<<i<<" "<< natlength(2,i)<< endl;
        
        R_out<<"$observed.number.of.males.greater.than.101.mm"<<endl<<obsSrvLegalNum_n<<endl;
        R_out<<"$observed.biomass.of.males.greater.than.101.mm"<<endl<<obsSrvLegalBio_n<<endl;
        R_out<<"$pop.estimate.numbers.of.males.101"        <<endl<<numLegalMales_y<<endl;
        R_out<<"$estimated.population.biomass.of.males.101"<<endl<<bioLegalMales_y<<endl;
        R_out<<"$estimated.survey.numbers.of.males.101"    <<endl<<modSrvNumLegal_y<<endl;
        R_out<<"$estimated.survey.biomass.of.males.101"    <<endl<<modSrvBioLegal_y<<endl;
        R_out<<"$estimated.biomass.of.males.101.fishtime"  <<endl<<bio_males_gt101<<endl;
        R_out<<"$Observed.survey.biomass"                  <<endl<< obsSrvBio_y(yrsObsSrvBio_n(1),endyr)<<endl;
        R_out<<"$predicted.survey.biomass"                 <<endl<<modSrvBio_xy(1)+modSrvBio_xy(2)<<endl;
        
        //survey numbers
        for(int k=1;k<=2;k++){
            for (int i=styr;i<=endyr;i++){
                 tmpo(k,i)=sum(obsSrvNum_xyz(k,i));
                 tmpp(k,i)=sum(modSrvNum_xyz(k,i));
            }
        }
        R_out << "$Observed.survey.numbers.female"<<endl;
        R_out << tmpo(1)<<endl;
        R_out << "$Observed.survey.numbers.male"<<endl;
        R_out << tmpo(2)<<endl;
        R_out << "$predicted.survey.numbers.female"<<endl;
        R_out << tmpp(1)<<endl;
        R_out << "$predicted.survey.numbers.male"<<endl;
        R_out << tmpp(2)<<endl;
        R_out << "$Observed.survey.female.spawning.biomass"<<endl;
        R_out << obsSrvMatBio_xy(1)<<endl;
        R_out << "$Observed.survey.female.spawning.biomass.cv"<<endl;
        R_out << obsSrvCV_xn(1)<<endl;
        R_out << "$Observed.survey.male.spawning.biomass"<<endl;
        R_out << obsSrvMatBio_xy(2)<<endl;
        R_out << "$Observed.survey.male.spawning.biomass.cv"<<endl;
        R_out << obsSrvCV_xn(2)<<endl;
        R_out << "$Observed.survey.female.new.spawning.numbers"<<endl;
        R_out << obsSrvMatNum_sxy(1,1)<<endl;
        R_out << "$Observed.survey.female.old.spawning.numbers"<<endl;
        R_out << obsSrvMatNum_sxy(2,1)<<endl;
        R_out << "$Observed.survey.male.new.spawning.numbers"<<endl;
        R_out << obsSrvMatNum_sxy(1,2)<<endl;
        R_out << "$Observed.survey.male.old.spawning.numbers"<<endl;
        R_out << obsSrvMatNum_sxy(2,2)<<endl;
        R_out << "$Observed.survey.female.biomass"<<endl;
        R_out << obsSrvBio_xy(1)<<endl;
        R_out << "$Observed.survey.male.biomass"<<endl;
        R_out << obsSrvBio_xy(2)<<endl;
        R_out << "$natural.mortality.immature.females.males" << endl;
        R_out << M_msx(IMMATURE,NEW_SHELL) << endl;
        R_out << "$natural.mortality.mature.new.shell.females.males" << endl;
        R_out << M_msx(  MATURE,NEW_SHELL) << endl;
        R_out << "$natural.mortality.mature.old.shell.females.males" << endl;
        R_out << M_msx(  MATURE,OLD_SHELL) << endl;
        R_out << "$Predicted.Biomass" << endl;
        R_out << pred_bio << endl;
        R_out << "$Predicted.total.population.numbers"<<endl;
        R_out << modPopNum_y<<endl;
        R_out << "$Female.Spawning.Biomass" << endl;
        R_out << fspbio << endl;
        R_out << "$Male.Spawning.Biomass" << endl;
        R_out << mspbio << endl;
        R_out << "$Total.Spawning.Biomass" << endl;
        R_out << fspbio+mspbio << endl;
        R_out << "$Female.Spawning.Biomass.at.fish.time" << endl;
        R_out << fspbio_fishtime << endl;
        R_out << "$Male.Spawning.Biomass.at.fish.time" << endl;
        R_out << mspbio_fishtime << endl;
        R_out << "$Total.Spawning.Biomass.at.fish.time" << endl;
        R_out << fspbio_fishtime+mspbio_fishtime << endl;
        R_out << "$Mating.time.Female.Spawning.Biomass" << endl;
        R_out << fspbio_matetime << endl;
        R_out << "$Mating.time.Male.Spawning.Biomass" << endl;
        R_out << mspbio_matetime << endl;
        R_out << "$Mating.time.Male.old.shell.Spawning.Biomasss" << endl;
        R_out << mspbio_old_matetime << endl;
        R_out << "$Mating.time.female.new.shell.Spawning.Biomass" << endl;
        R_out << fspbio_new_matetime << endl;
        R_out << "$Mating.time.Total.Spawning.Biomass" << endl;
        R_out << fspbio_matetime+mspbio_matetime << endl;
        R_out << "$Mating.time.effective.Female.Spawning.Biomass" << endl;
        R_out << efspbio_matetime << endl;
        R_out << "$Mating.time.effective.Male.Spawning.Biomass.old.shell.only" << endl;
        R_out << emspbio_matetime << endl;
        R_out << "$Mating.time.Total.effective.Spawning.Biomass" << endl;
        R_out << efspbio_matetime+emspbio_matetime << endl;
        R_out << "$Mating.time.male.Spawning.numbers" << endl;
        R_out << mspnum_matetime << endl;
        R_out << "$Mating.time.Female.Spawning.numbers" << endl;
        R_out << fspnum_matetime << endl;
        R_out << "$Mating.time.Male.Spawning.numbers.old.shell.only" << endl;
        R_out << mspnum_old_matetime << endl;
        R_out << "$Mating.time.effective.Female.new.shell.Spawning.biomass" << endl;
        R_out <<efspbio_new_matetime << endl;
        R_out << "$Mating.time.Female.new.shell.Spawning.numbers" << endl;
        R_out << fspnum_new_matetime << endl;
        R_out << "$Predicted.Female.survey.Biomass" << endl;
        R_out << modSrvBio_xy(1) << endl;
        R_out << "$Predicted.Male.survey.Biomass" << endl;
        R_out << modSrvBio_xy(2)<< endl;
        R_out << "$Predicted.Female.survey.mature.Biomass" << endl;
        R_out << modSrvMatBio_xy(FEMALE) << endl;
        R_out << "$Predicted.Male.survey.mature.Biomass" << endl;
        R_out << modSrvMatBio_xy(  MALE)<< endl;
        R_out << "$Predicted.total.survey.mature.Biomass" << endl;
        R_out << modSrvMatBio_xy(FEMALE)+modSrvMatBio_xy(  MALE)<< endl;
        R_out << "$Predicted.Female.survey.new.mature.numbers" << endl;
        R_out << modSrvMatNum_xsy(FEMALE,NEW_SHELL) << endl;
        R_out << "$Predicted.Female.survey.old.mature.numbers" << endl;
        R_out << modSrvMatNum_xsy(FEMALE,OLD_SHELL) << endl;
        R_out << "$Predicted.Male.survey.new.mature.numbers" << endl;
        R_out << modSrvMatNum_xsy(  MALE,NEW_SHELL)<< endl;
        R_out << "$Predicted.Male.survey.old.mature.numbers" << endl;
        R_out << modSrvMatNum_xsy(  MALE,OLD_SHELL)<< endl;
        
        R_out << "$Observed.Prop.fishery.ret.new.males"<< endl;
        for (int i=1; i<=nObsRetZCsTCF; i++) R_out << yrsObsRetZCsTCF_n(i) << " " << obsPrNatZ_TCFR_sn(1,i)<< endl;
        R_out << "$Predicted.length.prop.fishery.ret.new.males" << endl;
        for (int i=1; i<=nObsRetZCsTCF; i++) {
            ii=yrsObsRetZCsTCF_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_TCFR_syz(1,ii)  << endl;
        }
        R_out << "$Observed.Prop.fishery.ret.old.males"<< endl;
        for (int i=1; i<=nObsRetZCsTCF; i++) R_out << yrsObsRetZCsTCF_n(i) << " " << obsPrNatZ_TCFR_sn(2,i)<< endl;
        R_out << "$Predicted.length.prop.fishery.ret.old.males" << endl;
        for (int i=1; i<=nObsRetZCsTCF; i++) {
            ii=yrsObsRetZCsTCF_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_TCFR_syz(2,ii)  << endl;
        }
        
        R_out << "$Observed.Prop.fishery.total.new.males"<< endl;
        for (int i=1; i<=nObsZCsTCFM; i++) R_out << yrsObsZCsTCFM_n(i) << " " << obsPrNatZ_TCFM_snz(1,i) << endl;
        R_out << "$Predicted.length.prop.fishery.total.new.males" << endl;
        for (int i=1; i<=nObsZCsTCFM; i++) {
            ii=yrsObsZCsTCFM_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_TCFM_syz(1,ii)  << endl;
        }
        R_out << "$Observed.Prop.fishery.total.old.males"<< endl;
        for (int i=1; i<=nObsZCsTCFM; i++) R_out << yrsObsZCsTCFM_n(i) << " " << obsPrNatZ_TCFM_snz(2,i) << endl;
        R_out << "$Predicted.length.prop.fishery.total.old.males" << endl;
        for (int i=1; i<=nObsZCsTCFM; i++) {
            ii=yrsObsZCsTCFM_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_TCFM_syz(2,ii)  << endl;
        }
        R_out << "$Observed.Prop.fishery.discard.new.males"<< endl;
        for (int i=1; i<=nObsZCsTCFM; i++) R_out << yrsObsZCsTCFM_n(i) << " " << obs_p_fish_discm(1,i) << endl;
        R_out << "$Observed.Prop.fishery.discard.old.males"<< endl;
        for (int i=1; i<=nObsZCsTCFM; i++) R_out << yrsObsZCsTCFM_n(i) << " " << obs_p_fish_discm(2,i)<< endl;
        
        R_out << "$Observed.length.prop.fishery.discard.all.females" << endl;
        for (int i=1; i<=nObsZCsTCFF; i++) R_out <<  yrsObsZCsTCFF_n(i)  <<  " "  <<  obsPrNatZ_TCFF_nz(i)  << endl;
        R_out << "$Predicted.length.prop.fishery.discard.all.females" << endl;
        for (int i=1; i<=nObsZCsTCFF; i++) {
            ii=yrsObsZCsTCFF_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_TCFF_yz(ii)  << endl;
        }
        R_out << "$Observed.length.prop.snow.fishery.females" << endl;
        for (int i=1; i<=nObsZCsSCF; i++) {
            R_out <<  yrsObsZCsSCF_n(i)  <<  " "  <<  obsPrNatZ_SCF_xnz(1,i)  << endl;
        }
        R_out << "$Predicted.length.prop.snow.fishery.females" << endl;
        for (int i=1; i<=nObsZCsSCF; i++) {
            ii=yrsObsZCsSCF_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_SCF_xyz(1,ii)  << endl;
        }
        R_out << "$Observed.length.prop.snow.fishery.males" << endl;
        for (int i=1; i<=nObsZCsSCF; i++) R_out <<  yrsObsZCsSCF_n(i)  <<  " "  <<  obsPrNatZ_SCF_xnz(2,i)  << endl;
        R_out << "$Predicted.length.prop.snow.fishery.males" << endl;
        for (int i=1; i<=nObsZCsSCF; i++) {
            ii=yrsObsZCsSCF_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_SCF_xyz(2,ii)  << endl;
        }
        R_out << "$Observed.length.prop.redk.fishery.females" << endl;
        for (int i=1; i<=nObsZCsRKF; i++) R_out <<  yrsObsZCsRKF_n(i)  <<  " "  <<  obsPrNatZ_RKF_xnz(1,i)  << endl;
        R_out << "$Predicted.length.prop.redk.fishery.females" << endl;
        for (int i=1; i<=nObsZCsRKF; i++) {
            ii=yrsObsZCsRKF_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_RKF_xyz(1,ii)  << endl;
        }
        R_out << "$Observed.length.prop.redk.fishery.males" << endl;
        for (int i=1; i<=nObsZCsRKF; i++) R_out <<  yrsObsZCsRKF_n(i)  <<  " "  <<  obsPrNatZ_RKF_xnz(2,i)  << endl;
        R_out << "$Predicted.length.prop.redk.fishery.males" << endl;
        for (int i=1; i<=nObsZCsRKF; i++) {
            ii=yrsObsZCsRKF_n(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_RKF_xyz(2,ii)  << endl;
        }
        
        R_out << "$Predicted.length.prop.trawl.females" << endl;
        for (int i=1; i<=nObsZCsGTF; i++) {
            ii=yrsObsZCsGTF(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_GTF_xyz(1,ii)  << endl;
        }
        R_out << "$Observed.length.prop.trawl.females" << endl;
        for (int i=1; i<=nObsZCsGTF; i++) R_out <<  yrsObsZCsGTF(i)  <<  " "  <<  obsPrNatZ_GTF_xnz(1,i)  << endl;
        R_out << "$Predicted.length.prop.trawl.males" << endl;
        for (int i=1; i<=nObsZCsGTF; i++) {
            ii=yrsObsZCsGTF(i);  
            R_out <<  ii  <<  " "  <<  modPrNatZ_GTF_xyz(2,ii)  << endl;
        }
        R_out << "$Observed.length.prop.trawl.males" << endl;
        for (int i=1; i<=nObsZCsGTF; i++) R_out <<  yrsObsZCsGTF(i)  <<  " "  <<  obsPrNatZ_GTF_xnz(2,i)  << endl;
        
        R_out << "$Observed.Length.Prop.survey.immature.new.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);
            R_out << ii <<" " <<obsPrNatZ_Srv_msxnz(1,1,1,i) << endl;
        }
        R_out << "$Predicted.length.prop.survey.immature.new.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_NS_mxyz(1,1,ii) << endl;
        }
        R_out << "$Observed.Length.Prop.survey.immature.old.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);
            R_out << ii <<" " <<obsPrNatZ_Srv_msxnz(1,2,1,i) << endl;
        }
        R_out << "$Predicted.length.prop.survey.immature.old.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_OS_mxyz(1,1,ii) << endl;
        }
        
        R_out << "$Observed.Length.Prop.survey.immature.new.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) R_out << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(1,1,2,i) << endl;
        R_out << "$Predicted.length.prop.survey.immature.new.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_NS_mxyz(1,2,ii) << endl;
        }
        R_out << "$Observed.Length.Prop.survey.immature.old.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) R_out << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(1,2,2,i) << endl;
        R_out << "$Predicted.length.prop.survey.immature.old.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_OS_mxyz(1,2,ii) << endl;
        }
        R_out << "$Observed.Length.Prop.survey.mature.new.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);
            R_out << ii <<" " <<obsPrNatZ_Srv_msxnz(2,1,1,i) << endl;
        }
        R_out << "$Predicted.length.prop.survey.mature.new.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_NS_mxyz(2,1,ii) << endl;
        }
        R_out << "$Observed.Length.Prop.survey.mature.old.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);
            R_out << ii <<" " <<obsPrNatZ_Srv_msxnz(2,2,1,i) << endl;
            }
        R_out << "$Predicted.length.prop.survey.mature.old.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_OS_mxyz(2,1,ii) << endl;
        }
        
        R_out << "$Observed.Length.Prop.survey.mature.new.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) R_out << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(2,1,2,i) << endl;
        R_out << "$Predicted.length.prop.survey.mature.new.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_NS_mxyz(2,2,ii) << endl;
        }
        R_out << "$Observed.Length.Prop.survey.mature.old.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) R_out << yrsObsZCsSrv_n(i) <<" " <<obsPrNatZ_Srv_msxnz(2,2,2,i) << endl;
        R_out << "$Predicted.length.prop.survey.mature.old.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_OS_mxyz(2,2,ii) << endl;
        }
         R_out << "$Observed.Length.Prop.survey.all.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);
            R_out << ii <<" " <<obsPrNatZ_Srv_msxnz(1,1,1,i)+obsPrNatZ_Srv_msxnz(2,1,1,i)+obsPrNatZ_Srv_msxnz(2,2,1,i)<< endl;
                      tmpp4+=obsPrNatZ_Srv_msxnz(1,1,1,i)+obsPrNatZ_Srv_msxnz(2,1,1,i)+obsPrNatZ_Srv_msxnz(2,2,1,i);
        }
        R_out << "$Predicted.length.prop.survey.all.females" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_NS_mxyz(1,1,ii)+modSrvPrNatZ_NS_mxyz(2,1,ii)+modSrvPrNatZ_OS_mxyz(2,1,ii) << endl;
            tmpp1+=modSrvPrNatZ_NS_mxyz(1,1,ii)+modSrvPrNatZ_NS_mxyz(2,1,ii)+modSrvPrNatZ_OS_mxyz(2,1,ii);
        }
        R_out << "$Observed.Length.Prop.survey.all.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);
            R_out << ii <<" " <<obsPrNatZ_Srv_msxnz(1,1,2,i)+obsPrNatZ_Srv_msxnz(1,2,2,i)+obsPrNatZ_Srv_msxnz(2,1,2,i)+obsPrNatZ_Srv_msxnz(2,2,2,i)<< endl;
                 tmpp2+=obsPrNatZ_Srv_msxnz(1,1,2,i)+obsPrNatZ_Srv_msxnz(1,2,2,i)+obsPrNatZ_Srv_msxnz(2,1,2,i)+obsPrNatZ_Srv_msxnz(2,2,2,i);
        }
        R_out << "$Predicted.length.prop.survey.all.males" << endl;
        for (int i=1; i<=nObsZCsSrv; i++) {
            ii=yrsObsZCsSrv_n(i);  
            R_out << ii << " " << modSrvPrNatZ_NS_mxyz(1,2,ii)+modSrvPrNatZ_NS_mxyz(2,2,ii)+modSrvPrNatZ_OS_mxyz(2,2,ii) << endl;
            tmpp3+=modSrvPrNatZ_NS_mxyz(1,2,ii)+modSrvPrNatZ_NS_mxyz(2,2,ii)+modSrvPrNatZ_OS_mxyz(2,2,ii);
        }
        R_out << "$Sum.of.predicted.prop.survey.all.females" << endl;
                R_out <<tmpp1<<endl;
        R_out << "$Sum.of.predicted.prop.survey.all.males" << endl;
                R_out <<tmpp3<<endl;
        R_out << "$Sum.of.Observed.prop.survey.all.females" << endl;
                R_out <<tmpp4<<endl;
        R_out << "$Sum.of.Observed.prop.survey.all.males" << endl;
                R_out <<tmpp2<<endl;
        
        R_out << "$Predicted.mean.postmolt.length.females"<< endl;
        R_out << mean_length(1) << endl;
        R_out << "$Predicted.mean.postmolt.length.males"<< endl;
        R_out << mean_length(2)<<endl; 
        R_out << "$af1" << endl;
        R_out << af1 << endl;
        R_out << "$am1" << endl;
        R_out << am1 << endl;
        R_out << "$bf1" << endl;
        R_out << bf1 << endl;
        R_out << "$bm1" << endl;
        R_out << bm1 << endl;
        R_out << "$Predicted.probability.of.maturing.females"<<endl<<maturity_est(FEMALE)<<endl;
        R_out << "$Predicted.probability.of.maturing.males"  <<endl<<maturity_est(MALE)<<endl;
        R_out << "$molting.probs.female"<<endl<<moltp(FEMALE)<<endl;
        R_out << "$molting.probs.male"  <<endl<<moltp(MALE)<<endl;
        R_out << "$Molting.probability.mature.males"<<endl<<moltp_mat(MALE)<<endl;
        
        R_out << "$observed.TCF.years.retained.catch"         <<endl<<1965<<":"<<endyr-1<<endl;
        R_out << "$observed.retained.catch.biomass"           <<endl<<obsRetCatchBio(1965,endyr-1)<<endl;
        R_out << "$predicted.retained.catch.biomass"          <<endl<<predRetBioMortTCFM_y<<endl;
        R_out << "$predicted.retained.new.catch.biomass"      <<endl<<(catch_male_ret_new*wtm_z)<<endl;
        R_out << "$predicted.retained.old.catch.biomass"      <<endl<<(catch_male_ret_old*wtm_z)<<endl;
        
        R_out << "$observed.TCF.years.discard.catch"             <<endl<<yrsObsDscTCF_n<<endl;
        R_out << "$observed.TCF.male.tot.biomass.mortality"      <<endl<<obsTotBioMortTCFM<<endl;
        R_out << "$predicted.TCF.male.tot.biomass.mortality"     <<endl<<predTotBioMortTCFM_y<<endl;
        R_out << "$predicted.TCF.new.male.tot.biomass.mortality" <<endl<<(catch_lmale_new*wtm_z)(styr,endyr-1)<<endl;
        R_out << "$predicted.TCF.old.male.to.biomass.mortality"  <<endl<<(catch_lmale_old*wtm_z)(styr,endyr-1) <<endl;
        
        R_out << "$observed.TCF.male.discard.mortality.biomass1"     <<endl<<(obsTotBioMortTCFM-obsRetCatchBio(yrsObsDscTCF_n))<<endl;
        R_out << "$observed.TCF.male.discard.mortality.biomass"      <<endl<< obsDscBioMortTCF(  MALE) << endl;
        R_out << "$observed.TCF.female.discard.mortality.biomass"    <<endl<< obsDscBioMortTCF(FEMALE) << endl;
        R_out << "$predicted.TCF.male.discard.mortality.biomass"     <<endl<< predTotBioMortTCFM_y-predRetBioMortTCFM_y<<endl;
        R_out << "$predicted.TCF.female.discard.mortality.biomass"   <<endl<< predDscBioMortTCFF_y << endl;
        
        R_out << "$observed.SCF.years.discard.catch"                 <<endl<<yrsObsDscSCF<<endl;
        R_out << "$observed.SCF.male.discard.mortality.biomass"      <<endl<< obsDscBioMortSCF(  MALE) << endl;
        R_out << "$observed.SCF.female.discard.mortality.biomass"    <<endl<< obsDscBioMortSCF(FEMALE) << endl;
        R_out << "$predicted.SCF.male.discard.mortality.biomass"     <<endl<< predDscBioMortSCF_xy(  MALE) << endl;  
        R_out << "$predicted.SCF.female.discard.mortality.biomass"   <<endl<< predDscBioMortSCF_xy(FEMALE) << endl;
        
        R_out << "$observed.RKF.years.discard.catch"                 <<endl<<yrsObsDscRKF<<endl;
        R_out << "$observed.RKF.male.discard.mortality.biomass"      <<endl<< obsDscBioMortRKF(  MALE) << endl;
        R_out << "$observed.RKF.female.discard.mortality.biomass"    <<endl<< obsDscBioMortRKF(FEMALE) << endl;
        R_out << "$predicted.RKF.male.discard.mortality.biomass"     <<endl<< predDscBioMortRKF_xy(  MALE) << endl;
        R_out << "$predicted.RKF.female.discard.mortality.biomass"   <<endl<< predDscBioMortRKF_xy(FEMALE) << endl;
        
        R_out << "$observed.GTF.years.discard.catch"                 <<endl<< yrsObsDscGTF<<endl;
        R_out << "$observed.GTF.discard.mortality.biomass"           <<endl<< obsDscBioMortGTF<<endl;
        R_out << "$predicted.GTF.discard.mortality.biomass"          <<endl<< predDscBioMortGTF_y<<endl;
        
        R_out << "$predicted.total.male.catch.biomass" << endl;
        R_out<<predTotBioMortTCFM_y+predDscBioMortRKF_xy(  MALE)+predDscBioMortSCF_xy(  MALE)+predDscBioMortGTF_y/2.0<<endl;//WHY /2.0??
        R_out << "$predicted.total.female.catch.biomass" << endl;
        R_out<<predDscBioMortTCFF_y+predDscBioMortRKF_xy(FEMALE)+predDscBioMortSCF_xy(FEMALE)+predDscBioMortGTF_y/2.0<<endl;//why /2.0??
        //  R_out<<"$Estimated total catch div. by male spawing biomass at fishtime"<<endl;
        //  R_out<<elem_div(predTotBioMortTCFM_y(styr,endyr-1)+predDscBioMortRKFM(styr,endyr-1)+predDscBioMortSCFM(styr,endyr-1)+predDscBioMortGTF(styr,endyr-1)/2.0,mspbio_fishtime(styr,endyr-1))<<endl;
        //  R_out << "$estimated retained catch div. by male spawning biomass at fishtime" << endl;
        //  R_out <<elem_div(predRetBioMortTCFM_y,mspbio_fishtime)(styr,endyr-1) << endl;
        R_out << "$estimated.total.catch.divided.by.male.spawning.biomass.at.fishtime" << endl;
        totcatch=predTotBioMortTCFM_y(styr,endyr-1)+predDscBioMortRKF_xy(  MALE)(styr,endyr-1)+predDscBioMortSCF_xy(  MALE)(styr,endyr-1)+predDscBioMortGTF_y(styr,endyr-1)/2.0;//why /2.0??
        for (int i=styr;i<=(endyr-1);i++) R_out <<totcatch(i)/mspbio_fishtime(i) <<" "; R_out << endl;
        R_out << "$estimated.total.catch.of.legal.males.divided.by.legal.males.at.fishtime" << endl;
        for (int i=styr;i<=(endyr-1);i++) R_out <<predRetBioMortTCFM_y(i)/bio_males_gt101(i) <<" "; R_out << endl;
        
        adstring ftype;
        if (optFM==1) ftype="capture"; else ftype="mortality";
        R_out << "$estimated.TCF.fully.selected."<<ftype<<".rate"<<endl;
        R_out << fTCF_xy(MALE)(styr,endyr-1) << endl;
        R_out << "$estimated.SCF.fully.selected."<<ftype<<".rate"<<endl;
        R_out << fSCF_xy(MALE)(styr,endyr-1) << endl;
        R_out << "$estimated.RKF.fully.selected."<<ftype<<".rate"<<endl;
        R_out << fRKF_xy(MALE)(styr,endyr-1) << endl;
        R_out << "$estimated.GTF.fully.selected."<<ftype<<".rate"<<endl;
        R_out << fGTF_xy(MALE)(styr,endyr-1) <<endl;

        //wts: 20150601: fc's are CAPTURE rates (ONLY calculated if using gmacs calculations)
        if (optFM==1){
            //max rates
            R_out<<"$max.TCFM.new.shell.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$max.TCFM.old.shell.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$max.TCF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcTCFF_yz(i)) <<" "; R_out<< endl;
            R_out<<"$max.SCF.male.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcSCF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$max.SCF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$max.RKF.male.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcRKF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$max.RKF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$max.GTF.male.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcGTF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$max.GTF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            //mean rates
            R_out<<"$mean.TCFM.new.shell.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$mean.TCFM.old.shell.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$mean.TCF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcTCFF_yz(i)) <<" "; R_out<< endl;
            R_out<<"$mean.SCF.male.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcSCF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$mean.SCF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$mean.RKF.male.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcRKF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$mean.RKF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$mean.GTF.male.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcGTF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$mean.GTF.female.capture.rate"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        }

        //max fishing MORTALITY RATES (changed f... to fm... to clarify: 20150601)
        R_out<<"$max.TCF.male.new.shell.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$max.TCF.male.old.shell.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$max.TCF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTCFF_yz(i)) <<" "; R_out<< endl;
        R_out<<"$max.SCF.male.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmSCF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$max.SCF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$max.RKF.male.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmRKF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$max.RKF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$max.GTF.male.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmGTF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$max.GTF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$max.TOT.male.NS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(MALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$max.TOT.male.OS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(MALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$max.TOT.female.NS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(FEMALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$max.TOT.female.OS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(FEMALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        //mean fishing MORTALITY RATES (changed f... to fm... to clarify: 20150601)
        R_out<<"$mean.TCF.male.new.shell.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$mean.TCF.male.old.shell.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$mean.TCF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTCFF_yz(i)) <<" "; R_out<< endl;
        R_out<<"$mean.SCF.male.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmSCF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$mean.SCF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$mean.RKF.male.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmRKF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$mean.RKF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$mean.GTF.male.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmGTF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$mean.GTF.female.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$mean.TOT.male.NS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(MALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$mean.TOT.male.OS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(MALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$mean.TOT.female.NS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(FEMALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$mean.TOT.female.OS.mortality.rate"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(FEMALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        
        R_out <<"$max.retained.mortality.rate" << endl;
        for (int i=styr;i<=(endyr-1);i++) R_out <<max(fmTCFR_syz(NEW_SHELL,i))<<" "; R_out<<endl; //same as old shell
        R_out <<"$mean.retained.mortality.rate" << endl;
        for (int i=styr;i<=(endyr-1);i++) R_out <<mean(fmTCFR_syz(NEW_SHELL,i))<<" "; R_out<<endl; //same as old shell
        
        if (optFM==1){
            dvar_vector ratio1(1,nZBs);
            
            //total numbers-at-size captured
            for (int yr=styr;yr<=(endyr-1);yr++) {
                //numbers of males captured in TCF
                ratio1 = elem_prod(elem_div(fcTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(MALE,yr));
                //numbers of females captured in TCF
                ratio1 = elem_prod(elem_div(fmTCFF_yz(yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcTCFF_yz(yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcTCFF_yz(yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcTCFF_yz(yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(FEMALE,yr));

                //numbers of males captured in SCF
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(MALE,yr));
                //numbers of females captured in SCF
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(FEMALE,yr));

                //numbers of males captured in RKF
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(MALE,yr));
                //numbers of females captured in RKF
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(FEMALE,yr));

                //numbers of males captured in GTF
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(MALE,yr));
                //numbers of females captured in GTF
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(FEMALE,yr));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(FEMALE,yr));
            }//yr

            //write out fishery captures for TCF
            R_out<<"$numbers.captured.TCF.male.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iTCF,i,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.TCF.male.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iTCF,i,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.TCF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,MALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.TCF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iTCF,i,MALE,m,s)*wtm_z;///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$numbers.captured.TCF.female.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iTCF,i,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.TCF.female.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iTCF,i,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.TCF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,FEMALE))+sum(cpN_fyxmsz(iTCF,i,FEMALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.TCF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iTCF,i,FEMALE,m,s)*wtf_mz(m);///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;

            //write out fishery captures for SCF
            R_out<<"$numbers.captured.SCF.male.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iSCF,i,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iSCF,i,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.SCF.male.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iSCF,i,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iSCF,i,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.SCF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iSCF,i,MALE))+sum(cpN_fyxmsz(iSCF,i,MALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.SCF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iSCF,i,MALE,m,s)*wtm_z;///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$numbers.captured.SCF.female.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iSCF,i,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iSCF,i,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.SCF.female.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iSCF,i,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iSCF,i,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.SCF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iSCF,i,FEMALE))+sum(cpN_fyxmsz(iSCF,i,FEMALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.SCF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iSCF,i,FEMALE,m,s)*wtf_mz(m);///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;
            
            //write out fishery captures for RKF
            R_out<<"$numbers.captured.RKF.male.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iRKF,i,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iRKF,i,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.RKF.male.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iRKF,i,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iRKF,i,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.RKF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iRKF,i,MALE))+sum(cpN_fyxmsz(iRKF,i,MALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.RKF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iRKF,i,MALE,m,s)*wtm_z;///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$numbers.captured.RKF.female.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iRKF,i,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iRKF,i,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.RKF.female.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iRKF,i,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iRKF,i,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.RKF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iRKF,i,FEMALE))+sum(cpN_fyxmsz(iRKF,i,FEMALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.RKF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iRKF,i,FEMALE,m,s)*wtf_mz(m);///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;

            //write out fishery captures for GTF
            R_out<<"$numbers.captured.GTF.male.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iGTF,i,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iGTF,i,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.GTF.male.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iGTF,i,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iGTF,i,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.GTF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iGTF,i,MALE))+sum(cpN_fyxmsz(iGTF,i,MALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.GTF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iGTF,i,MALE,m,s)*wtm_z;///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$numbers.captured.GTF.female.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iGTF,i,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iGTF,i,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.GTF.female.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iGTF,i,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iGTF,i,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.captured.GTF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iGTF,i,FEMALE))+sum(cpN_fyxmsz(iGTF,i,FEMALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.captured.GTF.female"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += cpN_fyxmsz(iGTF,i,FEMALE,m,s)*wtf_mz(m);///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;
            
        
            //discard rates
            for (int s=NEW_SHELL;s<=OLD_SHELL;s++) {
                for (int i=styr;i<=(endyr-1);i++) fdTCFM_syz(s,i) = elem_prod((1.0-retFcn(s,i)),fcTCFM_syz(s,i));  //discard mortality rate
            }//shell category
            
            //total numbers-at-size discarded
            for (int yr=styr;yr<=(endyr-1);yr++) {
                //numbers of males discarded in TCF
                ratio1 = elem_prod(elem_div(fdTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,natl_inew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fdTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,natl_iold_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fdTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,natl_mnew_fishtime(MALE,yr));
                ratio1 = elem_prod(elem_div(fdTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,natl_mold_fishtime(MALE,yr));
            }
            
            //write out fishery discards for TCF
            R_out<<"$numbers.discarded.TCF.male.new.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iTCF,i,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.discarded.TCF.male.old.shell"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iTCF,i,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$numbers.discarded.TCF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << sum(cpN_fyxmsz(iTCF,i,MALE))+sum(cpN_fyxmsz(iTCF,i,MALE)) <<" "; R_out<< endl;
            R_out<<"$biomass.discarded.TCF.male"<<endl;
            for (int i=styr;i<=(endyr-1);i++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) wt += dsN_fyxmsz(iTCF,i,MALE,m,s)*wtm_z;///dot product on z
                }
                R_out << wt <<" "; 
            }  R_out<< endl;
       }
        
        
        R_out << "$estimated.number.of.recruits.female" << endl;
        R_out << rec_y <<endl;          //was endyr-1
        R_out<< "$estimated.number.of.recruits.male" << endl;
        R_out << rec_y <<endl;          //was endyr-1
        
        R_out<<"$distribution.of.recruits.to.length.bins"<<endl;
        R_out<<rec_len<<endl;
        R_out << "$selectivity.fishery.total.new.males"<< endl;
        R_out << selTCFM(NEW_SHELL) << endl;
        R_out << "$selectivity.fishery.total.old.males"<< endl;
        R_out << selTCFM(OLD_SHELL) << endl;
        R_out << "$selectivity.fishery.ret.new.males"<< endl;
        R_out << selTCFR(NEW_SHELL) << endl;
        R_out << "$selectivity.fishery.ret.old.males"<< endl;
        R_out << selTCFR(OLD_SHELL) << endl;
        //wts 2014-04-25: changed to output all years for retention functions
        R_out <<"$retention.curve.males.new"<< endl;
        R_out <<retFcn(NEW_SHELL)<<endl;
        R_out <<"$retention.curve.males.old"<< endl;
        R_out <<retFcn(OLD_SHELL)<<endl;
        R_out << "$selectivity.discard.females"<< endl;
        R_out <<selTCFF<<endl;
        R_out << "$selectivity.trawl.females"<< endl;
        //  R_out <<selGTF(1)<<endl;
        R_out <<selGTF(1,FEMALE)<<endl;
        R_out <<selGTF(2,FEMALE)<<endl;
        R_out <<selGTF(3,FEMALE)<<endl;  
        R_out << "$selectivity.trawl.males"<< endl;
        //  R_out <<selGTF(2)<<endl;
        R_out <<selGTF(1,MALE)<<endl;
        R_out <<selGTF(2,MALE)<<endl;
        R_out <<selGTF(3,MALE)<<endl;
        R_out << "$selectivity.snow.females"<< endl;
        R_out <<selSCF(1,FEMALE)<<endl;
        R_out <<selSCF(2,FEMALE)<<endl;
        R_out <<selSCF(3,FEMALE)<<endl;  
        R_out << "$selectivity.snow.males"<< endl;
        R_out <<selSCF(1,MALE)<<endl;
        R_out <<selSCF(2,MALE)<<endl;
        R_out <<selSCF(3,MALE)<<endl;
        R_out << "$selectivity.redk.females"<< endl;
        R_out <<selRKF(1,FEMALE)<<endl;
        R_out <<selRKF(2,FEMALE)<<endl;
        R_out <<selRKF(3,FEMALE)<<endl;  
        R_out << "$selectivity.redk.males"<< endl;
        R_out <<selRKF(1,MALE)<<endl;
        R_out <<selRKF(2,MALE)<<endl;
        R_out <<selRKF(3,MALE)<<endl;  
        R_out << "$selectivity.survey.females.1974.to.1981"<< endl;
        R_out << selSrv2(FEMALE) << endl;
        R_out << "$selectivity.survey.males.1974.to.1981"<< endl;
        R_out << selSrv2(MALE) << endl;
        R_out << "$selectivity.survey.females.1982.to.1987"<< endl;
        R_out << selSrv2a(FEMALE) << endl;
        R_out << "$selectivity.survey.males.1982.to.1987"<< endl;
        R_out << selSrv2a(MALE) << endl;
        R_out << "$selectivity.survey.females.1988.to.endyr"<< endl;
        R_out << selSrv3(FEMALE) << endl;
        R_out << "$selectivity.survey.males.1988.to.endyr"<< endl;
        R_out << selSrv3(MALE) << endl;
        R_out << "$Observed.Length.Prop.survey.all.females.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(MATURE,NEW_SHELL,FEMALE) << endl;
        R_out << "$Observed.Length.Prop.survey.all.males.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(MATURE,NEW_SHELL,MALE) << endl;
        R_out << "$Observed.Length.Prop.survey.immature.new.females.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,FEMALE) << endl;
        R_out << "$Observed.Length.Prop.survey.immature.new.males.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,MALE) << endl;
        R_out << "$Observed.Length.Prop.survey.immature.old.males.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(IMMATURE,OLD_SHELL,MALE) << endl;
        R_out << "$Observed.Length.Prop.survey.mature.new.females.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(MATURE,NEW_SHELL,FEMALE) << endl;
        R_out << "$Observed.Length.Prop.survey.mature.new.males.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(MATURE,NEW_SHELL,MALE) << endl;
        R_out << "$Observed.Length.Prop.survey.mature.old.females.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(MATURE,OLD_SHELL,FEMALE) << endl;
        R_out << "$Observed.Length.Prop.survey.mature.old.males.sampsize"<< endl;
        R_out << ssObsZCsSrv_msxn(MATURE,OLD_SHELL,MALE) << endl;
        R_out << "$Observed.Prop.fishery.ret.new.males.sampsize"<< endl;
        R_out << ssRetZCsTCF_sn(NEW_SHELL) << endl;
        R_out << "$Observed.Prop.fishery.ret.old.males.sampsize"<< endl;
        R_out << ssRetZCsTCF_sn(OLD_SHELL) << endl;
        R_out << "$Observed.Prop.fishery.total.new.males.sampsize"<< endl;
        R_out << ssTotZCsTCFM_sn(NEW_SHELL) << endl;  
        R_out << "$Observed.Prop.fishery.total.old.males.sampsize"<< endl;
        R_out << ssTotZCsTCFM_sn(OLD_SHELL) << endl;  
        R_out << "$Observed.Prop.fishery.discard.all.females.sampsize"<< endl;
        R_out << ssZCsTCFF_n << endl;  
        R_out <<"$effectiveN.survey.immature.newshell.female "<<endl;
        R_out <<effn_srv1(IMMATURE,NEW_SHELL,FEMALE)<<endl;
        R_out <<"$effectiveN survey.mature.newshell.female "<<endl;
        R_out <<effn_srv1(MATURE,NEW_SHELL,FEMALE)<<endl;
        R_out <<"$effectiveN.SURVEY.mature.oldshell.female "<<endl;
        R_out <<effn_srv1(MATURE,OLD_SHELL,FEMALE)<<endl;
        R_out <<"$effectiveN.survey.immature.newshell.male "<<endl;
        R_out <<effn_srv1(IMMATURE,NEW_SHELL,MALE)<<endl;
        R_out <<"$effectiveN.survey.immature.oldshell.male "<<endl;
        R_out <<effn_srv1(IMMATURE,OLD_SHELL,MALE)<<endl;
        R_out <<"$effectiveN.survey.mature.newshell.male "<<endl;
        R_out <<effn_srv1(MATURE,NEW_SHELL,MALE)<<endl;
        R_out <<"$effectiveN.survey.mature.old.shell.male "<<endl;
        R_out <<effn_srv1(MATURE,OLD_SHELL,MALE)<<endl;
        R_out <<"$effectiveN.retained.retained.byshell"<<endl;
        R_out <<effn_fish_ret<<endl;
        R_out <<"$effectiveN.total.TCF.byshell.male"<<endl;
        R_out <<effn_fish_tot<<endl;

//    cout<<"done writeToR"<<endl;

// ==========================================================================
FUNCTION void myWriteParamsToR(ostream& os)
    cout<<"starting myWriteParamsToR"<<endl;
    adstring strp;
    os<<"params=list(";
        os<<"growth=list(";
            os<<"af1="<<af1<<cc;
            os<<"bf1="<<bf1<<cc;
            os<<"am1="<<am1<<cc;
            os<<"bm1="<<bm1;
        os<<")"<<cc;
        os<<"natMort=list(";
            os<<"Mmult.imm="<<Mmult_imat<<cc;
            os<<"Mmult.m="<<Mmultm<<cc;
            os<<"Mmult.f="<<Mmultf<<cc;
            os<<"big.mort="; wts::writeToR(os,value(mat_big),dmX); 
        os<<")"<<cc;
        os<<"molting=list(";
            os<<"a="<<moltp_ammat<<cc;
            os<<"b="<<moltp_bmmat; 
        os<<")"<<cc;
        os<<"recruitment=list(";
            strp = "y="+str(mnYrRecCurr)+":"+str(endyr);
            os<<"pMnLnRec="<<pMnLnRec<<cc;
            os<<"pRecDevs="; wts::writeToR(os,value(pRecDevs),strp);  os<<cc;
            strp = "y="+str(styr)+":"+str(mnYrRecCurr-1);
            os<<"pMnLnRecHist="<<pMnLnRecHist<<cc;
            os<<"pRecDevsHist="; wts::writeToR(os,value(pRecDevsHist),strp); 
        os<<")"<<cc;
        os<<"fishery.mortality=list(";
            os<<"tcf=list(";
                strp = "y=c("+wts::to_qcsv(ptrMDS->pTCFR->yrsCatch)+")";
                os<<"pAvgLnFM="<<pAvgLnF_TCF<<cc;
                os<<"pAvgLnFMF="<<pAvgLnF_TCFF<<cc;
                os<<"pFmDevs="; wts::writeToR(os,value(pF_DevsTCF),strp); 
            os<<"),";
            os<<"scf=list(";
                strp = "y=c("+wts::to_qcsv(ptrMDS->pSCF->yrsCatch)+")";
                os<<"pAvgLnFM="<<pAvgLnF_SCF<<cc;
                os<<"pAvgLnFMF="<<pAvgLnF_SCFF<<cc;
                os<<"pFmDevs="; wts::writeToR(os,value(pF_DevsSCF),strp); 
            os<<"),";
            os<<"rkf=list(";
                strp = "y=c("+wts::to_qcsv(ptrMDS->pRKF->yrsCatch)+")";
                os<<"pAvgLnFM="<<pAvgLnF_RKF<<cc;
                os<<"pAvgLnFMF="<<pAvgLnF_RKFF<<cc;
                os<<"pFmDevs="; wts::writeToR(os,value(pF_DevsRKF),strp); 
            os<<"),";
            os<<"gtf=list(";
                strp = "y=c("+wts::to_qcsv(ptrMDS->pGTF->yrsCatch)+")";
                os<<"pAvgLnFM="<<pAvgLnF_GTF<<cc;
                os<<"pAvgLnFMF="<<pAvgLnF_GTFF<<cc;
                os<<"pFmDevs="; wts::writeToR(os,value(pF_DevsGTF),strp); 
            os<<")";
        os<<")"<<cc;
        os<<"fishery.selectivity=list(";
            {
            os<<"tcf=list(";
                dvar_vector sel(1,2); dvar_vector slp(1,2); 
                sel(1) = fish_fit_sel50_mn1; sel(2) = fish_fit_sel50_mn2;
                slp(1) = fish_fit_slope_mn1; slp(2) = fish_fit_slope_mn2;
                os<<"retention=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                os<<"male=list(";
                    if (phase_logistic_sel<0){
                      os<<"=list(z50="<<fish_sel50_1<<cc<<"slope="<<fish_slope_1<<"),";
                    } else {
                        sel(1) = fish_sel50_1; sel(2) = mfexp(log_avg_sel50_3);
                        slp(1) = fish_slope_1; slp(2) = fish_slope_yr_3;
                        strp = qt+str(1)+":"+str(nlog_sel50_dev_3)+qt;
                        os<<"z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<cc<<"devs.lnSel50="; wts::writeToR(os,value(log_sel50_dev_3),strp);
                        if (phase_fishsel) os<<",descending.limb=list(z50="<<fish_sel50_mn2<<cc<<"slope="<<fish_slope_mn2<<")";
                    }
                os<<"),";
                os<<"female=list(z50="<<selTCFF_z50<<cc<<"slope="<<selTCFF_slp<<")";
            os<<"),";
            }
            {
            os<<"scf=list(";
                dvar_vector sel(1,3); dvar_vector slp(1,3);
                sel(1) = selSCFF_z50A1; sel(2) = selSCFF_z50A2; sel(3) = selSCFF_z50A3;
                slp(1) = selSCFF_slpA1; slp(2) = selSCFF_slpA2; slp(3) = selSCFF_slpA3;
                os<<"female=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                os<<"male=list(";
                    sel(1) = selSCFM_z50A1; sel(2) = selSCFM_z50A2; sel(3) = selSCFM_z50A3;
                    slp(1) = selSCFM_slpA1; slp(2) = selSCFM_slpA2; slp(3) = selSCFM_slpA3;
                    os<<"ascending.limb=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                    sel(1) = selSCFM_lnZ50D1; sel(2) = selSCFM_lnZ50D2; sel(3) = selSCFM_lnZ50D3;
                    slp(1) = selSCFM_slpD1; slp(2) = selSCFM_slpD2; slp(3) = selSCFM_slpD3;
                    os<<"descending.limb=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<")";
                os<<")";
            os<<"),"<<endl;
            }
            {
            os<<"rkf=list(";
                dvar_vector sel(1,3); dvar_vector slp(1,3);
                sel(1) = selRKFF_z50A1; sel(2) = selRKFF_z50A2; sel(3) = selRKFF_z50A3;
                slp(1) = selRKFF_slpA1; slp(2) = selRKFF_slpA2; slp(3) = selRKFF_slpA3;
                os<<"female=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                sel(1) = selRKFM_z50A1; sel(2) = selRKFM_z50A2; sel(3) = selRKFM_z50A3;
                slp(1) = selRKFM_slpA1; slp(2) = selRKFM_slpA2; slp(3) = selRKFM_slpA3;
                os<<"male=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<")";
            os<<"),";
            }
            {
            os<<"gtf=list(";
                dvar_vector sel(1,3); dvar_vector slp(1,3);
                sel(1) = selGTFF_z50A1; sel(2) = selGTFF_z50A2; sel(3) = selGTFF_z50A3;
                slp(1) = selGTFF_slpA1; slp(2) = selGTFF_slpA2; slp(3) = selGTFF_slpA3;
                os<<"female=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                sel(1) = selGTFM_z50A1; sel(2) = selGTFM_z50A2; sel(3) = selGTFM_z50A3;
                slp(1) = selGTFM_slpA1; slp(2) = selGTFM_slpA2; slp(3) = selGTFM_slpA3;
                os<<"male=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<")";
            os<<")";
            }
        os<<")";
    os<<")";
    cout<<"finished myWriteParamsToR"<<endl;
    
// ==========================================================================
FUNCTION void myWriteModPopInfoToR(ostream& os)
    cout<<"starting myWriteModPopInfoToR"<<endl;
    
    d5_array nAtZ(1,nSXs,1,nSCs,1,nMSs,styr,endyr,1,nZBs);
    for (int x=1;x<=nSXs;x++){
        for (int y=styr;y<=endyr;y++){
            nAtZ(x,NEW_SHELL,IMMATURE,y) = value(natlength_inew(x,y));
            nAtZ(x,NEW_SHELL,  MATURE,y) = value(natlength_mnew(x,y));
            nAtZ(x,OLD_SHELL,IMMATURE,y) = value(natlength_iold(x,y));
            nAtZ(x,OLD_SHELL,  MATURE,y) = value(natlength_mold(x,y));
        }
    }
    
    os<<"mod.pop=list("<<endl;
        os<<"rec="; wts::writeToR(os,value(rec_y),dmY);                             os<<cc<<endl;
        os<<"MMB="; wts::writeToR(os,value(mspbio_matetime),dmYm1);                 os<<cc<<endl;
        os<<"nAtZ="<<endl; wts::writeToR(os,nAtZ,dmX,dmS,dmM,dmY,dmZ); os<<endl;
    os<<")";
    
    cout<<"finished myWriteModPopInfoToR"<<endl;
    
// ==========================================================================
FUNCTION void myWriteModFshInfoToR(ostream& os)
    cout<<"starting myWriteModFshInfoToR"<<endl;

    dvector zsTMB_TCFM_y(styr,endyr-1);
    dvector zsDMB_TCFF_y(styr,endyr-1);
    dmatrix zsDMB_SCF_xy(1,nSXs,styr,endyr-1);
    dmatrix zsDMB_RKF_xy(1,nSXs,styr,endyr-1);
    dvector zsDMB_GTF_y(styr,endyr-1);
    
    zsTMB_TCFM_y.initialize();
    zsDMB_TCFF_y.initialize();
    zsDMB_SCF_xy.initialize();
    zsDMB_RKF_xy.initialize();
    zsDMB_GTF_y.initialize();
    cout<<"TCF"<<endl;
    for (int n=1;n<=nObsDscTCF;n++){
        int y = yrsObsDscTCF_n(n);
        cout<<n<<tb<<y<<endl;
        if (optTCFMfit==0){
            zsTMB_TCFM_y(y) = value(zsTotMortBio_TCFM_n(n));
        } else {
            zsTMB_TCFM_y(y) = value(zsDscMortBio_TCFM_n(n));
        }
        zsDMB_TCFF_y(y) = value(zsDscMortBio_TCFF_n(n));
    }
    cout<<"SCF"<<endl;
    for (int n=1;n<=nObsDscSCF;n++){
        int y = yrsObsDscSCF(n);
        cout<<n<<tb<<y<<endl;
        for (int x=1;x<=nSXs;x++) zsDMB_SCF_xy(x,y) = value(zsDscMortBio_SCF_xn(x,n));
    }
    cout<<"RKF"<<endl;
    for (int n=1;n<=nObsDscRKF;n++){
        int y = yrsObsDscRKF(n);
        cout<<n<<tb<<y<<endl;
        for (int x=1;x<=nSXs;x++) zsDMB_RKF_xy(x,y) = value(zsDscMortBio_RKF_xn(x,n));
    }
    cout<<"GTF"<<endl;
    for (int n=1;n<=nObsDscGTF;n++){
        int y = yrsObsDscGTF(n);
        cout<<n<<tb<<y<<endl;
        zsDMB_GTF_y(y) = value(zsDscMortBio_GTF_n(n));
    }
    
    os<<"fsh=list("<<endl;
        os<<"TCFR=list("<<endl;
            os<<"fits=list("<<endl;
                os<<"zscr="; wts::writeToR(os,value(zsRetMortBio_TCFR_y),dmYm1); os<<cc<<endl;
                os<<"like="<<lkRetMortBio_TCFR<<cc<<"type='retained'"<<endl;
            os<<")"; os<<endl;
        os<<")"; os<<cc<<endl;
        os<<"TCFM=list("<<endl;
            os<<"fits=list("<<endl;
                os<<"zscr="; wts::writeToR(os,zsTMB_TCFM_y,dmYm1); os<<cc<<endl;
                if (optTCFMfit==0){
                    os<<"like="<<lkTotMortBio_TCFM<<cc<<"type='total'"<<endl;
                } else {
                    os<<"like="<<lkDscMortBio_TCFM<<cc<<"type='bycatch'"<<endl;
                }
            os<<")"; os<<endl;
        os<<")"; os<<cc<<endl;
        os<<"TCFF=list("<<endl;
            os<<"fits=list("<<endl;
                os<<"zscr="; wts::writeToR(os,zsDMB_TCFF_y,dmYm1); os<<cc<<endl;
                os<<"like="<<lkDscMortBio_TCFF<<cc<<"type='bycatch'"<<endl;
            os<<")"; os<<endl;
        os<<")"; os<<cc<<endl;
        os<<"SCF=list("<<endl;
            os<<"fits=list("<<endl;
                os<<"zscr="; wts::writeToR(os,zsDMB_SCF_xy,dmX,dmYm1);        os<<cc<<endl;
                os<<"like="; wts::writeToR(os,value(lkDscMortBio_SCF_x),dmX); os<<cc<<"type='bycatch'"<<endl;
            os<<")"; os<<endl;
        os<<")"; os<<cc<<endl;
        os<<"RKF=list("<<endl;
            os<<"fits=list("<<endl;
                os<<"zscr="; wts::writeToR(os,zsDMB_RKF_xy,dmX,dmYm1);        os<<cc<<endl;
                os<<"like="; wts::writeToR(os,value(lkDscMortBio_RKF_x),dmX); os<<cc<<"type='bycatch'"<<endl;
            os<<")"; os<<endl;
        os<<")"; os<<cc<<endl;
        os<<"GTF=list("<<endl;
            os<<"fits=list("<<endl;
                os<<"zscr="; wts::writeToR(os,zsDMB_GTF_y,dmYm1); os<<cc<<endl;
                os<<"like="<<lkDscMortBio_GTF<<cc<<"type='bycatch'"<<endl;
            os<<")"; os<<endl;
        os<<")"; os<<endl;
    os<<")";
    
    cout<<"finished myWriteModFshInfoToR"<<endl;
    
// ==========================================================================
FUNCTION void myWriteModSrvInfoToR(ostream& os)
    cout<<"starting myWriteModSrvInfoToR"<<endl;

    dmatrix zsc(1,nSXs,styr,endyr);
    zsc.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int n=1;n<=nObsSrvBio;n++){
            int y = yrsObsSrvBio_n(n);
            zsc(x,y) = value(zsSrvMatBio_xn(x,n));
        }
    }
    
    os<<"srv=list("<<endl;
        os<<"fits=list("<<endl;
            os<<"zscr="; wts::writeToR(os,zsc,dmX,dmY); os<<cc<<endl;
            os<<"like="; wts::writeToR(os,value(lkSrvMatBio_x), dmX);     os<<endl;
        os<<")"; os<<endl;
    os<<")";
    
    cout<<"finished myWriteModSrvInfoToR"<<endl;
    
// ==========================================================================
FUNCTION void myWriteToR(ostream& os)
    cout<<"starting myWriteToR"<<endl;
    os<<"res<-list("<<endl;
        ptrMDS->writeToR(os,adstring("model.data")); os<<cc<<endl;
        myWriteParamsToR(os);                        os<<cc<<endl;
        myWriteModPopInfoToR(os);                    os<<cc<<endl;
        myWriteModFshInfoToR(os);                    os<<cc<<endl;
        myWriteModSrvInfoToR(os);                    os<<endl;
    os<<")"<<endl;
    
    cout<<"finished myWriteToR"<<endl;
    
// ==========================================================================
FUNCTION void myWriteMCMCToR(ostream& os)
    cout<<"starting myWriteMCMCToR"<<endl;
    
    d4_array nAtZlast(1,nSXs,1,nSCs,1,nMSs,1,nZBs);
    for (int x=1;x<=nSXs;x++){
        for (int y=styr;y<=endyr;y++){
            nAtZlast(x,NEW_SHELL,IMMATURE) = value(natlength_inew(x,endyr));
            nAtZlast(x,NEW_SHELL,  MATURE) = value(natlength_mnew(x,endyr));
            nAtZlast(x,OLD_SHELL,IMMATURE) = value(natlength_iold(x,endyr));
            nAtZlast(x,OLD_SHELL,  MATURE) = value(natlength_mold(x,endyr));
        }
    }
    
    os<<"list("<<endl;
        myWriteParamsToR(os); os<<cc<<endl;
        os<<"rec="; wts::writeToR(os,value(rec_y),dmY); os<<cc<<endl;
        os<<"MMB="; wts::writeToR(os,value(mspbio_matetime),dmYm1); os<<cc<<endl;
        os<<"nAtZ.last="; wts::writeToR(os,nAtZlast,dmX,dmS,dmM,dmZ); os<<cc<<endl;
        
    os<<"dummy=0),"<<endl;
    
    cout<<"finished myWriteMCMCToR"<<endl;
    
// ==========================================================================
FUNCTION void writeLikelihoodComponents(ostream& os, int toR)
    os<<"----------------------------------------------------------------------------"<<endl;
    os<<"Likelihood components"<<endl;
    os<<"----------------------------------------------------------------------------"<<endl;
    os<<f<<cc<<tb<<"objective function value"<<endl;
    os<<"idx,   weight,      likelihood,      objFun,    description"<<endl;
    for (int i=1;i<=NUM_FOUT;i++){
        os<<i<<cc<<wgtsOut(i)<<cc<<likeOut(i)<<cc<<objfOut(i)<<cc<<tb<<strFOUT(i)<<endl;
    }
  
// ==========================================================================
// ==========================================================================
REPORT_SECTION
    cout<<"starting REPORT_SECTION"<<endl;
    
    if (active(log_sel50_dev_3)) { 
        double llw = 0.0;
        double red = 0.01;
        int max_number_phases = 8;
        int phs = log_sel50_dev_3.get_phase_start();
        llw = 1.0; llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(max_number_phases-phs)))*llw;
        cout<<"llw for log_sel50_dev_3"<<endl;
        cout<<current_phase()<<tb<<phs<<tb<<max_number_phases<<endl;
        cout<<(current_phase()-phs)<<tb<<max(1.0,1.0*(max_number_phases-phs))<<endl;
        cout<<pow(red,(current_phase()-phs)/max(1.0,1.0*(max_number_phases-phs)))<<endl;
     }
    if (last_phase()) {
        ofstream os("TCSAM_WTS.final.R", ios::trunc);    
        myWriteToR(os);
        os.close();

        writeReport(report);

        os.open("TCSAM_OLDSTYLE.R");
        writeToR(os);
        os.close();
        
        os.open("TCProjMod2013.dat");
        writeMyProjectionFile(os);
        os.close();
    
        os.open("TCSAM_WTS.final_params.active.csv");
        writeParameters(os,0,1);
        os.close();
        os.open("TCSAM_WTS.final_params.all.csv");
        writeParameters(os,0,0);
        os.close();
        
        os.open("TCSAM_WTS.final_likelihood_components.csv");
        writeLikelihoodComponents(os,0);
        os.close();
    }
//  report << offset << endl;
    cout<<"finished REPORT_SECTION"<<endl;

// ===============================================================================
// ===============================================================================
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 500,1000,3000,3000,5000,5000,10000
  convergence_criteria 1,1,.01,.001,.001,.001,1e-3,1e-4

// ===============================================================================
// ===============================================================================
TOP_OF_MAIN_SECTION
  arrmblsize = 9000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(4000000); // this may be incorrect in the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(150000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(6000);
  time(&start);
  CheckFile.open("CheckFile.dat");
  

// ===============================================================================
// ===============================================================================
FINAL_SECTION
  
    if (option_match(ad_comm::argc,ad_comm::argv,"-jitter")>-1) {
        ofstream fs("jitterInfo.csv");
        fs<<"seed"<<cc<<"objfun"<<endl;
        fs<<iSeed<<cc<<f<<endl;
    }
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        closeMCMCFile();
    }    

    time(&finish); 
    elapsed_time = difftime(finish,start);
    hour = long(elapsed_time)/3600;
    minute = long(elapsed_time)%3600/60;
    second = (long(elapsed_time)%3600)%60;
    cout << endl << endl << "Starting time: " << ctime(&start);
    cout << "Finishing time: " << ctime(&finish);
    cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
    CheckFile << endl << endl << "Starting time: " << ctime(&start);
    CheckFile<<endl<<"----------------------------------------------------------"<<endl;
    CheckFile << "Finishing time: " << ctime(&finish);
    CheckFile << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
    CheckFile<<"----------------------------------------------------------"<<endl;
    CheckFile.close();


