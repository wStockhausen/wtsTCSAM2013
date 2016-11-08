//Bering Sea Tanner crab model
//********
//--Rec devs indexing changed so that rec_devs(y) enter n-at-size in year y. Previously, rec_devs(y) entered
//      population in year y+1.
//--Changed fishery-related devs so last year is endyr-1, not endyr. This is consistent with year being "survey year". Assessment in 2013
//      uses 2013 survey but last possible fishery year is 2012/13 (so 2012 by survey year).
//--20140419: changed to read model control file name from model configuration, rather than hard-wired name
//--20140425: changed output of selTCFR_syz to retFcn_syz in a number of places
//--20140506: added parameter jitter functionality
//--20140523: updated to use wtsADMB library.
//--20140602: updated to use writeParameter and jitterParameter functions in wtsADMB library.
//--20140805: decreased lower bound on pSelTCFM_devsZ50 (directed fishery log-scale selectivity devs) from -0.5 to -1.0
//--20140814: changed bounds on pSelTCFM_devsZ50 (directed fishery log-scale selectivity devs) from [-1.0,0.5] to [-5,5]
//--20140821: put prior on log_sel50_dev3, implemented decrease with in nll weighting on fishing-related devs 
//--20140821: added flag (doPenRed) to control file to reduce penalties on fishing-related devs with phase. 
//--20140823: added FmRKF phase and pSelTCFM_devsZ50 bounds to control file as inputs, 
//            changed jitter factor on devs to 0.1*fac
//--20140830: corrected penalty reduction for pF_DevsGTF, corrected objfOut for pSelTCFM_devsZ50 and pFmTCF
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
//--20150817: 1. Replaced hard-wired prior mean, variance (0.88, 0.05) used for pSrv2_QM and pSrv2_QF priors
//               with srv3_qPriorMean, Var and srv3_qFemPriorMean, Var, with values read in from the
//               control file.
//            2. Changed q_prior_switch to srv3_qPriorWgt and srv3_qFemPriorWgt to have separate
//               switches for male, female q's. In addition, these also function as likelihood multipliers (weights).
//            3. Changed ...femQ... to ...qFem...
//            4. Now reading phsM from control file to control estimation phase for natural mortality components
//            5. Now reading initial values for pMfac_Imm, pMfac_MatM, pMfac_MatF, pMfac_Big from control file
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
//            3. added zLegal (as legalSize), recLag, mnYrRecDevsHist and mnYrRecCurr to R output
//            4. Modified stock/recruit-related sdreport variables to run styr:(endyr-recLag)
//            5. Added option to reparameterize prob of molt to maturity based on logit-scale
//                  parameters, as well as to set estimation phase, in the control file.
//            6. Renamed M multipliers.
//            7. Added option to ESTIMATE the scaling factor from effort to fishing mortality
//                  used to extrapolate effort to fishing mortality in years where (by)catch data
//                  is unavailable but effort data is.
//--20160401: 1. Fixed array size problems in writeToR_OLD when gmacs fishing mortality option selected.
//            2. Writing 2 sets of initial, final parameter values to csv files as
//                  TCSAM2013.params.XXX.init.csv and TCSAM2013.params.XXX.final.csv,
//                  where XXX = 'active' and 'all'.
//--20160608: 1. Revised output to "old" R file.
//            2. Converted all internal calculations and most output to units of:
//                  MILLIONS for abundance
//                  1000's t for biomass
//                Note that this changes the estimated ln-scale mean recruitment by Log(1000)
//                Units for the Projection Model file remain as before.
//--20160609: 1. Corrected legal abundance/biomass calculations
//            2. Corrected sample size calculations
//--20160614: 1. Removed old report in favor of output to R
//--20160615: 1. Revised recruitment parameters and rec_y to refer to total recruitment, 
//                 not recruitment by sex (new rec_y = 2 * old rec_y).
////            Not done (but should be--encountered estimation problems):
////            2. Revised rec_y and associated sdr variables to run styr to endyr-1,
////                 so rec_y(y) enters population at start of y+1 (as in Jack's original approach).
////                 This will allow direct comparisons with TCSAM2015 when starting population from 0.
//--20160705: 1. Removed all molt parameters no longer involved in any calculations (pPrMoltFA, pPrMoltFB, 
//                  pPrMoltMA, pPrMoltMB, pPrMoltMatMA, pPrMoltMatMB). Modified get_moltingp() to reflect
//                  assumption that only immature crab molt, and they molt annually.
//            2. Removed unused parameters srv2a_q, srv2a_seldiff, srv2a_sel50, 
//                  srv2a_qFem, srv2a_seldiff_f, srv2a_sel50_f, and pPrNewShellRecruits.
//            3. Removed unused parameters fish_slope_mn, log_avg_sel50_mn, log_sel50_dev_mn, 
//                  fish_sel50_1, fish_slope_mn2, fish_sel50_mn2.
//            4. Changed a number of variable names defined in PARAMETER_SECTION to clarify dimensions.
//            5. Changed/removed some variable names in DATA_SECTION. 
//            6. Changed a number of survey-related parameter names to better standardize naming.
//            7. Changed pQFshEff_* to pLnEffXtr_*.
//            8. Changed fishery-related parameter names to better standardize naming.
//--20160706: 1. Changed some control file variable names.
//            2. Updated version to 20160706 and verModelControlFile to 20160622.
//            3. Added CHECK1 macro.
//            4. Reconfigured to read new (20160622) control file format.
//            5. Removed sel_50m_penal and penal_sexr penalties (were always 0).
//--20160708: 1. Switched to using "inp" values from control file to set initial parameter values.
//                  The orig.chk test yielded same final objective function value (2049.13) 
//                  but the max gradient was smaller than previously.
//--20160709: 1. Removed all INITIALIZATION_SECTION values EXCEPT for pPrM2MF, pPrM2MM. 
//                  orig.chk passed as in 20160708.
//--20160711: 1. Changed estimation phases for all parameters to those from control file.
//                  orig.chk passed as in 20160709 (and 08).
//            2. Added comments, moved some definitions of dvar-type quantities (NOT parameters) around.
//                  orig.chk passed as previously.
//            3. Revised order in "writeParameters".
//                  orig.chk passed as previously.
//            4. Converted "exp" to "mfexp" in all instances. <-ROLLED THIS BACK!!
//                  orig.chk resulted in obj fun = 2423.12!
//            5. Converted "exp" to "mfexp" in all instances EXCEPT gamma function calc.s (prGr_xzz, prRec_z).
//                  orig.chk passed as previously.
//            6. Rearranged survey selectivity calculations PRIOR to using useSomertonOtto flags.
//                  orig.chk passed as previously.
//--20160712: 1. Added checks for command line inputs '-ainp' and '-binp' indicating non-standard pin files.
//            2. Implemented useSomertonOtto flags. New approach to orig.chk testing (using
//--20160715: 1. Started to implement Francis method for re-weighting size comps.
//            2. Started revising output names in "oldstyle" R output
//--20160719: 1. Continued revising output names in "oldstyle" R output.
//            2. Combined modSrvPrNatZ_NS_mxyz & modSrvPrNatZ_NS_mxyz as modSrvPrNatZ_msxyz.
//            3. Renamed obsPrNatZ_Srv_msxnz to obsSrvPrNatZ_msxnz.
//            4. Revised effSS calc.s for survey size comps (again) after realizing
//                  normalization of obs(mod)SrvPrNatZ_msxn(y)z is over m,s,x,z, although
//                  the fit combines Pr's over shell condition. Because of the normalization,
//                  only an annual sample size is appropriate (not individual ones for m,x components).
//            5. Added output of Pearson's residuals to "oldstyle" R output to simplify R side.
//            6. Removed lkSrvZCs_msx and calc.s (did not contribute to objfun and was not output).
//--20160720: 1. Continued revising output names in "oldstyle" R output.
//            2. Eliminated OLD_SHELL output for selTCFM_syz, selrTCFM_syz, retTCFM_syz
//                  as redundant.
//--20160726: 1. Added description string to writeParameter(...) functions in writePrameters(...).
//            2. Updated model version to 20160726.
//            3. Renamed pMnLnRecHist and pRecDevsHist as pMnLnRecInit and pRecDevsInit.
//--20160727: 1. Revised sd_report variables (sdr...) to reduce duplication. Added sdrMeanGrowthF/M output.
//--20160730: 1. Changed csv file for objective function components to "TCSAM2013.final_likelihood_components.csv".
//--20160731: 1. Added "category" column to csv file for objective function components.
//            2. Added AR1 penalty on pSelTCFM_devsZ50, with weight set to llwSelTCFM_devsZ50. 
//                  There had been an AR1 penalty on 'fish_sel50_mn' in old code, but it was never
//                  calculated because estimation for fish_sel50_mn was never turned on in recent assessments. 
//--20160814: 1. Added pAvgLn_XXXF parameters, ln-scale female offsets to male fishing mortality, to projection model file
//            2. Changed output projection model file from "TCProjMod2013.dat" to "TCSAM2013ProjMod.dat".
//--20160820: 1. Added 20160820 as a possible model control version (verModelControlFile).
//            2. Added ability to switch off minimum F's using control file flag optMinFs.
//            3. Incremented model version to 20160820.
//            4. Changed zLegal to 125 mm CL.
//--20160825: 1. Revised penalty reductions for F devs so llw = 0 in last phase.
//            2. Revised parameterization and options for effort extrapolation.
//--20160828: 1. Revised penalty reductions for F devs so llw = 0 in phase 8 (can be modified setting -maxph to something else).
//            2. Revised (again) implementation for effort extrapolation.
//            3. Incremented model version to 20160828.
//            4. Pre-data F levels are (again) turned on, regardless of optMinFs, which
//                  now works only for RKC Fs where minF check is made.
//            5. Revised sdrNatMort_ output to run styr to endyr-1.
//--20160905: 1. Writing initial likelihood components to csv file "TCSAM2013.init_likelihood_components.csv".
//            2. Revised PARAMETER_SECTION to better handle modPrM2M when pin file is used.
//            3. Incremented model version to 20160905.
//            4. Corrected labels in writeParameters for male survey selectivity parameters.
//--20160907: 1. writing retFcn_syz to projection model file if optFM=1 instead of selTCFR_syz.
//--20161011: 1. writing survey, fishery numbers/biomass-at-size arrays to rep file
//            2. Incremented version to 20161011.
//--20161013: 1. writing pop numbers/biomass-at-size arrays to rep file
//--20161107: 1. Added internal OFL calculations via OFL_Calcs.hpp/cpp.
//--20161108: 1. Corrected improper use of spmo in OFL calculations.
//
//IMPORTANT: 2013-09 assessment model had RKC params for 1992+ discard mortality TURNED OFF. 
//           THE ESTIMATION PHASE FOR RKC DISCARD MORTALITY IS NOW SET IN THE CONTROLLER FILE!
//
//  Output units (except where explicitly noted):
//                  MILLIONS for abundance
//                  1000's t for biomass
//
//********
//to run mcmc 
//tcsam2013alta -nox -mcmc 1000000 -mcsave 200
// then have to run  tcsam2013alta -mceval to get output
//whatever is in sd report file will have a distribution and output will go to eval.csv
//for whatever have written to post later in program in the mcmc function part
//
// ===============================================================================
// ===============================================================================
GLOBALS_SECTION
    #include <math.h>
    #include <time.h>
    #include <fenv.h> // must appear before admodel.h
    #include <admodel.h>
    #include "wtsADMB.hpp" 
    #include "ModelConstants.hpp"
    #include "ModelConfiguration.hpp"
    #include "FisheryData.hpp"
    #include "ModelData.hpp"
    #include "OFLCalcs.hpp"
    
    adstring version = "20161011";//model version
    ivector verModelControlFile(1,2);  //model control file version
    
    int maxPhase = 8;//default max phase for penalty reduction
    
    double zLegal = 125;//current (2015/16) legal size
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
    
    int recLag = 5;       //default lag from fertilization to recruitment (yrs)
    
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
    
    int doOFL    = 0;
    int debugOFL = 0;
    OFLResults oflResults;
    
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
    
    #undef REP2R1
    #define REP2R1(o) R_out<<"$" #o <<endl<<o<<endl

    #undef REP2R2
    #define REP2R2(o1,o2) R_out<<"$" #o1 <<endl<<o2<<endl

    #undef REP2RTS
    #define REP2RTS(o1,o2,o3) R_out<<"$" #o1 <<endl; for (int i=o3.indexmin();i<=o3.indexmax();i++) {R_out<<o2(o3(i))<<endl;}

    #undef CHECK1
    #define CHECK1(o) CheckFile<<o<<tb<<"#" #o<<endl

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
    verModelControlFile[1] = 20160622;
    verModelControlFile[2] = 20160820;
    
    int kf = 1;
    strFOUT(kf++) = "penalty, recruitment penalty";
    strFOUT(kf++) = "penalty, sex ratio penalty";
    strFOUT(kf++) = "penalty, natural mortality penalty (immatures)";
    strFOUT(kf++) = "penalty, natural mortality penalty (mature males)";
    strFOUT(kf++) = "penalty, natural mortality penalty (immature females)";

    strFOUT(kf++) = "priors, survey q penalty";
    strFOUT(kf++) = "priors, female survey q penalty";
    strFOUT(kf++) = "priors, female growth parameter a";
    strFOUT(kf++) = "priors, female growth parameter b";
    strFOUT(kf++) = "priors, male growth parameter a";
    strFOUT(kf++) = "priors, male growth parameter b";
    strFOUT(kf++) = "penalty, maturity curve smoothness (females)";
    strFOUT(kf++) = "penalty, maturity curve smoothness (males)";
    strFOUT(kf++) = "penalty, z50 devs for male selectivity in TCF (AR1)";
    strFOUT(kf++) = "penalty, penalty on F-devs in directed fishery";
    strFOUT(kf++) = "penalty, penalty on F-devs in snow crab fishery";
    strFOUT(kf++) = "penalty, penalty on F-devs in BBRKC fishery";
    strFOUT(kf++) = "penalty, penalty on F-devs in groundfish fishery";

    strFOUT(kf++) = "likelihood: size comps, fishery: TCF retained males";
    strFOUT(kf++) = "likelihood: size comps, fishery: TCF total males";
    strFOUT(kf++) = "likelihood: size comps, fishery: TCF discarded females";
    strFOUT(kf++) = "likelihood: size comps, fishery: SCF males";
    strFOUT(kf++) = "likelihood: size comps, fishery: SCF females";
    strFOUT(kf++) = "likelihood: size comps, fishery: RKC males";
    strFOUT(kf++) = "likelihood: size comps, fishery: RKC females";
    strFOUT(kf++) = "likelihood: size comps, fishery: GTF males+females";
    strFOUT(kf++) = "likelihood: size comps, survey: immature males";
    strFOUT(kf++) = "likelihood: size comps, survey: mature males";
    strFOUT(kf++) = "likelihood: size comps, survey: immature females";
    strFOUT(kf++) = "likelihood: size comps, survey: mature females";

    strFOUT(kf++) = "likelihood: catch biomass, survey: mature crab";
    strFOUT(kf++) = "likelihood: catch biomass, fishery: TCF retained males";
    strFOUT(kf++) = "likelihood: catch biomass, fishery: TCF male total catch biomass";
    strFOUT(kf++) = "likelihood: catch biomass, fishery: TCF female catch biomass";
    strFOUT(kf++) = "likelihood: catch biomass, fishery: SCF total catch biomass";
    strFOUT(kf++) = "likelihood: catch biomass, fishery: RKF total catch biomass";
    strFOUT(kf++) = "likelihood: catch biomass, fishery: GTF total catch biomass";

    strFOUT(kf++) = "penalty, z50 devs for male selectivity in TCF (norm2)";
 END_CALCS
 
 LOCAL_CALCS  
    //process command line options
    int on = 0;
    int flg = 0;
    //pin file use
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin=1;
        echo<<"#using default pin file"<<endl;
        flg = 1;
    }
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ainp"))>-1) {
        usePin=1;
        echo<<"#using pin file "<<ad_comm::argv[on+1]<<endl;
        flg = 1;
    }
   if ((on=option_match(ad_comm::argc,ad_comm::argv,"-binp"))>-1) {
        usePin=1;
        echo<<"#using pin file "<<ad_comm::argv[on+1]<<endl;
        flg = 1;
    }
    //recruitment lag
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-lag"))>-1) {
        recLag=atoi(ad_comm::argv[on+1]);
        echo<<"#assumed lag for recruitment changed to: "<<recLag<<endl;
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
    //maxPhase
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-imaxph"))>-1) {
        maxPhase=atoi(ad_comm::argv[on+1]);
        rng.reinitialize(iSeed);
        flg = 1;
    }
    echo<<maxPhase<<"  #maxPhase"<<endl;
    //calcOFL
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcOFL"))>-1) {
        doOFL=1;
        echo<<"#OFL calculations turned ON"<<endl;
        echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugOFL
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugOFL"))>-1) {
        debugOFL=1;
        echo<<"#debugOFL turned ON"<<endl;
        echo<<"#-------------------------------------------"<<endl;
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
    
    asmtYr  = ptrMC->asmtYr;//assessment year
    mnYr    = ptrMC->mnYr;  //start year of model
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
    int nObsZCsSCF                                // number of years of snow crab fishery female bycatch length data
    !!nObsZCsSCF = ptrMDS->pSCF->nyNatZ;
    ivector yrsObsZCsSCF_n(1,nObsZCsSCF)          // years which have snow crab fishery bycatch length data
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
    int nObsZCsRKF                                // number of years of BBRKC fishery bycatch size comps
    !!nObsZCsRKF = ptrMDS->pRKF->nyNatZ;
    ivector yrsObsZCsRKF_n(1,nObsZCsRKF)          // years which have BBRKC fishery bycatch size comps
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
    int nObsDscSCF                        // number of years of male and female bycatch weight data
    !!nObsDscSCF = ptrMDS->pSCF->nyCatch;
    ivector yrsObsDscSCF(1,nObsDscSCF)    // years which have male and female bycatch weight data
 LOCAL_CALCS
    CheckFile<<"nObsDscSCF = "<<nObsDscSCF<<endl;
    yrsObsDscSCF = ptrMDS->pSCF->yrsCatch;
    CheckFile<<"yrsObsDscSCF "<<yrsObsDscSCF<<endl;
 END_CALCS
    
    //BBRKC bycatch weight
    int nObsDscRKF                          // number of years of male and female bycatch weight data
    !!nObsDscRKF = ptrMDS->pRKF->nyCatch;
    ivector yrsObsDscRKF(1,nObsDscRKF)      // years which have male and female bycatch weight data
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
    int nObsZCsGTF                          // number of years of trawl bycatch length comps
    !!nObsZCsGTF = ptrMDS->pGTF->nyNatZ;
    ivector yrsObsZCsGTF_n(1,nObsZCsGTF)    // years which have trawl bycatch length data
    matrix ssZCsGTF_xn(1,nSXs,1,nObsZCsGTF) // sample sizes for trawl bycatch length comps (by year and sex)
    vector ssObsZCsGTF_n(1,nObsZCsGTF)      // combined sample sizes for trawl bycatch length comps (by year)
 LOCAL_CALCS
    CheckFile<<"nObsZCsGTF     = "<<nObsZCsGTF<<endl;
    yrsObsZCsGTF_n = ptrMDS->pGTF->yrsNatZ;
    for (int x=1;x<=nSXs;x++) ssZCsGTF_xn(x) = ptrMDS->pGTF->ssNatZ_xsy(x,ALL_SHELL);
    ssObsZCsGTF_n.initialize();
    for (int n=1;n<=nObsZCsGTF;n++) {
        for (int x=1;x<=nSXs;x++) ssObsZCsGTF_n(n) += ssZCsGTF_xn(x,n);
    }
    CheckFile<<"yrsObsZCsGTF_n = "<<yrsObsZCsGTF_n<<endl;
    CheckFile<<"ssZCsGTF_xn  = "<<endl<<ssZCsGTF_xn<<endl;
    CheckFile<<"ssObsZcsGTF_x   = "<<endl<<ssObsZCsGTF_n<<endl;
 END_CALCS   
    
    //survey aggregate data
    int nObsSrvBio                           // number of years of survey biomass data
    !!nObsSrvBio = ptrMDS->pTSD->nyAbund;
    ivector yrsObsSrvBio_n(1,nObsSrvBio)     // years which have survey biomass estimates
    vector obsSrvNum_n(1,nObsSrvBio)         // survey numbers (total) in millions of crab
    matrix obsSrvCV_xn(1,nSXs,1,nObsSrvBio)  // survey cv by sex x year
 LOCAL_CALCS
    obsSrvNum_n = ptrMDS->pTSD->abund_y;
    CheckFile<<"obsSrvNum_n = "<<endl<<tb<<obsSrvNum_n<<endl;
    obsSrvCV_xn = ptrMDS->pTSD->cvsAbund_xy;
    CheckFile<<"obsSrvCV_xn = "<<endl<<obsSrvCV_xn<<endl;
 END_CALCS   
            
    //survey numbers-at-size data
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
    vector obsTotSrvNum_n(1,nObsZCsSrv);                                ///<total survey abundance from survey size comps (normalization factor for probabilities)
    5darray obsSrvNatZs_msxnz(1,nMSs,1,nSCs,1,nSXs,1,nObsZCsSrv,1,nZBs);///<survey size comps (numbers at size by xms)
 LOCAL_CALCS
    obsTotSrvNum_n.initialize();
    obsSrvNatZs_msxnz.initialize();
    for (int m=1;m<=nMSs;m++) {
        for (int s=1;s<=nSCs;s++) {
            for (int x=1;x<=nSXs;x++) {
                for (int n=1;n<=nObsZCsSrv;n++) {
                    obsSrvNatZs_msxnz(m,s,x,n) = ptrMDS->pTSD->nAtZ_xsmyz(x,s,m,n);
                    obsTotSrvNum_n(n) += sum(obsSrvNatZs_msxnz(m,s,x,n));
                }
                CheckFile<<"obsSrvNatZs_msxnz("<<m<<cc<<s<<cc<<x<<") = "<<endl<<obsSrvNatZs_msxnz(m,s,x)<<endl;
            }//x
        }//s
    }//m
    CheckFile<<"obsTotSrvNum_n (millions)"<<endl<<obsTotSrvNum_n<<endl;
    CheckFile<<"obsSrvNum_n (millions)"<<endl<<tb<<obsSrvNum_n<<endl;
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
    int nSelTCFM_devsZ50                      //number of years of directed fishery post 1990
    ivector hasDirectedFishery_y(styr,endyr-1); //flags indicating if directed fishery is prosecuted (>0)
    vector obsRetCatchNum_y(styr,endyr-1)       //retained catch, numbers                 (IMPORTANT CHANGE: used to be "1965,endyr")
    vector obsRetCatchBio_y(styr,endyr-1)       //retained catch, millions of lbs of crab (IMPORTANT CHANGE: used to be "1965,endyr")
 LOCAL_CALCS
    nYrsTCF = 0;
    nSelTCFM_devsZ50    = 0;
    hasDirectedFishery_y  = 0;//set all years to "no directed fishery"
    obsRetCatchNum_y.initialize();
    obsRetCatchBio_y.initialize();
    for (int i=1;i<=ptrMDS->pTCFR->nyCatch;i++) {
        int y = ptrMDS->pTCFR->yrsCatch(i);
        if ((styr<=y)&&(y<endyr)){
            nYrsTCF++;
            hasDirectedFishery_y(y) = 1;
            obsRetCatchNum_y(y)     = ptrMDS->pTCFR->catch_ty(1,i);
            obsRetCatchBio_y(y)     = ptrMDS->pTCFR->catch_ty(2,i);
        }
        if (y>1990) nSelTCFM_devsZ50++;
    }
    CheckFile<<"number of directed fishery years after 1990: "<<nSelTCFM_devsZ50<<endl;
    CheckFile<<"retained numbers:        obsRetCatchNum_y"<<endl<<obsRetCatchNum_y<<endl;
    CheckFile<<"retained biomass (mlbs): obsRetCatchBio_y"<<endl<<obsRetCatchBio_y<<endl;
    obsRetCatchBio_y /= 2.2045; // convert from millions lbs to 1000's tons   
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
    
    3darray wt_xmz(1,nSXs,1,nMSs,1,nZBs)  // weight at length (from kodiak program) in kg (??)
 LOCAL_CALCS
    wt_xmz = ptrMDS->pBio->wAtZ_xmz;
    CheckFile<<"wt_xmz(FEMALE,IMMATURE)= "<<endl<<tb<<wt_xmz(FEMALE,IMMATURE)<<endl;
    CheckFile<<"wt_xmz(MALE,  MATURE)  = "<<endl<<tb<<wt_xmz(FEMALE,  MATURE)<<endl;
    CheckFile<<"wt_xmz(MALE,IMMATURE)  = "<<endl<<tb<<wt_xmz(  MALE,IMMATURE)<<endl;
    CheckFile<<"wt_xmz(MALE,  MATURE)  = "<<endl<<tb<<wt_xmz(  MALE,  MATURE)<<endl;
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
      
    vector zBs(1,nZBs)             // Midpoints of length bins
    !!zBs = ptrMDS->pBio->zBins-0.5;//IMPORTANT: subtract 0.5 to match up with old numbers
    !!CheckFile <<"zBs"<<endl<<tb<<zBs<<endl;
    
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
    bool tst = false;
    for (int i=1;i<=verModelControlFile.indexmax();i++) tst = tst|(inpVerMCF==verModelControlFile[i]);
    if (!tst){
        cout<<"Model Control File version inconsistent with model."<<endl;
        cout<<"Current versions = "<<verModelControlFile<<endl;
        cout<<"Version in file = "<<inpVerMCF<<endl;
        cout<<"Please use one of the current versions."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
 END_CALCS
    //fixed values
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    !!CheckFile<<"#fixed values"<<endl;
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    !!CheckFile<<"##--handling mortality"<<endl;
    init_number hm_pot      // fraction of pot discards that die (.5)
    init_number hm_trawl    // fraction of trawl discards that die(.8)
    !!CHECK1(hm_pot);
    !!CHECK1(hm_trawl);
    
    !!CheckFile<<"##--surveys"<<endl;
    init_number multQ                 // Q  mult by pop biomass to get survey biomass
    init_vector sel_som(1,5)          // parameters for somerton-otto selectivity curve
    !!CHECK1(multQ);
    !!CHECK1(sel_som);
         
    !!CheckFile<<"##--natural mortality base values"<<endl;
    init_vector baseM_Imm(1,nSXs)       // base value for natural mortality on immature crab by sex
    init_vector baseM_MatNS(1,nSXs)     // base value for natural mortality on mature new shell crab by sex
    init_vector baseM_MatOS(1,nSXs)     // base value for natural mortality on mature old shell crab by sex
    3darray baseM_msx(1,nMSs,1,nSCs,1,nSXs) //convenience array for all base M's
    3darray inpMfac_msx(1,nMSs,1,nSCs,1,nSXs) //convenience array for all initial M multipliers
 LOCAL_CALCS
    baseM_msx(IMMATURE,NEW_SHELL,FEMALE) = baseM_Imm(FEMALE);
    baseM_msx(IMMATURE,NEW_SHELL,  MALE) = baseM_Imm(  MALE);
    baseM_msx(IMMATURE,OLD_SHELL,FEMALE) = baseM_Imm(FEMALE);
    baseM_msx(IMMATURE,OLD_SHELL,  MALE) = baseM_Imm(  MALE);
    baseM_msx(  MATURE,NEW_SHELL,FEMALE) = baseM_MatNS(FEMALE);
    baseM_msx(  MATURE,NEW_SHELL,  MALE) = baseM_MatNS(  MALE);
    baseM_msx(  MATURE,OLD_SHELL,FEMALE) = baseM_MatOS(FEMALE);
    baseM_msx(  MATURE,OLD_SHELL,  MALE) = baseM_MatOS(  MALE);
    CHECK1(baseM_Imm);
    CHECK1(baseM_MatNS);
    CHECK1(baseM_MatOS);
 END_CALCS
            
    !!CheckFile<<"##--miscellaneous"<<endl;
    init_number mate_ratio            // mating ratio (USED)
    init_number fraction_new_error    // accounts for shell error
    init_int nages                    // number of ages to track for mature old shell 
    !!CHECK1(mate_ratio);
    !!CHECK1(fraction_new_error);
    !!CHECK1(nages);
    
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    !!CheckFile<<"#likeilihood weights"<<endl;
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    init_vector like_wght(1,NUM_LEN_LIKE) // likelihood weights for size composition data
    init_number like_wght_CatchBio        // likelihood weight for fishery catch biomass fits [was like_wght_(6)]
    init_number like_wght_fbio            // likelihood weight for female biomass fit [was like_wght(5)]
    init_number like_wght_mbio            // likelihood weight for male biomass fit
    init_number like_wght_rec         // ??
    init_number llwSelTCFM_devsZ50 //likelihood weight for penalty on pSelTCFM_devsZ50 (originally 0))
    init_number bndSelTCFM_devsZ50 //upper/lower bounds on pSelTCFM_devsZ50 deviations (originally 0.5)
 LOCAL_CALCS
    CheckFile<<"##--Size composition weights"<<endl;
    CHECK1(like_wght);
    CheckFile<<"##--aggregated abundance/biomass weights"<<endl;
    CHECK1(like_wght_CatchBio);
    CHECK1(like_wght_fbio);
    CHECK1(like_wght_mbio);
    CheckFile<<"##--penalty weights"<<endl;
    CHECK1(like_wght_rec);
    CHECK1(llwSelTCFM_devsZ50);    
    CHECK1(bndSelTCFM_devsZ50);  
 END_CALCS
    
    !!CheckFile<<"##--priors"<<endl;
    init_number srv3_qPriorWgt        //likelihood multiplier for survey q prior (off=0)
    init_number srv3_qPriorMean       //prior mean
    init_number srv3_qPriorStD        //prior variance
    init_number srv3_qFemPriorWgt     //likelihood multiplier for survey female q prior (off=0)
    init_number srv3_qFemPriorMean    //prior mean
    init_number srv3_qFemPriorStD     //prior variance
    !!CHECK1(srv3_qPriorWgt);
    !!CHECK1(srv3_qPriorMean);
    !!CHECK1(srv3_qPriorStD);
    !!CHECK1(srv3_qFemPriorWgt);
    !!CHECK1(srv3_qFemPriorMean);
    !!CHECK1(srv3_qFemPriorStD);
    
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    !!CheckFile<<"#options"<<endl;
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    //population process-related
    !!CheckFile<<"##--options for estimating recruitment"<<endl;
    init_int mnYrRecDevsHist; //year to start estimating "historic" rec devs
    init_int mnYrRecCurr;     //year to start estimating "current" recruitment
    !!if (mnYrRecDevsHist<=0) mnYrRecDevsHist = mnYr;
    !!CHECK1(mnYrRecDevsHist);
    !!CHECK1(mnYrRecCurr);
    !!CheckFile<<"##--options for estimating natural mortality"<<endl;
    init_int mort_switch              // extra mort on 1, off 0
    init_int mort_switch2             // apply mort_switch correctly (1), or as in 2012 (0)
    init_int lyr_mort                 // start yr extra mort
    init_int uyr_mort                 //end yr extra mort
    !!CHECK1(mort_switch);
    !!CHECK1(mort_switch2);
    !!CHECK1(lyr_mort);
    !!CHECK1(uyr_mort);
    !!CheckFile<<"##--options for using maturity info and estimating prM2M"<<endl;
    init_int optInitM2M  ///<option to initialize prM2M(MALE)
    !!CHECK1(optInitM2M);
    init_int optMatForSrvOSM  ///<option to use maturity ogive to determine maturity status
    !!CHECK1(optMatForSrvOSM);
    init_int optPrM2M  ///< option for parameterizing prM2M, the probability of molt to maturity
    number lbPrM2M     ///< lower bound on PrM2M parameters (if optPrM2M==0)
    number ubPrM2M     ///< lower bound on PrM2M parameters (if optPrM2M==0)
 LOCAL_CALCS
    CHECK1(optPrM2M);
    if (optPrM2M==0) {
        CheckFile<<"###--ln-scale parameterization selected for pPrM2M's"<<endl;
        lbPrM2M = -15.0;
        ubPrM2M =   0.0;
    }
    if (optPrM2M==1){
        CheckFile<<"###--logit-scale parameterization selected for pPrM2M's"<<endl;
        lbPrM2M = -15.0;
        ubPrM2M =  15.0;
    }
 END_CALCS
    //selectivity-related options
    !!CheckFile<<"#---options for survey selectivity"<<endl;
    init_int useSomertonOtto1  //use Somerton & Otto selectivity curve for survey period 1
    init_int useSomertonOtto2  //use Somerton & Otto selectivity curve for survey period 2
    init_int useSomertonOtto3  //use Somerton & Otto selectivity curve for survey period 3
 LOCAL_CALCS
    CHECK1(useSomertonOtto1);
    CHECK1(useSomertonOtto2);
    CHECK1(useSomertonOtto3);
 END_CALCS
    !!CheckFile<<"##--asymptotic selectivity options"<<endl;
    init_int optFshSel;//option to force asymptote=1 for logistic fishery selectivity functions
    init_int optSrvSel;//option to force asymptote=1 for logistic survey selectivity functions
    !!if (optFshSel!=1) optFshSel = 0;//if not 1 (on), use old approach to normalization
    !!if (optSrvSel!=1) optSrvSel = 0;//if not 1 (on), use old approach to normalization
    !!CHECK1(optFshSel);
    !!CHECK1(optSrvSel);
    //fisheries-related options
    !!CheckFile<<"##--fishing mortality model options"<<endl;
    init_int optFM;//original (0) or gmacs (1) -like fishing mortality
    !!if (optFM!=1) optFM = 0;//if gmacs not specified, then use original
    !!CHECK1(optFM);
    !!CheckFile<<"##--fishery catch likelihood options"<<endl;
    init_int optFshNLLs  //flag indicating error model for fishery catch data (0=norm2, 1=lognormal)
    !!CHECK1(optFshNLLs);    
    //the following obsErr are sd for normal error, cv for lognormal error
    init_number obsErrTCFR  //assumed observation error level for retained catch data
    init_number obsErrTCFD  //assumed observation error level for discard catch data in TCF
    init_number obsErrSCF   //assumed observation error level for discard catch data in SCF
    init_number obsErrRKF   //assumed observation error level for discard catch data in RKF
    init_number obsErrGTF   //assumed observation error level for discard catch data in GTF
    !!CHECK1(obsErrTCFR);   
    !!CHECK1(obsErrTCFD);   
    !!CHECK1(obsErrSCF);   
    !!CHECK1(obsErrRKF);   
    !!CHECK1(obsErrGTF);   
    !!CheckFile<<"##--options for fitting GTF size comps"<<endl;
    init_int optPrNatZ_GTF; //
    !!CHECK1(optPrNatZ_GTF);
    //new 20160323: options for fitting male mortality in TCF
    !!CheckFile<<"##--options for fitting male mortality in TCF"<<endl;
    init_int optTCFMfit;
    !!CHECK1(optTCFMfit);
    !!CheckFile<<"##--options for penalty reduction on F devs"<<endl;
    init_int doPenRed      //flag (0/1) to reduce penalties on fishing-related devs by phase
    !!CHECK1(doPenRed);    
    !!CheckFile<<"##--options for minimum F's"<<endl;
    int optMinFs;
    !!optMinFs = 1;//min F's used (old style)
    !!if (inpVerMCF>=20160820) (*(ad_comm::global_datafile))>>optMinFs;
    !!if (optMinFs>0) {optMinFs = 1;} else {optMinFs = 0;}//make sure this is 1 if turned on, 0 if off
    !!CHECK1(optMinFs);
    !!CheckFile<<"##--options for effort extrapolation"<<endl;
    init_int optEffXtr_TCF
    init_int optEffXtr_SCF
    init_int optEffXtr_RKF
    init_int optEffXtr_GTF
LOCAL_CALCS
    optEffXtr_TCF = 0; //option not implemented yet
    optEffXtr_GTF = 0; //option not implemented yet
    optEffXtr_SCF = 1; //only option implemented for SCF is 1
    if (optEffXtr_RKF==0) optEffXtr_RKF = 2;//old style is now 2 
    CHECK1(optEffXtr_TCF);
    CHECK1(optEffXtr_SCF);
    CHECK1(optEffXtr_RKF);
    CHECK1(optEffXtr_GTF);
 END_CALCS
         
        
    //parameter estimation phases
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    !!CheckFile<<"#parameter phases"<<endl;
    !!CheckFile<<"#---------------------------------------------------"<<endl;
    !!CheckFile<<"#-----recruitment"<<endl;
    init_int phsMnLnRec      // mean ln-scale total recruitment
    init_int phsRecDevs      // recruitment devs
    !!CHECK1(phsMnLnRec);    
    !!CHECK1(phsRecDevs);
    //--natural mortality
    !!CheckFile<<"#-----natural mortality"<<endl;
    init_int phsM            // ordinary mortality
    init_int phsBigM         // multipliers for high mortality period (1980-1984)                  
    !!CHECK1(phsM);    
    !!CHECK1(phsBigM);
    //--growth
    !!CheckFile<<"#-----growth"<<endl;
    init_int phsGr           // growth
    !!CHECK1(phsGr);
    //--pr(molt-to-maturity)
    !!CheckFile<<"#-----pr(molt-to-maturity)"<<endl;
    init_int phsPrM2M        //molt to maturity
    !!CHECK1(phsPrM2M);
    //--surveys
    !!CheckFile<<"#-----surveys"<<endl;
    ////catchabilities
    !!CheckFile<<"#-------catchabilities"<<endl;
    init_int phsSrvQM1  //males
    init_int phsSrvQM2  //males
    init_int phsSrvQF1  //females
    init_int phsSrvQF2  //females
    !!CHECK1(phsSrvQM1);    
    !!CHECK1(phsSrvQM2);
    !!CHECK1(phsSrvQF1);    
    !!CHECK1(phsSrvQF2);
    ////selectivities
    !!CheckFile<<"#-------selectivities"<<endl;
    init_int phsSelSrvM1  //males
    init_int phsSelSrvM2  //males
    init_int phsSelSrvF1  //females
    init_int phsSelSrvF2  //females
    !!CHECK1(phsSelSrvM1);    
    !!CHECK1(phsSelSrvM2);
    !!CHECK1(phsSelSrvF1);    
    !!CHECK1(phsSelSrvF2);
    //--fisheries
    !!CheckFile<<"#-----fisheries"<<endl;
    ////fishing mortality or capture rates
    !!CheckFile<<"#-------fishing mortality or capture rates"<<endl;
    //--males
    init_int phsTCFM // directed Tanner crab fishery
    init_int phsSCFM // snow crab fishery bycatch
    init_int phsRKFM // BBRKC fishery bycatch
    init_int phsGTFM // groundfish fisheries bycatch
    !!CHECK1(phsTCFM);    
    !!CHECK1(phsSCFM);
    !!CHECK1(phsRKFM);    
    !!CHECK1(phsGTFM);
    //--females
    init_int phsTCFF ///< directed fishery bycatch
    init_int phsSCFF ///< snow crab fishery bycatch
    init_int phsRKFF ///< BBRKC fishery bycatch
    init_int phsGTFF ///< groundfish fisheries bycatch
    !!CHECK1(phsTCFF);    
    !!CHECK1(phsSCFF);
    !!CHECK1(phsRKFF);    
    !!CHECK1(phsGTFF);
    ////selectivity/retention functions
    !!CheckFile<<"#-------retention function/selectivities"<<endl;
    init_int phsRet_TCFM //retention
    init_int phsSel_TCFM //directed fishery male total catch
    init_int phsSel_TCFF //directed fishery female bycatch
    init_int phsSel_SCFM //snow crab fishery male bycatch
    init_int phsSel_SCFF //snow crab fishery female bycatch
    init_int phsSel_RKFM //BBRKC fishery male bycatch
    init_int phsSel_RKFF //BBRKC fishery female bycatch
    init_int phsSel_GTFM //groundfish fisheries male bycatch
    init_int phsSel_GTFF //groundfish fisheries female bycatch
    !!CHECK1(phsRet_TCFM);
    !!CHECK1(phsSel_TCFM);    
    !!CHECK1(phsSel_TCFF);
    !!CHECK1(phsSel_SCFM);    
    !!CHECK1(phsSel_SCFF);
    !!CHECK1(phsSel_RKFM);    
    !!CHECK1(phsSel_RKFF);
    !!CHECK1(phsSel_GTFM);    
    !!CHECK1(phsSel_GTFF);
    ////effort extrapolation
    init_int phsLnEffXtr_TCF  ///< TCF effort extrapolation parameter
    init_int phsLnEffXtr_SCF  ///< SCF effort extrapolation parameter
    init_int phsLnEffXtr_RKF  ///< RKF effort extrapolation parameter
    init_int phsLnEffXtr_GTF  ///< GTF effort extrapolation parameter
    !!CHECK1(phsLnEffXtr_TCF);    
    !!CHECK1(phsLnEffXtr_SCF);    
    !!CHECK1(phsLnEffXtr_RKF);    
    !!CHECK1(phsLnEffXtr_GTF);    
 LOCAL_CALCS
    if (useSomertonOtto1) {
        //turn off associated logistic function parameters
        phsSrvQM1   = -1;
        phsSelSrvM1 = -1;
        CheckFile<<"#--Using Somerton & Otto for survey period 1. TURNING male parameter estimation OFF--"<<endl;
    }
    if ((useSomertonOtto2)&&(useSomertonOtto3)) {
        //turn off associated logistic function parameters
        phsSrvQM2   = -1;
        phsSelSrvM2 = -1;
        CheckFile<<"#--Using Somerton & Otto for survey periods 2 & 3. TURNING male parameter estimation OFF--"<<endl;
    }
 END_CALCS
        
    //initial parameter values
    !!CheckFile<<"#---parameter initial values"<<endl;
    ////recruitment
    !!CheckFile<<"#-----recruitment"<<endl;
    init_number inpMnLnRecInit   // Mean ln-scale total "historic" recruitment
    init_number inpMnLnRec       // Mean ln-scale total "current" recruitment
    init_number inpRecAlpha      // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    init_number inpRecBeta       // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    !!CHECK1(inpMnLnRecInit);    
    !!CHECK1(inpMnLnRec);    
    !!CHECK1(inpRecAlpha);    
    !!CHECK1(inpRecBeta);    
    ////natural mortality
    !!CheckFile<<"#-----natural mortality multipliers"<<endl;
    init_number inpMfac_Imm      // natural mortality multiplier for immature females and males
    init_number inpMfac_MatM     // natural mortality multiplier for mature males
    init_number inpMfac_MatF     // natural mortality multiplier for mature females
    init_vector inpMfac_Big(1,nSXs)  // natural mortality multipliers for mature crab, 1980-1984
    !!CHECK1(inpMfac_Imm);    
    !!CHECK1(inpMfac_MatM);    
    !!CHECK1(inpMfac_Big);    
    ////growth
    !!CheckFile<<"#-----growth"<<endl;
    init_number inpGrAF1       // Female growth-increment
    init_number inpGrBF1       // Female growth-increment
    init_number inpGrAM1       // Male growth-increment
    init_number inpGrBM1       // Male growth-increment
    init_vector inpGrBeta_x(1,nSXs)  // Growth beta
    !!CHECK1(inpGrAF1);    
    !!CHECK1(inpGrBF1);    
    !!CHECK1(inpGrAM1);    
    !!CHECK1(inpGrBM1);    
    !!CHECK1(inpGrBeta_x);    
    ////pr(molt-to-maturity)
    !!CheckFile<<"#-----pr(molt-to-maturity)"<<endl;
    !!CheckFile<<"#<none>"<<endl;
    ////surveys
    !!CheckFile<<"#-----surveys"<<endl;
    //--catchabilities
    !!CheckFile<<"#-------catchabilities"<<endl;
    ////--males
    init_number inpSrv1_QM
    init_number inpSrv2_QM
    !!CHECK1(inpSrv1_QM);    
    !!CHECK1(inpSrv2_QM);    
    ////--females
    init_number inpSrv1_QF
    init_number inpSrv2_QF
    !!CHECK1(inpSrv1_QF);    
    !!CHECK1(inpSrv2_QF);    
    ////--selectivities
    !!CheckFile<<"#-------selectivities"<<endl;
    //////--males
    init_number inpSrv1M_z50
    init_number inpSrv1M_dz5095
    !!CHECK1(inpSrv1M_z50);    
    !!CHECK1(inpSrv1M_dz5095);    
    //////--1982+
    init_number inpSrv2M_z50   
    init_number inpSrv2M_dz5095
    !!CHECK1(inpSrv2M_z50);    
    !!CHECK1(inpSrv2M_dz5095);    
    //--females
    ////--1974-1981
    init_number inpSrv1F_z50
    init_number inpSrv1F_dz5095    
    !!CHECK1(inpSrv1F_z50);    
    !!CHECK1(inpSrv1F_dz5095);    
    ////--1982+
    init_number inpSrv2F_z50
    init_number inpSrv2F_dz5095
    !!CHECK1(inpSrv2F_z50);    
    !!CHECK1(inpSrv2F_dz5095);    
    //fisheries
    !!CheckFile<<"#-----fisheries"<<endl;
    //--log-scale fishing mortality or capture rates
    !!CheckFile<<"#-------fishing mortality/capture rates"<<endl;
    init_number inpAvgLnF_TCF // directed Tanner crab fishery
    init_number inpAvgLnF_SCF // snow crab fishery bycatch
    init_number inpAvgLnF_RKF // BBRKC fishery bycatch
    init_number inpAvgLnF_GTF // groundfish fisheries bycatch
    !!CHECK1(inpAvgLnF_TCF);    
    !!CHECK1(inpAvgLnF_SCF);
    !!CHECK1(inpAvgLnF_RKF);    
    !!CHECK1(inpAvgLnF_GTF);
    ////--log-scale offsets for females
    init_number inpAvgLnF_TCFF  ///< female offset to ln-scale mean fishing mortality in directed fishery
    init_number inpAvgLnF_SCFF  ///< female offset to ln-scale mean fishing mortality in snow crab fishery
    init_number inpAvgLnF_RKFF  ///< female offset to ln-scale mean fishing mortality in BBRKC fishery
    init_number inpAvgLnF_GTFF  ///< female offset to ln-scale mean fishing mortality in groundfish trawl fisheries
    !!CHECK1(inpAvgLnF_TCFF);    
    !!CHECK1(inpAvgLnF_SCFF);
    !!CHECK1(inpAvgLnF_RKFF);    
    !!CHECK1(inpAvgLnF_GTFF);
    //directed fishery selectivity and retention
    // Retention functions
    !!CheckFile<<"#-------retention functions"<<endl;
    //-- styr-1990
    init_number inpRetTCFM_slpA1
    init_number inpRetTCFM_z50A1
    //-- 1991+  
    init_number inpRetTCFM_slpA2
    init_number inpRetTCFM_z50A2
    !!CHECK1(inpRetTCFM_slpA1);    
    !!CHECK1(inpRetTCFM_z50A1);
    !!CHECK1(inpRetTCFM_slpA2);    
    !!CHECK1(inpRetTCFM_z50A2);
    // Selectivity functions
    !!CheckFile<<"#-------selectivities"<<endl;
    !!CheckFile<<"###--TCF"<<endl;
    //--males, styr-1996
    init_number inpSelTCFM_slpA1      
    //--males, 2005+
    init_number inpSelTCFM_slpA2      
    init_number inpSelTCFM_mnLnZ50A2
    !!CHECK1(inpSelTCFM_slpA1);    
    !!CHECK1(inpSelTCFM_slpA2);    
    !!CHECK1(inpSelTCFM_mnLnZ50A2);
    //--females, styr+
    init_number inpSelTCFF_slp
    init_number inpSelTCFF_z50
    !!CHECK1(inpSelTCFF_slp);    
    !!CHECK1(inpSelTCFF_z50);
    // snow crab fishery bycatch selectivity
    !!CheckFile<<"###--SCF"<<endl;
    //--males styr-1996
    init_number inpSelSCFM_slpA1 
    init_number inpSelSCFM_z50A1
    init_number inpSelSCFM_slpD1
    init_number inpSelSCFM_lnZ50D1
    !!CHECK1(inpSelSCFM_slpA1);    
    !!CHECK1(inpSelSCFM_z50A1);
    !!CHECK1(inpSelSCFM_slpD1);    
    !!CHECK1(inpSelSCFM_lnZ50D1);
    //--males, 1997-2004
    init_number inpSelSCFM_slpA2  
    init_number inpSelSCFM_z50A2
    init_number inpSelSCFM_slpD2
    init_number inpSelSCFM_lnZ50D2
    !!CHECK1(inpSelSCFM_slpA2);    
    !!CHECK1(inpSelSCFM_z50A2);
    !!CHECK1(inpSelSCFM_slpD2);    
    !!CHECK1(inpSelSCFM_lnZ50D2);
    //--males, 2005+
    init_number inpSelSCFM_slpA3
    init_number inpSelSCFM_z50A3
    init_number inpSelSCFM_slpD3
    init_number inpSelSCFM_lnZ50D3  
    !!CHECK1(inpSelSCFM_slpA3);    
    !!CHECK1(inpSelSCFM_z50A3);
    !!CHECK1(inpSelSCFM_slpD3);    
    !!CHECK1(inpSelSCFM_lnZ50D3);
    //females, styr-1996
    init_number inpSelSCFF_slpA1 
    init_number inpSelSCFF_z50A1
    //--females, 1997-2004
    init_number inpSelSCFF_slpA2 
    init_number inpSelSCFF_z50A2
    //--females, 2005+
    init_number inpSelSCFF_slpA3 
    init_number inpSelSCFF_z50A3
    !!CHECK1(inpSelSCFF_slpA1);    
    !!CHECK1(inpSelSCFF_z50A1);
    !!CHECK1(inpSelSCFF_slpA2);    
    !!CHECK1(inpSelSCFF_z50A2);
    !!CHECK1(inpSelSCFF_slpA3);    
    !!CHECK1(inpSelSCFF_z50A3);
    // BBRKC fishery bycatch selectivity
    !!CheckFile<<"###--RKF"<<endl;
    //--males, styr-1996
    init_number inpSelRKFM_slpA1  
    init_number inpSelRKFM_z50A1
    //--males, 1997-2004
    init_number inpSelRKFM_slpA2  
    init_number inpSelRKFM_z50A2
    //--males, 2005+
    init_number inpSelRKFM_slpA3  
    init_number inpSelRKFM_z50A3
    //--females, styr-1996
    init_number inpSelRKFF_slpA1 
    init_number inpSelRKFF_z50A1 
    //--females, 1997-2004
    init_number inpSelRKFF_slpA2 
    init_number inpSelRKFF_z50A2 
    //--females, 2005+
    init_number inpSelRKFF_slpA3
    init_number inpSelRKFF_z50A3 
    !!CHECK1(inpSelRKFM_slpA1);    
    !!CHECK1(inpSelRKFM_z50A1);
    !!CHECK1(inpSelRKFM_slpA2);    
    !!CHECK1(inpSelRKFM_z50A2);
    !!CHECK1(inpSelRKFM_slpA3);    
    !!CHECK1(inpSelRKFM_z50A3);
    !!CHECK1(inpSelRKFF_slpA1);    
    !!CHECK1(inpSelRKFF_z50A1);
    !!CHECK1(inpSelRKFF_slpA2);    
    !!CHECK1(inpSelRKFF_z50A2);
    !!CHECK1(inpSelRKFF_slpA3);    
    !!CHECK1(inpSelRKFF_z50A3);
    // groundfish fisheries bycatch selectivity 
    !!CheckFile<<"###--GTF"<<endl;
    //--males, 1973-1987
    init_number inpSelGTFM_slpA1
    init_number inpSelGTFM_z50A1
    //--males, 1988-1996
    init_number inpSelGTFM_slpA2
    init_number inpSelGTFM_z50A2
    //--males, 1997+
    init_number inpSelGTFM_slpA3
    init_number inpSelGTFM_z50A3
    //--females, 1973-1987
    init_number inpSelGTFF_slpA1
    init_number inpSelGTFF_z50A1
    //--females, 1988-1996
    init_number inpSelGTFF_slpA2
    init_number inpSelGTFF_z50A2 
    //females, 1997+
    init_number inpSelGTFF_slpA3
    init_number inpSelGTFF_z50A3
    !!CHECK1(inpSelGTFM_slpA1);    
    !!CHECK1(inpSelGTFM_z50A1);
    !!CHECK1(inpSelGTFM_slpA2);    
    !!CHECK1(inpSelGTFM_z50A2);
    !!CHECK1(inpSelGTFM_slpA3);    
    !!CHECK1(inpSelGTFM_z50A3);
    !!CHECK1(inpSelGTFF_slpA1);    
    !!CHECK1(inpSelGTFF_z50A1);
    !!CHECK1(inpSelGTFF_slpA2);    
    !!CHECK1(inpSelGTFF_z50A2);
    !!CHECK1(inpSelGTFF_slpA3);    
    !!CHECK1(inpSelGTFF_z50A3);

    //effort extrapolation parameters
    init_number inpLnEffXtr_TCF  ///< TCF effort extrapolation parameter
    init_number inpLnEffXtr_SCF  ///< SCF effort extrapolation parameter
    init_number inpLnEffXtr_RKF  ///< RKF effort extrapolation parameter
    init_number inpLnEffXtr_GTF  ///< GTF effort extrapolation parameter
    !!CHECK1(inpLnEffXtr_TCF);
    !!CHECK1(inpLnEffXtr_SCF);
    !!CHECK1(inpLnEffXtr_RKF);
    !!CHECK1(inpLnEffXtr_GTF);
         
    !!CheckFile<<"#--end of file check!"<<endl;
    init_int chkEOF
    !!CHECK1(chkEOF);
 LOCAL_CALCS
    if (chkEOF!=999){
        cout<<"Problem reading control file!!"<<endl;
        cout<<"chkEOF = "<<chkEOF<<endl;
        cout<<"Should be '999'"<<endl;
        cout<<"Exiting..."<<endl;
        exit(-1);
    }
 END_CALCS
    //Finished reading control file
    !!CheckFile<<"--Finished reading control file."<<endl;
    !!CheckFile<<"--------------------------------------"<<endl;
    //TEMPORARY VARIABLES
    int phase_logistic_sel
    int survsel1_phase
    int survsel_phase
    int maturity_switch
 LOCAL_CALCS
    phase_logistic_sel = phsRet_TCFM;
    survsel1_phase = phsSelSrvM1;
    survsel_phase = phsSelSrvM2;
    maturity_switch = optInitM2M;
    CHECK1(phase_logistic_sel);
    CHECK1(survsel1_phase);
    CHECK1(survsel_phase);
    CHECK1(maturity_switch);
 END_CALCS
    
    // the rest are working variables 
    
    matrix obsDscBioMortTCF_xn(1,nSXs,1,nObsDscTCF)   // observed discard MORTALITY in TCF (1000's tons)
    matrix obsDscBioMortSCF_xn(1,nSXs,1,nObsDscSCF)   // observed discard MORTALITY in SCF (1000's tons)
    matrix obsDscBioMortRKF_xn(1,nSXs,1,nObsDscRKF)   // observed discard MORTALITY in RKF (1000's tons)
    vector obsDscBioMortGTF_n(1,nObsDscGTF)    // observed discard MORTALITY in GTF (1000's tons)
    vector obsTotBioMortTCFM_n(1,nObsDscTCF)   // observed total male MORTALITY in TCF (1000's tons)
 LOCAL_CALCS
      //nums are now (20160608) kept in millions, so n*wt(kg) gives biomass in 1000's t
//    obsSrvNum_n=obsSrvNum_n*1000000;       // survey numbers read in are millions of crab (but /1000 above so obsSrvNum_n in thousands of crab, now)
//    wt_xmz(FEMALE) = wt_xmz(FEMALE)*0.001; // change weights from kg to tons
//    wt_xmz(  MALE) = wt_xmz(  MALE)*0.001; // change weights from kg  to tons
    //apply discard mortality to observed discards
    obsDscBioMortTCF_xn = hm_pot*obsDscBioTCF_xn;//apply discard mortality to catches
    obsDscBioMortSCF_xn = hm_pot*obsDscBioSCF_xn;
    obsDscBioMortRKF_xn = hm_pot*obsDscBioRKF_xn;
    obsDscBioMortGTF_n = hm_trawl*obsDscBioGTF_n;
    // Calculate TOTAL male mortality in Tanner crab fishery
    obsTotBioMortTCFM_n.initialize();
    obsTotBioMortTCFM_n = obsDscBioMortTCF_xn(MALE)+obsRetCatchBio_y(yrsObsDscTCF_n);
//    for (int i=1;i<=nObsDscTCF;i++) {
//        obsTotBioMortTCFM_n(i) = obsDscBioMortTCF_xn(MALE,i)+obsRetCatchBio_y(yrsObsDscTCF_n(i));
//    }
    CheckFile<<"yrsObsDscTCF_n"<<endl<<yrsObsDscTCF_n<<endl;
    CheckFile<<"obsDscBioMortTCF_xn(  MALE)"<<endl<<obsDscBioMortTCF_xn(  MALE)<<endl;
    CheckFile<<"obsDscBioMortTCF_xn(FEMALE)"<<endl<<obsDscBioMortTCF_xn(FEMALE)<<endl;
    CheckFile<<"obsRetBioMortTCFM"<<endl<<obsRetCatchBio_y(yrsObsDscTCF_n)<<endl;
    CheckFile<<"obsTotBioMortTCFM_n"<<endl<<obsTotBioMortTCFM_n<<endl;
    CheckFile<<"yrsObsDscSCF"<<endl<<yrsObsDscSCF<<endl;
    CheckFile<<"obsDscBioMortSCF_xn(  MALE)"<<endl<<obsDscBioMortSCF_xn(  MALE)<<endl;
    CheckFile<<"obsDscBioMortSCF_xn(FEMALE)"<<endl<<obsDscBioMortSCF_xn(FEMALE)<<endl;
    CheckFile<<"yrsObsDscRKF"<<endl<<yrsObsDscRKF<<endl;
    CheckFile<<"obsDscBioMortRKF_xn(  MALE)"<<endl<<obsDscBioMortRKF_xn(  MALE)<<endl;
    CheckFile<<"obsDscBioMortRKF_xn(FEMALE)"<<endl<<obsDscBioMortRKF_xn(FEMALE)<<endl;
    CheckFile<<"yrsObsDscGTF"<<endl<<yrsObsDscGTF<<endl;
    CheckFile<<"obsDscBioMortGTF_n"<<endl<<obsDscBioMortGTF_n<<endl;
    
 END_CALCS
    
    3darray obsPrNatZ_TCFR_sn(1,nSCs,1,nObsRetZCsTCF,1,nZBs)    // length-frequency of retained catch?  
    3darray obsPrNatZ_TCFM_snz(1,nSCs,1,nObsZCsTCFM,1,nZBs)     // length-frequency of total male catch in directed fishery
    matrix obsPrNatZ_TCFF_nz(1,nObsZCsTCFF,1,nZBs)              // female discards
    3darray obsPrNatZ_GTF_xnz(1,nSXs,1,nObsZCsGTF,1,nZBs)       // bycatch in trawl fishery
    
    4darray obsSrvPrNatZ_mxnz(1,nMSs,1,nSXs,1,nObsZCsSrv,1,nZBs)         // normalized survey size comps by maturity state, sex
    5darray obsSrvPrNatZ_msxnz(1,nMSs,1,nSCs,1,nSXs,1,nObsZCsSrv,1,nZBs) // normalized survey size comps by maturity state, shell condition, sex
    matrix obsSrvImmBio_xy(1,nSXs,styr,endyr)                            // Survey immature biomass (by sex)
    matrix obsSrvMatBio_xy(1,nSXs,styr,endyr)                            // Survey mature biomass (by sex)
    3darray obsSrvImmNum_sxy(1,nSCs,1,nSXs,styr,endyr)                   // Survey immature numbers (by shell condition, sex)
    3darray obsSrvMatNum_sxy(1,nSCs,1,nSXs,styr,endyr)                   // Survey mature numbers (by shell condition, sex)
    
    3darray obsPrNatZ_SCF_xnz(1,nSXs,1,nObsZCsSCF,1,nZBs)  
    3darray obsPrNatZ_RKF_xnz(1,nSXs,1,nObsZCsRKF,1,nZBs)    
    
    vector obsSrvNumLegal_n(1,nObsZCsSrv)            // Legal male numbers, in survey
    vector obsSrvBioLegal_n(1,nObsZCsSrv)            // Legal male biomass, in survey
    
    3darray obsSrvNum_xyz(1,nSXs,styr,endyr,1,nZBs)  // Survey numbers, by sex
    vector obsSrvNum_y(styr,endyr)                   // Total survey numbers, all years
    matrix obsSrvBio_xy(1,nSXs,styr,endyr)           // Survey biomass, by sex
    vector obsSrvBio_y(styr,endyr)                   // Total survey biomass
    
    number mnEff_SCF; ///< mean effort in SCF over period with observed bycatch data
    number mnEff_RKF; ///< mean effort in RKF over period with observed bycatch data
    
    vector selSO_z(1,nZBs); //Somerton & Otto survey selectivity function
    !!selSO_z = sel_som(1)/(1.+sel_som(2)*mfexp(-1.*sel_som(3)*zBs));
    !!CheckFile<<"selSO_z"<<endl<<selSO_z<<endl;
    
 LOCAL_CALCS
    dmX   = "x=c("+qt+STR_FEMALE+qt    +cc+ qt+STR_MALE+qt+")";
    dmS   = "s=c("+qt+STR_NEW_SHELL+qt +cc+ qt+STR_OLD_SHELL+qt+")";
    dmM   = "m=c("+qt+STR_IMMATURE+qt  +cc+ qt+STR_MATURE+qt+")";
    dmY   = "y="+str(styr)+":"+str(endyr);
    dmYm1 = "y="+str(styr)+":"+str(endyr-1);  
    dmZ   = "z=c("+wts::to_qcsv(zBs)+")";
 END_CALCS
            
    !!CheckFile<<"End of DATA_SECTION-------------------------"<<endl;
    !!CheckFile<<"--------------------------------------------"<<endl<<endl;

// =======================================================================
// =======================================================================
INITIALIZATION_SECTION
    pPrM2MF -1.0
    pPrM2MM -1.0
 
// =======================================================================
// =======================================================================
PARAMETER_SECTION
    
    //growth
    init_bounded_number pGrAF1(0.4,0.7,phsGr)      // Female growth-increment
    init_bounded_number pGrBF1(0.6,1.2,phsGr)      // Female growth-increment
    init_bounded_number pGrAM1(0.3,0.6,phsGr)      // Male growth-increment
    init_bounded_number pGrBM1(0.7,1.2,phsGr)      // Male growth-increment
    init_bounded_vector pGrBeta_x(1,nSXs,0.5,1.0,-1)  // Growth beta:  NOT estimated

    //natural mortality
    init_bounded_number pMfac_Imm(0.2,2.0,phsM)                // natural mortality multiplier for immature females and males
    init_bounded_number pMfac_MatM(0.1,1.9,phsM)               // natural mortality multiplier for mature males
    init_bounded_number pMfac_MatF(0.1,1.9,phsM)               // natural mortality multiplier for mature females
    init_bounded_vector pMfac_Big(1,nSXs,0.1,10.0,phsBigM)     // mult. on 1980-1984 M for mature males and females  

    //recruitment
    init_bounded_number pRecAlpha(11.0,12.0,-1)              // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    init_bounded_number pRecBeta(3.0,5.0,-1)                 // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    init_number pMnLnRec(phsMnLnRec)                                                     // Mean ln-scale total "current" recruitment
    init_bounded_dev_vector pRecDevs(mnYrRecCurr,endyr,-15,15,phsRecDevs)                // "current" ln-scale total recruitment devs
    init_number pMnLnRecInit(phsMnLnRec)                                                 // Mean ln-scale total "historic" recruitment
    init_bounded_dev_vector pRecDevsHist(mnYrRecDevsHist,mnYrRecCurr-1,-15,15,phsRecDevs)// "historic" ln-scale total recruitment devs
    
    //20150601: changed ...Fm... to ...F_... because of ambiguity as to whether
    //they represent fishing mortality (original model) of capture (gmacs) rates
    init_number pAvgLnF_TCF(phsTCFM)                                 //log-scale mean directed fishing mortality
    init_bounded_dev_vector pF_DevsTCF(1,nYrsTCF,-15,15,phsTCFM+1)   //log-scale directed fishing mortality devs IMPORTANT CHANGE: USED TO BE "1966,endyr-12"
    init_number pAvgLnF_GTF(phsGTFM)                                 // fishing mortality (trawl)
    init_bounded_dev_vector pF_DevsGTF(1973,endyr-1,-15,15,phsGTFM+1)// trawl fishery f-devs       (IMPORTANT CHANGE: used to be "endyr") 1973 seems OK
    init_number pAvgLnF_SCF(phsSCFM)                                 // fishing mortality snow crab fishery discards
    init_bounded_dev_vector pF_DevsSCF(1992,endyr-1,-15,15,phsSCFM+1)// snow crab fishery f-devs   (IMPORTANT CHANGE: used to be "endyr")  1992 is OK
    init_bounded_number pAvgLnF_RKF(-10,5,phsRKFM)                   // fishing mortality red king crab fishery discards
    init_bounded_dev_vector pF_DevsRKF(1,nObsDscRKF,-15,15,phsRKFM+1)// 
    
    // Retention function
    //-- styr-1990
    init_bounded_number pRetTCFM_slpA1(00.250,001.001,phsRet_TCFM)
    init_bounded_number pRetTCFM_z50A1(85.000,160.000,phsRet_TCFM)
    //-- 1991+  
    init_bounded_number pRetTCFM_slpA2(00.250,002.001,phsRet_TCFM)
    init_bounded_number pRetTCFM_z50A2(85.000,160.000,phsRet_TCFM)
    
    //Directed fishery male selectivity pattern
    //-- styr-1996
    init_bounded_number pSelTCFM_slpA1(00.05,000.75,phsSel_TCFM)      
    //-- 2005+
    init_bounded_number pSelTCFM_slpA2(0.1,0.4,phsSel_TCFM)      
    init_bounded_number pSelTCFM_mnLnZ50A2(4.0,5.0,phsSel_TCFM)
    //-- 1991+
    init_bounded_dev_vector pSelTCFM_devsZ50(1,nSelTCFM_devsZ50,-bndSelTCFM_devsZ50,bndSelTCFM_devsZ50,phsSel_TCFM)
    
    // Directed fishery female bycatch selectivity pattern
    //-- styr+
    init_bounded_number pSelTCFF_slp(00.1,000.4,phsSel_TCFF)
    init_bounded_number pSelTCFF_z50(80.0,150.0,phsSel_TCFF)
    
    // Snow crab fishery female bycatch selectivity pattern
    //-- 1989-1996
    init_bounded_number pSelSCFF_slpA1(00.05,000.5,phsSel_SCFF)
    init_bounded_number pSelSCFF_z50A1(50.00,150.0,phsSel_SCFF)
    //-- 1997-2004
    init_bounded_number pSelSCFF_slpA2(00.05,000.5,phsSel_SCFF)
    init_bounded_number pSelSCFF_z50A2(50.00,120.0,phsSel_SCFF)
    //-- 2005+
    init_bounded_number pSelSCFF_slpA3(00.05,000.5,phsSel_SCFF)
    init_bounded_number pSelSCFF_z50A3(50.00,120.0,phsSel_SCFF)
    
    // Snow crab fishery male bycatch selectivity pattern
    //-- 1989-1996
    init_bounded_number pSelSCFM_slpA1(00.1,000.5,phsSel_SCFM)  
    init_bounded_number pSelSCFM_z50A1(40.0,140.0,phsSel_SCFM)
    init_bounded_number pSelSCFM_slpD1(00.1,000.5,phsSel_SCFM)
    init_bounded_number pSelSCFM_lnZ50D1(2,4.5,phsSel_SCFM)
    //-- 1997-2004
    init_bounded_number pSelSCFM_slpA2(00.1,000.5,phsSel_SCFM)
    init_bounded_number pSelSCFM_z50A2(40.0,140.0,phsSel_SCFM)
    init_bounded_number pSelSCFM_slpD2(00.1,000.5,phsSel_SCFM)
    init_bounded_number pSelSCFM_lnZ50D2(2,4.5,phsSel_SCFM)
    //-- 2005+
    init_bounded_number pSelSCFM_slpA3(00.1,000.5,phsSel_SCFM)
    init_bounded_number pSelSCFM_z50A3(40.0,140.0,phsSel_SCFM)
    init_bounded_number pSelSCFM_slpD3(00.1,000.5,phsSel_SCFM)
    init_bounded_number pSelSCFM_lnZ50D3(2,4.5,phsSel_SCFM)  
    
    // Red king crab fishery female bycatch selectivity pattern
    //-- 
    init_bounded_number pSelRKFF_slpA1(00.05,000.5,phsSel_RKFF) 
    init_bounded_number pSelRKFF_z50A1(50.00,150.0,phsSel_RKFF) 
    //-- 
    init_bounded_number pSelRKFF_slpA2(00.05,000.5,phsSel_RKFF) 
    init_bounded_number pSelRKFF_z50A2(50.00,150.0,phsSel_RKFF) 
    //-- 
    init_bounded_number pSelRKFF_slpA3(00.05,000.5,phsSel_RKFF) 
    init_bounded_number pSelRKFF_z50A3(50.00,170.0,phsSel_RKFF) 
    
    // Red king crab fishery male bycatch selectivity pattern
    //-- 
    init_bounded_number pSelRKFM_slpA1(.01,.50,phsSel_RKFM)          
    init_bounded_number pSelRKFM_z50A1(95.0,150.0,phsSel_RKFM)
    //-- 
    init_bounded_number pSelRKFM_slpA2(.01,.50,phsSel_RKFM)          
    init_bounded_number pSelRKFM_z50A2(95.0,150.0,phsSel_RKFM)
    //-- 
    init_bounded_number pSelRKFM_slpA3(.01,.50,phsSel_RKFM)          
    init_bounded_number pSelRKFM_z50A3(95.0,150.0,phsSel_RKFM)
    
    // Groundfish fisheries female bycatch selectivity pattern
    //-- 1973-1987
    init_bounded_number pSelGTFF_slpA1(0.01,0.5,phsSel_GTFF)
    init_bounded_number pSelGTFF_z50A1(40.0,125.01,phsSel_GTFF)
    //-- 1988-1996
    init_bounded_number pSelGTFF_slpA2(0.005,0.5,phsSel_GTFF)
    init_bounded_number pSelGTFF_z50A2(40.0,250.01,phsSel_GTFF) 
    //-- 1997+
    init_bounded_number pSelGTFF_slpA3(0.01,0.5,phsSel_GTFF)
    init_bounded_number pSelGTFF_z50A3(40.0,150.01,phsSel_GTFF)

    // Groundfish fisheries male bycatch selectivity pattern
    ///-- styr-1987
    init_bounded_number pSelGTFM_slpA1(0.01,0.5,phsSel_GTFM)
    init_bounded_number pSelGTFM_z50A1(40.0,120.01,phsSel_GTFM)
    //-- 1988-1996
    init_bounded_number pSelGTFM_slpA2(0.01,0.5,phsSel_GTFM)
    init_bounded_number pSelGTFM_z50A2(40.0,120.01,phsSel_GTFM)
    //-- 1997+
    init_bounded_number pSelGTFM_slpA3(0.01,0.5,phsSel_GTFM)
    init_bounded_number pSelGTFM_z50A3(40.0,120.01,phsSel_GTFM)

    //Survey-related parameters for males
    //-- males: pre-1982 
    init_bounded_number pSrv1_QM(0.50,1.001,phsSrvQM1)
    init_bounded_number pSrv1M_dz5095(0.0,100.0,phsSelSrvM1)
    init_bounded_number pSrv1M_z50(0.0,90.0,phsSelSrvM1)
    //--males: 1982+
    init_bounded_number pSrv2_QM(0.2,2.0,phsSrvQM2)
    init_bounded_number pSrv2M_dz5095(0.0,100.0,phsSelSrvM2)
    init_bounded_number pSrv2M_z50(0.0,69.0,phsSelSrvM2)    

    //Probability of molt-to-maturity parameters 
    init_bounded_vector pPrM2MF(1,16,lbPrM2M,ubPrM2M,phsPrM2M)  //females
    init_bounded_vector pPrM2MM(1,nZBs,lbPrM2M,ubPrM2M,phsPrM2M)//males
    
    //Survey-related parameters for females
    //-- females: pre-1982 
    init_bounded_number pSrv1_QF(0.5,1.001,phsSrvQF1)
    init_bounded_number pSrv1F_dz5095(0.0,100.0,phsSelSrvF1)    
    init_bounded_number pSrv1F_z50(-200.0,100.01,phsSelSrvF1)
    //--females: 1982+
    init_bounded_number pSrv2_QF(0.2,1.0,phsSrvQF2)
    init_bounded_number pSrv2F_dz5095(0.0,100.0,phsSelSrvF2)
    init_bounded_number pSrv2F_z50(-50.0,69.0,phsSelSrvF2)
    
    //Offsets to fully-selected male fishing mortality/fishery capture rate for females
    init_bounded_number pAvgLnF_TCFF(-5.0,5.0,phsTCFF)  ///< female offset to ln-scale mean fishing mortality in directed fishery
    init_bounded_number pAvgLnF_SCFF(-5.0,5.0,phsSCFF)  ///< female offset to ln-scale mean fishing mortality in snow crab fishery
    init_bounded_number pAvgLnF_RKFF(-5.0,5.0,phsRKFF)  ///< female offset to ln-scale mean fishing mortality in BBRKC fishery
    init_bounded_number pAvgLnF_GTFF(-5.0,5.0,phsGTFF)  ///< female offset to ln-scale mean fishing mortality in groundfish trawl fisheries

    //Effort extrapolation parameters
    init_bounded_number pLnEffXtr_TCF(-20.0,0.0,phsLnEffXtr_TCF)  ///< TCF effort extrapolation parameter
    init_bounded_number pLnEffXtr_SCF(-20.0,0.0,phsLnEffXtr_SCF)  ///< SCF effort extrapolation parameter
    init_bounded_number pLnEffXtr_RKF(-20.0,0.0,phsLnEffXtr_RKF)  ///< RKF effort extrapolation parameter
    init_bounded_number pLnEffXtr_GTF(-20.0,0.0,phsLnEffXtr_GTF)  ///< GTF effort extrapolation parameter
    ////end of estimated parameters///////////////
    
    3darray retFcn_syz(1,nSCs,styr,endyr-1,1,nZBs)    // Retention curve for males caught in directed fishery    (IMPORTANT CHANGE: used to be "endyr")
    3darray selTCFR_syz(1,nSCs,styr,endyr-1,1,nZBs)   // full selectivity for retained males in directed fishery (IMPORTANT CHANGE: used to be "endyr")
    3darray selTCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)   // selectivity for all males in directed fishery           (IMPORTANT CHANGE: used to be "endy1")
    vector  selTCFF_z(1,nZBs)                       // selectivity for females in directed fishery             
    3darray selSCF_cxz(1,3,1,nSXs,1,nZBs)      // 3D array to accommodate 3 selectivity periods 
    3darray selRKF_cxz(1,3,1,nSXs,1,nZBs)      // 3D array to accommodate 3 selectivity periods 
    3darray selGTF_cxz(1,3,1,nSXs,1,nZBs)      // 3D array to accommodate 3 selectivity periods
    
    matrix selSrv1_xz(1,nSXs,1,nZBs) // Survey selectivity 1 pre 1982
    matrix selSrv2_xz(1,nSXs,1,nZBs) // Survey selectivity 2 1982-1987
    matrix selSrv3_xz(1,nSXs,1,nZBs) // Survey selectivity 3 1988+
        
    3darray M_msx(1,nMSs,1,nSCs,1,nSXs);//natural mortality rates
    
    vector modPopNum_y(styr,endyr) // Total population numbers on July 1, endyr (output)
    vector modPopBio_y(styr,endyr) // Predicted biomass (determines sdrDepletion) 
    
    vector fspbio(styr,endyr)                       // Predicted female spawning biomass on July 1
    vector mspbio(styr,endyr)                       // Predicted   male spawning biomass on July 1  
    matrix modSrvImmBio_xy(1,nSXs,styr,endyr)          // Predicted immature biomass at survey time, as seen by survey
    matrix modSrvMatBio_xy(1,nSXs,styr,endyr)          // Predicted mature biomass at survey time, as seen by survey
    3darray modSrvImmNum_xsy(1,nSXs,1,nSCs,styr,endyr) // Predicted survey numbers for immature crab by sex, shell condition (output)
    3darray modSrvMatNum_xsy(1,nSXs,1,nSCs,styr,endyr) // Predicted survey numbers for mature crab by sex, shell condition (output)
    
    vector modPopNumLegal_y(styr,endyr)              // Number of legal males at time of survey 
    vector modPopBioLegal_y(styr,endyr)              // Biomass of legal males in survey (output)          
    
    matrix modSrvNum_xy(1,nSXs,styr,endyr)             // Survey abundance, by sex          
    matrix modSrvBio_xy(1,nSXs,styr,endyr)             // Survey biomass, by sex      
    vector modSrvNumLegal_y(styr,endyr)                // Survey-selected legal-sized males   
    matrix modSrvNumLegal_sy(1,nSCs,styr,endyr)        // Survey-selected legal-sized males, by shell condition  
    vector modSrvBioLegal_y(styr,endyr)                // Survey-selected legal male biomass  
    3darray modSrvNum_xyz(1,nSXs,styr,endyr,1,nZBs)    // Predicted survey abundance by length bin       
    4darray modSrvPrNatZ_mxyz(1,nMSs,1,nSXs,styr,endyr,1,nZBs)         // Predicted size comps (integrated over shell condition)
    5darray modSrvPrNatZ_msxyz(1,nMSs,1,nSCs,1,nSXs,styr,endyr,1,nZBs) // Predicted size comps 
    
    3darray modPrNatZ_TCFR_syz(1,nSCs,styr,endyr-1,1,nZBs)  // Predicted retained catch proportions
    3darray modPrNatZ_TCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)  // Predicted male   proportions in directed fishery (total catch)      
    matrix  modPrNatZ_TCFF_yz(styr,endyr-1,1,nZBs)          // Predicted female proportions in directed fishery
    3darray modPrNatZ_GTF_xyz(1,nSXs,styr,endyr-1,1,nZBs)   // Predicted trawl proportions                      
    3darray modPrNatZ_SCF_xyz(1,nSXs,styr,endyr-1,1,nZBs)   // Predicted snow crab fishery  proportions         
    3darray modPrNatZ_RKF_xyz(1,nSXs,styr,endyr-1,1,nZBs)   // Predicted red king crab proportions              
    
    3darray modNum_xyz(1,nSXs,styr,endyr,1,nZBs)                 // Total numbers by year, sex, size
    5darray modNum_yxmsz(styr,endyr,1,nSXs,1,nMSs,1,nSCs,1,nZBs) // Total numbers by year, sex, maturity state, shell condition, size
    4darray modNumAtAZ_xyaz(1,nSXs,styr,endyr,1,nages,1,nZBs)    // Age- and length-structure          
    3darray natlength_old(1,nSXs,styr,endyr,1,nZBs)           // Old-shell numbers by sex, length, and year
    3darray natlength_new(1,nSXs,styr,endyr,1,nZBs)           // New-shell numbers by sex, length, and year
    3darray natlength_imm(1,nSXs,styr,endyr,1,nZBs)             // Immature numbers by sex, length, and year
    3darray natlength_mat(1,nSXs,styr,endyr,1,nZBs)           // Mature numbers by sex, length, and year
    3darray natl_new_fishtime(1,nSXs,styr,endyr,1,nZBs)           // Numbers-at-length (new shell)
    3darray natl_old_fishtime(1,nSXs,styr,endyr,1,nZBs)           // Numbers-at-length (old shell)    
    5darray modFT_PopNum_yxmsz(styr,endyr,1,nSXs,1,nMSs,1,nSCs,1,nZBs)  // Numbers-at-length at fishing time          
    
    //--growth-related quantities
    3darray prGr_xzz(1,nSXs,1,nZBs,1,nZBs)    // length to length growth array
    matrix prMoltImm_xz(1,nSXs,1,nZBs)        // molting probabilities for female, male by length bin 
    matrix prMoltMat_xz(1,nSXs,1,nZBs)        // molting probs for mature female, male by length bin
    matrix meanPostMoltSize(1,nSXs,1,nZBs)    // Predicted mean post-moult sizes
    
    //--recruitment-related quantities
    vector rec_y(styr,endyr)      //arithmetic-scale total recruitment (millions)
    vector prRec_z(1,nZBs)        // Recruitment size frequency
    vector modPopXR_y(styr,endyr) // Population sex-ratio - output
    
    //changed endyr to endyr-1
    //20150601: changed fm... to f... because these could be fishing mortality OR capture rates
    //fully-selected fishery capture (gmacs option) OR mortality (original) rates
    //20160225: changing to matrices to incorporate sex-specific rates
    matrix fTCF_xy(1,nSXs,styr,endyr-1)    //directed fishery
    matrix fSCF_xy(1,nSXs,styr,endyr-1)    //snow crab fishery
    matrix fRKF_xy(1,nSXs,styr,endyr-1)    //BBRKC fishery
    matrix fGTF_xy(1,nSXs,styr,endyr-1)    //groundfish fisheries
    
    //ratios for extrapolating effort to fishing mortality
    number qSCF  //snow crab fishery
    number qRKF  //BBRKC fishery
                      
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
    3darray predRetNumMortTCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)   ///< Predicted retained catch mortality, by shell condition and size
    matrix predRetNumMortTCFM_yz(styr,endyr-1,1,nZBs)            ///< Predicted retained catch mortality, by size
    vector predRetBioMortTCFM_y(styr,endyr-1)                    ///< Retained catch mortality (biomass)
    3darray predTotNumMortTCFM_syz(1,nSCs,styr,endyr-1,1,nZBs)   ///< Predicted total catch mortality for males, by shell condition and size
    matrix predTotNumMortTCFM_yz(styr,endyr-1,1,nZBs)            ///< Predicted total catch mortality, by size
    vector predTotBioMortTCFM_y(styr,endyr-1)                    ///< Total catch mortality (biomass)
    3darray predDscNumMortTCF_xyz(1,nSXs,styr,endyr-1,1,nZBs)     ///< Predicted discard catch mortality, by sex and size
    matrix predDscBioMortTCF_xy(1,nSXs,styr,endyr-1)             ///< Predicted discard catch mortality, by sex
    
    //SCF
    3darray predDscNumMortSCF_xyz(1,nSXs,styr,endyr-1,1,nZBs)
    matrix predDscBioMortSCF_xy(1,nSXs,styr,endyr-1)
    
    //BBRKC
    3darray predDscNumMortRKF_xyz(1,nSXs,styr,endyr-1,1,nZBs)
    matrix predDscBioMortRKF_xy(1,nSXs,styr,endyr-1)
    
    //GTF
    3darray predDscNumMortGTF_xyz(1,nSXs,styr,endyr-1,1,nZBs)
    matrix predDscBioMortGTF_xy(1,nSXs,styr,endyr-1)                  ///< Trawl catch mortality (in biomass and by sex)
    
    //changed endyr to endyr-1   
    vector effnTCF_ret_y(styr,endyr-1)            // Effective sample sizes
    matrix effnTCF_tot_xy(1,nSXs,styr,endyr-1)    // Effective sample sizes
    matrix effnSCF_tot_xy(1,nSXs,styr,endyr-1)    // Effective sample sizes
    matrix effnRKF_tot_xy(1,nSXs,styr,endyr-1)    // Effective sample sizes
    vector effnGTF_tot_y(styr,endyr-1)            // Effective sample sizes
    
    vector effnSrv_y(styr,endyr) // Effective sample sizes
    
    // Offsets
    vector offset(1,NUM_LEN_LIKE)  //<-these are constants. should be in DATA_SECTION!!
                                              
    vector objfOut(1,NUM_FOUT);      // objective function components (weighted likelihoods and penalties)
    vector likeOut(1,NUM_FOUT);      //unweighted likelihood components and penalties
    vector wgtsOut(1,NUM_FOUT);      //weights
    
    // Penalties
    number nat_penalty
    number penal_rec         // Recruitment
    number fpen              // Penalties (misc)
    number af_penal          // Prior on af 
    number srv3q_penalty
    number am_penal          // Prior on am
    number bf_penal          // Prior on bf
    number bm_penal          // Prior on bm
    
    // Likelihood components
    vector zsRetMortBio_TCFR_y(styr,endyr-1)        ///< z-scores for fit to retained male catch mortality in directed TCF 
    vector zsTotMortBio_TCFM_n(1,nObsDscTCF)        ///< z-scores for fit to total male catch mortality in directed TCF 
    vector zsDscMortBio_TCFM_n(1,nObsDscTCF)        ///< z-scores for fit to discard male catch mortality in directed TCF 
    vector zsDscMortBio_TCFF_n(1,nObsDscTCF)        ///< z-scores for fit to female bycatch mortality in directed TCF 
    matrix zsDscMortBio_SCF_xn(1,nSXs,1,nObsDscSCF) ///< z-scores for fit to bycatch mortality, by sex, in SCF 
    matrix zsDscMortBio_RKF_xn(1,nSXs,1,nObsDscRKF) ///< z-scores for fit to bycatch mortality, by sex, in RKFCF 
    vector zsDscMortBio_GTF_n(1,nObsDscGTF)         ///< z-scores for fit to total bycatch mortality in GTF 
    number lkRetMortBio_TCFR          ///< likelihood for retained male catch mortality in directed TCF 
    number lkTotMortBio_TCFM          ///< likelihood for total male catch mortality in directed TCF
    number lkDscMortBio_TCFM          ///< likelihood for discard male catch mortality in directed TCF
    number lkDscMortBio_TCFF          ///< likelihood for female bycatch mortality in directed TCF 
    vector lkDscMortBio_SCF_x(1,nSXs) ///< likelihood for bycatch in SCF, by sex
    vector lkDscMortBio_RKF_x(1,nSXs) ///< likelihood for bycatch in RKF, by sex
    number lkDscMortBio_GTF           ///< likelihood for bycatch in GTF
    vector lkZCs(1,NUM_LEN_LIKE)      ///< likelihood for fishery size compositions
    
    matrix zsSrvMatBio_xn(1,nSXs,1,nObsSrvBio) ///< z-scores for fits to survey mature biomass
    vector lkSrvMatBio_x(1,nSXs)               ///< likelihood for survey mature biomass data
    
    number lkPrM2M ///< likelihood associated with the probability of molting to maturity
    
    //IMPORTANT CHANGE: was "endyr".  cannot be calculated in endyr. 
    3darray modSpNumMateTime_xsy(1,nSXs,1,nSCs,styr,endyr-1)
    matrix modSpNumMateTime_xy(1,nSXs,styr,endyr-1)
    matrix modSpBioMateTime_xy(1,nSXs,styr,endyr-1)
    vector mspbio_old_matetime(styr,endyr-1)  
    vector fspbio_new_matetime(styr,endyr-1)
    
    vector emspbio_matetime(styr,endyr-1)     // Spawning biomass at mating time stuff
    vector efspbio_new_matetime(styr,endyr-1) //"effective" female spawning biomass at mating time
    vector efspbio_matetime(styr,endyr-1)   
    
    //can be calculated in endyr:
    matrix modSpBioFishTime_xy(1,nSXs,styr,endyr)
    
    matrix modPrM2M(1,nSXs,1,nZBs)                              // Maturity-at-length
    
    // Outputs (not in the likelihood function)
    vector modFT_PopNumLegal_y(styr,endyr)
    vector modFT_PopBioLegal_y(styr,endyr)
    vector modTotBioMortLegal_TCFM_y(styr,endyr-1)     //(IMPORTANT CHANGE: used to be "endyr")
    vector modTotNumMortLegal_TCFM_y(styr,endyr-1)  //(IMPORTANT CHANGE: used to be "endyr")
    
    sdreport_number sdrDepletion
    
    sdreport_vector sdrNatMort_INF(styr,endyr-1);//natural mortality by year on immature new shell females
    sdreport_vector sdrNatMort_INM(styr,endyr-1);//natural mortality by year on immature new shell   males
    sdreport_vector sdrNatMort_MNF(styr,endyr-1);//natural mortality by year on   mature new shell females
    sdreport_vector sdrNatMort_MNM(styr,endyr-1);//natural mortality by year on   mature new shell   males
    sdreport_vector sdrNatMort_MOF(styr,endyr-1);//natural mortality by year on   mature old shell females
    sdreport_vector sdrNatMort_MOM(styr,endyr-1);//natural mortality by year on   mature old shell   males
    
    sdreport_vector sdrPrM2M_F(1,nZBs);//female prM2M
    sdreport_vector sdrPrM2M_M(1,nZBs);//  male prM2M
    sdreport_vector sdrMnGrw_F(1,nZBs);//mean female growth
    sdreport_vector sdrMnGrw_M(1,nZBs);//mean   male growth
      
    sdreport_vector sdrMMB(styr,endyr-1);   //MMB (at time of mating, so final year not included)
    sdreport_vector sdrMFB(styr,endyr-1);   //MFB (at time of mating, so final year not included)
    sdreport_vector sdrLnRec(styr,endyr);   //recruitment (NOTE: no lag to fertilization year!)
    
    objective_function_value f
    
//========================================================================
//========================================================================
PRELIMINARY_CALCS_SECTION

    //set initial probabilities for molt to maturity (20160324))
    //need to do this BEFORE setting initial parameter values
    modPrM2M(FEMALE) = 1.0;
    modPrM2M(FEMALE)(1,16) = obsAvgMatNS_xz(FEMALE)(1,16);
    modPrM2M(MALE)         = obsAvgMatNS_xz(  MALE);
    if (maturity_switch > 0) {
        // use logistic maturity curve for new shell males instead of fractions by year
        // this would be for initial population not probability of moving to mature
//        obsAvgMatNS_xz(MALE) = obsPrMatureM_z;
        modPrM2M(MALE) = obsPrMatureM_z;
    }
    
    if (usePin) {
        //need to recalculate modPrM2M values based on pin
        //--prMoltToMaturity(female|size)
        modPrM2M(FEMALE) = 1.0; //--prM2M(females> zBs(16)) assumed = 1
        if (optPrM2M==0){
            modPrM2M(FEMALE)(1,16) = mfexp(pPrM2MF);
            modPrM2M(MALE) = mfexp(pPrM2MM);
        } else {
            modPrM2M(FEMALE)(1,16) = 1.0/(1.0+mfexp(-pPrM2MF));
            modPrM2M(MALE) = 1.0/(1.0+mfexp(-pPrM2MM));
        }
    } else {
        //set initial values from control file inputs 
        //recruitment
        pMnLnRecInit = inpMnLnRecInit;
        pMnLnRec     = inpMnLnRec;
        
        //natural mortality multipliers
        pMfac_Imm  = inpMfac_Imm;
        pMfac_MatM = inpMfac_MatM;
        pMfac_MatF = inpMfac_MatF;
        pMfac_Big  = inpMfac_Big;
        
        //growth
        pGrAF1    = inpGrAF1;    // Female growth-increment
        pGrBF1    = inpGrBF1;    // Female growth-increment
        pGrAM1    = inpGrAM1;    // Male growth-increment
        pGrBM1    = inpGrBM1;    // Male growth-increment
        pGrBeta_x = inpGrBeta_x; // Growth beta--NOT estimated
        
        //set pr(molt-to-maturity) parameters based on input ogives
        if (optPrM2M==0){
            //2015 approach: do nothing and use default initialization
        } else if (optPrM2M==1){
            //2015 approach: parameters are logit-scale
            pPrM2MF = log(elem_div(modPrM2M(FEMALE)(1,16),1.0-modPrM2M(FEMALE)(1,16)));
            pPrM2MM = log(elem_div(modPrM2M(  MALE)      ,1.0-modPrM2M(  MALE)      ));
        }
        CheckFile<<"#--Maturity parameters"<<endl;
        CheckFile<<"# pPrM2MF = "<<endl<<pPrM2MF<<endl;
        CheckFile<<"# pPrM2MM = "<<endl<<pPrM2MM<<endl;
        
        //surveys
        //--catchabilities
        ////--males
        pSrv1_QM = inpSrv1_QM;
        pSrv2_QM = inpSrv2_QM;
        ////--females
        pSrv1_QF = inpSrv1_QF;
        pSrv2_QF = inpSrv2_QF;
        ////--selectivities
        //////--males
        pSrv1M_z50 = inpSrv1M_z50;
        pSrv1M_dz5095 = inpSrv1M_dz5095;
        //////--1982+
        pSrv2M_z50 = inpSrv2M_z50;
        pSrv2M_dz5095 = inpSrv2M_dz5095;
        //--females
        ////--1974-1981
        pSrv1F_z50 = inpSrv1F_z50;
        pSrv1F_dz5095= inpSrv1F_dz5095;    
        ////--1982+
        pSrv2F_z50 = inpSrv2F_z50;
        pSrv2F_dz5095 = inpSrv2F_dz5095;
        //fisheries
        //--log-scale fishing mortality or capture rates
        pAvgLnF_TCF = inpAvgLnF_TCF; // directed Tanner crab fishery
        pAvgLnF_SCF = inpAvgLnF_SCF; // snow crab fishery bycatch
        pAvgLnF_RKF = inpAvgLnF_RKF; // BBRKC fishery bycatch
        pAvgLnF_GTF = inpAvgLnF_GTF; // groundfish fisheries bycatch
        ////--log-scale offsets for females
        pAvgLnF_TCFF = inpAvgLnF_TCFF;  ///< female offset to ln-scale mean fishing mortality in directed fishery
        pAvgLnF_SCFF = inpAvgLnF_SCFF;  ///< female offset to ln-scale mean fishing mortality in snow crab fishery
        pAvgLnF_RKFF = inpAvgLnF_RKFF;  ///< female offset to ln-scale mean fishing mortality in BBRKC fishery
        pAvgLnF_GTFF = inpAvgLnF_GTFF;  ///< female offset to ln-scale mean fishing mortality in groundfish trawl fisheries
        //directed fishery selectivity and retention
        // Retention functions
        //-- styr-1990
        pRetTCFM_slpA1 = inpRetTCFM_slpA1;
        pRetTCFM_z50A1 = inpRetTCFM_z50A1;
        //-- 1991+  
        pRetTCFM_slpA2 = inpRetTCFM_slpA2;
        pRetTCFM_z50A2 = inpRetTCFM_z50A2;
        // Selectivity functions
        //--males, styr-1996
        pSelTCFM_slpA1 = inpSelTCFM_slpA1; 
        //--males, 2005+
        pSelTCFM_slpA2 = inpSelTCFM_slpA2;  
        pSelTCFM_mnLnZ50A2 = inpSelTCFM_mnLnZ50A2;
        //--females, styr+
        pSelTCFF_slp = inpSelTCFF_slp;
        pSelTCFF_z50 = inpSelTCFF_z50;
        // snow crab fishery bycatch selectivity
        //--males styr-1996
        pSelSCFM_slpA1 = inpSelSCFM_slpA1;
        pSelSCFM_z50A1 = inpSelSCFM_z50A1;
        pSelSCFM_slpD1 = inpSelSCFM_slpD1;
        pSelSCFM_lnZ50D1 = inpSelSCFM_lnZ50D1;
        //--males, 1997-2004
        pSelSCFM_slpA2 = inpSelSCFM_slpA2;
        pSelSCFM_z50A2 = inpSelSCFM_z50A2;
        pSelSCFM_slpD2 = inpSelSCFM_slpD2;
        pSelSCFM_lnZ50D2 = inpSelSCFM_lnZ50D2;
        //--males, 2005+
        pSelSCFM_slpA3 = inpSelSCFM_slpA3;
        pSelSCFM_z50A3 = inpSelSCFM_z50A3;
        pSelSCFM_slpD3 = inpSelSCFM_slpD3;
        pSelSCFM_lnZ50D3 = inpSelSCFM_lnZ50D3;
        //--females, styr-1996
        pSelSCFF_slpA1 = inpSelSCFF_slpA1;
        pSelSCFF_z50A1 = inpSelSCFF_z50A1;
        //--females, 1997-2004
        pSelSCFF_slpA2 = inpSelSCFF_slpA2;
        pSelSCFF_z50A2 = inpSelSCFF_z50A2;
        //--females, 2005+
        pSelSCFF_slpA3 = inpSelSCFF_slpA3;
        pSelSCFF_z50A3 = inpSelSCFF_z50A3;
        // BBRKC fishery bycatch selectivity
        //--males, styr-1996
        pSelRKFM_slpA1 = inpSelRKFM_slpA1;
        pSelRKFM_z50A1 = inpSelRKFM_z50A1;
        //--males, 1997-2004
        pSelRKFM_slpA2 = inpSelRKFM_slpA2;
        pSelRKFM_z50A2 = inpSelRKFM_z50A2;
        //--males, 2005+
        pSelRKFM_slpA3 = inpSelRKFM_slpA3;
        pSelRKFM_z50A3 = inpSelRKFM_z50A3;
        //--females, styr-1996
        pSelRKFF_slpA1 = inpSelRKFF_slpA1;
        pSelRKFF_z50A1 = inpSelRKFF_z50A1;
        //--females, 1997-2004
        pSelRKFF_slpA2 = inpSelRKFF_slpA2;
        pSelRKFF_z50A2 = inpSelRKFF_z50A2;
        //--females, 2005+
        pSelRKFF_slpA3 = inpSelRKFF_slpA3;
        pSelRKFF_z50A3 = inpSelRKFF_z50A3;
        // groundfish fisheries bycatch selectivity 
        //--males, 1973-1987
        pSelGTFM_slpA1 = inpSelGTFM_slpA1;
        pSelGTFM_z50A1 = inpSelGTFM_z50A1;
        //--males, 1988-1996
        pSelGTFM_slpA2 = inpSelGTFM_slpA2;
        pSelGTFM_z50A2 = inpSelGTFM_z50A2;
        //--males, 1997+
        pSelGTFM_slpA3 = inpSelGTFM_slpA3;
        pSelGTFM_z50A3 = inpSelGTFM_z50A3;
        //--females, 1973-1987
        pSelGTFF_slpA1 = inpSelGTFF_slpA1;
        pSelGTFF_z50A1 = inpSelGTFF_z50A1;
        //--females, 1988-1996
        pSelGTFF_slpA2 = inpSelGTFF_slpA2;
        pSelGTFF_z50A2 = inpSelGTFF_z50A2;
        //females, 1997+
        pSelGTFF_slpA3 = inpSelGTFF_slpA3;
        pSelGTFF_z50A3 = inpSelGTFF_z50A3;

        //effort extrapolation parameters
        pLnEffXtr_TCF = inpLnEffXtr_TCF;
        pLnEffXtr_SCF = inpLnEffXtr_SCF;
        pLnEffXtr_RKF = inpLnEffXtr_RKF;
        pLnEffXtr_GTF = inpLnEffXtr_GTF;
    }
    if (jitter) jitterParameters(ptrMC->jitFrac);
 
//     CheckFile<<"catch_disc(1) "<<catch_disc(1)<<endl;
//     CheckFile<<"catch_disc(2) "<<catch_disc(2)<<endl;
    CheckFile<<"catch ret numbers (millions)"<<endl<<obsRetCatchNum_y<<endl;
    
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
    
    // Survey data  
    d5_array tmpObsPrNatZ_Srv_msxnz(1,nMSs,1,nSCs,1,nSXs,1,nObsZCsSrv,1,nZBs);    // temporary array for survey length frequency by maturity state, shell condition, sex
    tmpObsPrNatZ_Srv_msxnz.initialize();
    for (int m=1; m<=nMSs; m++) { //maturity
        for (int s=1; s<=nSCs; s++) { //shell condition
            for (int x=1; x<=nSXs;x++) { //sex
                for (int n=1; n <= nObsZCsSrv; n++){ //year counter
                    tmpObsPrNatZ_Srv_msxnz(m,s,x,n) = obsSrvNatZs_msxnz(m,s,x,n)/obsTotSrvNum_n(n);
                } //year counter
            } //sex
        } //shell condition
    } //maturity
    
    // use logistic maturity curve for new and old shell male survey data if switch>0 instead of yearly samples
    // old shell already uses ok maturity curve (AEP only applies to OLD SHELL?)
    if (maturity_switch > 0){
        for(int i=1; i <= nObsZCsSrv; i++){
//             tmps = (tmpObsPrNatZ_Srv_msxnz(1,2,2,i)+tmpObsPrNatZ_Srv_msxnz(2,2,2,i));
//             tmpObsPrNatZ_Srv_msxnz(2,2,2,i) = elem_prod(obsAvgMatOS_xz(2),tmps);
//             tmpObsPrNatZ_Srv_msxnz(1,2,2,i) = elem_prod(1.0-obsAvgMatOS_xz(2),tmps);
            dvector tmps = (tmpObsPrNatZ_Srv_msxnz(IMMATURE,OLD_SHELL,MALE,i)+tmpObsPrNatZ_Srv_msxnz(MATURE,OLD_SHELL,MALE,i));
            tmpObsPrNatZ_Srv_msxnz(  MATURE,OLD_SHELL,MALE,i) = elem_prod(obsAvgMatOS_xz(MALE),tmps);
            tmpObsPrNatZ_Srv_msxnz(IMMATURE,OLD_SHELL,MALE,i) = elem_prod(1.0-obsAvgMatOS_xz(MALE),tmps);
        }
    }
    
    // Store results
    obsSrvPrNatZ_msxnz(IMMATURE) = tmpObsPrNatZ_Srv_msxnz(IMMATURE);
    obsSrvPrNatZ_msxnz(  MATURE) = tmpObsPrNatZ_Srv_msxnz(  MATURE);
    CheckFile<<"obsSrvPrNatZ_msxnz(IMMATURE)"<<endl<<obsSrvPrNatZ_msxnz(IMMATURE)<<endl;
    CheckFile<<"obsSrvPrNatZ_msxnz(  MATURE)"<<endl<<obsSrvPrNatZ_msxnz(  MATURE)<<endl;
    
    // for maturity and shell condition together in survey length comp fits
    obsSrvPrNatZ_mxnz.initialize();
    {
        int sex;
        for (int i=1; i <= nObsZCsSrv; i++){
            sex = MALE;
            obsSrvPrNatZ_mxnz(IMMATURE,sex,i) = obsSrvPrNatZ_msxnz(IMMATURE,NEW_SHELL,sex,i)+obsSrvPrNatZ_msxnz(IMMATURE,OLD_SHELL,sex,i);
            obsSrvPrNatZ_mxnz(  MATURE,sex,i) = obsSrvPrNatZ_msxnz(  MATURE,NEW_SHELL,sex,i)+obsSrvPrNatZ_msxnz(  MATURE,OLD_SHELL,sex,i);
            offset(9)  -= ssObsZCsSrv_msxn(MATURE,NEW_SHELL,sex,i)*obsSrvPrNatZ_mxnz(IMMATURE,sex,i)*log(obsSrvPrNatZ_mxnz(IMMATURE,sex,i)+p_const);//DON'T THINK CORRECT NSAMPLES IS BEING APPLIED HERE!
            offset(10) -= ssObsZCsSrv_msxn(MATURE,OLD_SHELL,sex,i)*obsSrvPrNatZ_mxnz(  MATURE,sex,i)*log(obsSrvPrNatZ_mxnz(  MATURE,sex,i)+p_const);
            sex = FEMALE;
            obsSrvPrNatZ_mxnz(IMMATURE,sex,i) = obsSrvPrNatZ_msxnz(IMMATURE,NEW_SHELL,sex,i)+obsSrvPrNatZ_msxnz(IMMATURE,OLD_SHELL,sex,i);
            obsSrvPrNatZ_mxnz(  MATURE,sex,i) = obsSrvPrNatZ_msxnz(  MATURE,NEW_SHELL,sex,i)+obsSrvPrNatZ_msxnz(  MATURE,OLD_SHELL,sex,i);
            offset(11) -= ssObsZCsSrv_msxn(MATURE,NEW_SHELL,sex,i)*obsSrvPrNatZ_mxnz(IMMATURE,sex,i)*log(obsSrvPrNatZ_mxnz(IMMATURE,sex,i)+p_const);
            offset(12) -= ssObsZCsSrv_msxn(MATURE,OLD_SHELL,sex,i)*obsSrvPrNatZ_mxnz(  MATURE,sex,i)*log(obsSrvPrNatZ_mxnz(  MATURE,sex,i)+p_const);
        }
    }
    CheckFile<<"offset( 9) = "<<offset( 9)<< endl;  
    CheckFile<<"offset(10) = "<<offset(10)<< endl;  
    CheckFile<<"offset(11) = "<<offset(11)<< endl;  
    CheckFile<<"offset(12) = "<<offset(12)<< endl;  
                
    
    //survey numbers from aggregate input data
    obsSrvNum_y.initialize();//wts: now initializing this
    for(int i=1;i<=nObsSrvBio;i++) obsSrvNum_y(yrsObsSrvBio_n(i)) = obsSrvNum_n(i);//<-wts : yrsObsZCsSrv_n(i) used to index obsSrvNum_y below [so yrsObsSrvBio_n=yrsObsZCsSrv_n?]
    CheckFile<<"obsSrvNum_y"<<endl<<obsSrvNum_y<<endl;
    
    // Compute survey biomass
    obsSrvNum_xyz.initialize();
    obsSrvImmNum_sxy.initialize();
    obsSrvMatNum_sxy.initialize();
    obsSrvBio_y.initialize();
    obsSrvBio_xy.initialize();
    obsSrvImmBio_xy.initialize();
    obsSrvMatBio_xy.initialize();
    for (int m=1;m<=2;m++) { //maturity status
        for (int s=1;s<=2;s++) { //s condition
            for (int x=1;x<=2;x++) { //x
                for (int i=1; i <= nObsZCsSrv; i++) {
                    obsSrvNum_xyz(x,yrsObsZCsSrv_n(i)) += obsSrvPrNatZ_msxnz(m,s,x,i)*obsSrvNum_y(yrsObsZCsSrv_n(i));
                    obsSrvBio_xy(x,yrsObsZCsSrv_n(i))  += obsSrvPrNatZ_msxnz(m,s,x,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wt_xmz(x,m);
                    obsSrvBio_y(yrsObsZCsSrv_n(i))     += obsSrvPrNatZ_msxnz(m,s,x,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wt_xmz(x,m);
                    //  sum to get mature biomass by x (AEP index is mature animals only?)
                    if(m==MATURE) {
                        obsSrvMatNum_sxy(s,x,yrsObsZCsSrv_n(i)) += sum(obsSrvPrNatZ_msxnz(m,s,x,i)*obsSrvNum_y(yrsObsZCsSrv_n(i)));
                        obsSrvMatBio_xy(x,yrsObsZCsSrv_n(i))    +=     obsSrvPrNatZ_msxnz(m,s,x,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wt_xmz(x,m);
                    } else {
                        obsSrvImmNum_sxy(s,x,yrsObsZCsSrv_n(i)) += sum(obsSrvPrNatZ_msxnz(m,s,x,i)*obsSrvNum_y(yrsObsZCsSrv_n(i)));
                        obsSrvImmBio_xy(x,yrsObsZCsSrv_n(i))    +=     obsSrvPrNatZ_msxnz(m,s,x,i)*obsSrvNum_y(yrsObsZCsSrv_n(i))*wt_xmz(x,m);
                    }
                }
            }
        }
    }
    CheckFile<<"obsSrvNum_xyz"   <<endl<<obsSrvNum_xyz   <<endl;
    CheckFile<<"obsSrvImmNum_sxy"<<endl<<obsSrvImmNum_sxy<<endl;
    CheckFile<<"obsSrvMatNum_sxy"<<endl<<obsSrvMatNum_sxy<<endl;
    CheckFile<<"obsSrvBio_y"     <<endl<<obsSrvBio_y     <<endl;
    CheckFile<<"obsSrvBio_xy"    <<endl<<obsSrvBio_xy    <<endl;
    CheckFile<<"obsSrvImmBio_xy" <<endl<<obsSrvImmBio_xy <<endl;
    CheckFile<<"obsSrvMatBio_xy" <<endl<<obsSrvMatBio_xy <<endl;
    
    // Number of large males
    obsSrvNumLegal_n.initialize();
    obsSrvBioLegal_n.initialize();
    for(int i=1;i<=nObsZCsSrv;i++) {
//        // take 1/2 of the 100-104 bin, 
//        obsSrvNumLegal_n(i) = 0.5*obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i),23);            //<--hardwired index
//        obsSrvBioLegal_n(i) = obsSrvNumLegal_n(i)*wt_xmz(MALE,  MATURE)(23);                            //<--hardwired index
//        for(int j=24;j<=nZBs;j++) {                                          //<--hardwired index
//            obsSrvNumLegal_n(i) += obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i),j);
//            obsSrvBioLegal_n(i) += obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i),j)*wt_xmz(MALE,  MATURE)(j);
//        }
        obsSrvNumLegal_n(i) += sum(obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i))(iZLegal,nZBs));
        obsSrvBioLegal_n(i) += obsSrvNum_xyz(MALE,yrsObsZCsSrv_n(i))(iZLegal,nZBs)*wt_xmz(MALE,  MATURE)(iZLegal,nZBs);
    }
    CheckFile<<"obsSrvNumLegal_n (millions)"<<endl<<obsSrvNumLegal_n<<endl;
    CheckFile<<"obsSrvBioLegal_n (1000's t)"<<endl<<obsSrvBioLegal_n<<endl;
    
    //mean effort in bycatch fisheries over years with bycatch data
    mnEff_SCF = mean(effSCF_y(yrsObsDscSCF));
    mnEff_RKF = mean(effRKF_y(yrsObsDscRKF));

    //make sure these calculations get done at least once!!
    // Compute the moulting probabilities
    get_moltingp();
    //  cout<<"done moltingp"<<endl;
    // estimate growth function
    get_growth1();//only option now
    //  cout<<"done growth"<<endl;
    // Set maturity
    get_maturity(); //should do nothing because parameter estimation not yet active
    
    //run population mode with initial parameter values
    runPopMod();
    
    //evaluate the objective function for initial parameter values
    evaluate_the_objective_function();        
    
    //write reports for initial model configuration
    ofstream initReptToR("TCSAM2013.NEWSTYLE.init.R");
    writeToR_NEW(initReptToR);
    initReptToR.close();
    ofstream initReptToR1("TCSAM2013.OLDSTYLE.init.R");
    writeToR_OLD(initReptToR1);
    initReptToR1.close();
    ofstream initLLsToCSV("TCSAM2013.init_likelihood_components.csv");
    writeLikelihoodComponents(initLLsToCSV,0);
    initLLsToCSV.close();
    
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        openMCMCFile();
        cout<<"MCEVAL is on"<<endl;
    }
    
    CheckFile<<"------------------------------------------------------------------------"<<endl<<endl;
    CheckFile<<"Initial Parameter Settings"<<endl;
    writeParameters(CheckFile,0,0);
    ofstream os1("TCSAM2013.params.all.init.csv");
    writeParameters(os1,0,0);
    os1.close();
    ofstream os2("TCSAM2013.params.active.init.csv");
    writeParameters(os2,0,1);
    os2.close();

        if (doOFL&&debugOFL){
            cout<<"Test OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.init.txt", ios::trunc);
            echoOFL<<"----Testing calcOFL()"<<endl;
            calcOFL(asmtYr,debugOFL,echoOFL);//updates oflResults
            oflResults.writeCSVHeader(echoOFL); echoOFL<<endl;
            oflResults.writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL<<"----Finished testing calcOFL()!"<<endl;
            echoOFL.close();
            cout<<"Finished testing OFL calculations!"<<endl;
        }

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
//    writeToR_OLD(R_out);             //for testing
    
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
//     writeToR_NEW(myRout);
//     myRout.close();
//     exit(1);
    
    //for testing
//     CheckFile<<"fTCF_xy       = "<<endl<<tb<<fTCF_xy<<endl;
//     CheckFile<<"fmortdf     = "<<endl<<tb<<fmortdf<<endl;
//     CheckFile<<"fSCF_xy = "<<endl<<tb<<fSCF_xy<<endl;
//     CheckFile<<"fRKF_xy   = "<<endl<<tb<<fRKF_xy<<endl;
//     CheckFile<<"fGTF_xy      = "<<endl<<tb<<fGTF_xy<<endl;
//     CheckFile<<"selTCFM_syz(NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<selTCFM_syz(NEW_SHELL,iy)<<endl;
//     CheckFile<<"selTCFM_syz(OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<selTCFM_syz(OLD_SHELL,iy)<<endl;

//     CheckFile<<"retFcn_syz(NEW_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<retFcn_syz(NEW_SHELL,iy)<<endl;
//     CheckFile<<"retFcn_syz(OLD_SHELL) = "<<endl;
//     for (int iy=styr;iy<=endyr;iy++) CheckFile<<iy<<tb<<retFcn_syz(OLD_SHELL,iy)<<endl;
//     CheckFile<<"selSCF_cxz(1-3,FEMALE) = "<<endl<<tb<<selSCF_cxz(1,FEMALE)<<endl<<tb<<selSCF_cxz(2,FEMALE)<<endl<<tb<<selSCF_cxz(3,FEMALE)<<endl;
//     CheckFile<<"selSCF_cxz(1-3,  MALE) = "<<endl<<tb<<selSCF_cxz(1,  MALE)<<endl<<tb<<selSCF_cxz(2,  MALE)<<endl<<tb<<selSCF_cxz(3,  MALE)<<endl;
//     CheckFile<<"selRKF_cxz(1-3,FEMALE) = "<<endl<<tb<<selRKF_cxz(1,FEMALE)<<endl<<tb<<selRKF_cxz(2,FEMALE)<<endl<<tb<<selRKF_cxz(3,FEMALE)<<endl;
//     CheckFile<<"selRKF_cxz(1-3,  MALE) = "<<endl<<tb<<selRKF_cxz(1,  MALE)<<endl<<tb<<selRKF_cxz(2,  MALE)<<endl<<tb<<selRKF_cxz(3,  MALE)<<endl;
//     CheckFile<<"selGTF_cxz(1-3,FEMALE) = "<<endl<<tb<<selGTF_cxz(1,FEMALE)<<endl<<tb<<selGTF_cxz(2,FEMALE)<<endl<<tb<<selGTF_cxz(3,FEMALE)<<endl;
//     CheckFile<<"selGTF_cxz(1-3,  MALE) = "<<endl<<tb<<selGTF_cxz(1,  MALE)<<endl<<tb<<selGTF_cxz(2,  MALE)<<endl<<tb<<selGTF_cxz(3,  MALE)<<endl;

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
    // growth matrix initially calculated in prelimn calcs. 
    // recalculate only if growth parameters are estimated
    if(active(pGrAM1) || active(pGrBM1) || active(pGrAF1) || active(pGrBF1) || active(pGrBeta_x)) {
        get_growth1();//only option now
//        cout<<" growth "<<endl;
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
FUNCTION void writeParameters(ofstream& os,int toR, int willBeActive)
//NOTE: use '\\n' in labels to insert a newline character in a plot title

    os<<"index, phase, idx.mn, idx.mx, min, max, value, name, type, category, process, label, description"<<endl;
    wts::writeParameter(os,pMnLnRec,toR,willBeActive,    "population, recruitment, log-scale mean, log-scale mean");      
    wts::writeParameter(os,pMnLnRecInit,toR,willBeActive,"population, recruitment, initial log-scale mean, initial log-scale mean"); 
    wts::writeParameter(os,pRecAlpha,toR,willBeActive,   "population, recruitment, alpha, size distribution alpha parameter");      
    wts::writeParameter(os,pRecBeta,toR,willBeActive,    "population, recruitment, beta, size distribution beta parameter");      
    wts::writeParameter(os,pRecDevs,toR,willBeActive,    "population, recruitment devs, dev, log-scale deviation");      
    wts::writeParameter(os,pRecDevsHist,toR,willBeActive,"population, initial recruitment devs, dev, log-scale deviation"); 
    
    wts::writeParameter(os,pMfac_Imm,toR,willBeActive, "population, natural mortality multipliers, immature\\ncrab, multiplier for immature crab");      
    wts::writeParameter(os,pMfac_MatM,toR,willBeActive,"population, natural mortality multipliers, mature\\nmales, multiplier for mature male crab");      
    wts::writeParameter(os,pMfac_MatF,toR,willBeActive,"population, natural mortality multipliers, mature\\nfemales, multiplier for mature female crab");      
    wts::writeParameter(os,pMfac_Big,toR,willBeActive, "population, natural mortality multipliers, mature crab\\n1980-1984, multiplier for 1980-1984");      
    
    wts::writeParameter(os,pPrM2MF,toR,willBeActive,"population, molt-to-maturity: females, bin, female"); 
    wts::writeParameter(os,pPrM2MM,toR,willBeActive,"population, molt-to-maturity: males, bin, male"); 
    
    wts::writeParameter(os,pGrAF1,toR,willBeActive,"population, growth, female a, female mean growth a parameter");      
    wts::writeParameter(os,pGrBF1,toR,willBeActive,"population, growth, female b, female mean growth b parameter");      
    wts::writeParameter(os,pGrAM1,toR,willBeActive,"population, growth, male a, male mean growth a parameter");      
    wts::writeParameter(os,pGrBM1,toR,willBeActive,"population, growth, male b, male mean growth b parameter");      
    wts::writeParameter(os,pGrBeta_x,toR,willBeActive,"population, growth, beta, size transition beta parameter");  
    
    wts::writeParameter(os,pSrv1_QM,toR,willBeActive,"surveys, surveys, -1981 \\nmale Q, males [-1981]");       
    wts::writeParameter(os,pSrv2_QM,toR,willBeActive,"surveys, surveys, 1982+ \\nmale Q, males [1982+]");       
    wts::writeParameter(os,pSrv1_QF,toR,willBeActive,"surveys, surveys, -1981 \\nfemale Q, females [-1981]");      
    wts::writeParameter(os,pSrv2_QF,toR,willBeActive,"surveys, surveys, 1982+ \\nfemale Q, females [1982+]"); 
    
    wts::writeParameter(os,pSrv1M_z50,toR,willBeActive,   "surveys, survey selectivity, -1981 \\nmale z50,     male size at 50%-selected [-1981]");   
    wts::writeParameter(os,pSrv1M_dz5095,toR,willBeActive,"surveys, survey selectivity, -1981 \\nmale z95-z50, male offset to 95%-selected [-1981]"); 
    wts::writeParameter(os,pSrv2M_z50,toR,willBeActive,   "surveys, survey selectivity, 1982+ \\nmale z50,     male size at 50%-selected [1982+]");   
    wts::writeParameter(os,pSrv2M_dz5095,toR,willBeActive,"surveys, survey selectivity, 1982+ \\nmale z95-z50, male offset to 95%-selected [1982+]"); 
    
    wts::writeParameter(os,pSrv1F_z50,toR,willBeActive,   "surveys, survey selectivity, -1981 \\nfemale z50, female size at 50%-selected [-1981]");   
    wts::writeParameter(os,pSrv1F_dz5095,toR,willBeActive,"surveys, survey selectivity, -1981 \\nfemale z95-z50, female offset to 95%-selected [-1981]"); 
    wts::writeParameter(os,pSrv2F_z50,toR,willBeActive,   "surveys, survey selectivity, 1982+ \\nfemale z50, female size at 50%-selected [1982+]");
    wts::writeParameter(os,pSrv2F_dz5095,toR,willBeActive,"surveys, survey selectivity, 1982+ \\nfemale z95-z50, female offset to 95%-selected [1982+]"); 
    
    wts::writeParameter(os,pAvgLnF_TCF,toR,willBeActive,"fisheries, mortality/capture rate, TCF 1965+ \\nln-scale mean, TCF ln-scale mean [1965+]");   
    wts::writeParameter(os,pAvgLnF_SCF,toR,willBeActive,"fisheries, mortality/capture rate, SCF 1992+ \\nln-scale mean, SCF ln-scale mean [1992+]");   
    wts::writeParameter(os,pAvgLnF_RKF,toR,willBeActive,"fisheries, mortality/capture rate, RKF 1992+ \\nln-scale mean, RKF ln-scale mean [1992+]");   
    wts::writeParameter(os,pAvgLnF_GTF,toR,willBeActive,"fisheries, mortality/capture rate, GTF 1973+ \\nln-scale mean, GTF ln-scale mean [1973+]");   
    
    wts::writeParameter(os,pAvgLnF_TCFF,toR,willBeActive,"fisheries, mortality/capture rate, TCF female offset, TCF ln-scale female offset");
    wts::writeParameter(os,pAvgLnF_SCFF,toR,willBeActive,"fisheries, mortality/capture rate, SCF female offset, SCF ln-scale female offset");
    wts::writeParameter(os,pAvgLnF_RKFF,toR,willBeActive,"fisheries, mortality/capture rate, RKF female offset, RKF ln-scale female offset");
    wts::writeParameter(os,pAvgLnF_GTFF,toR,willBeActive,"fisheries, mortality/capture rate, GTF female offset, GTF ln-scale female offset");
    
    wts::writeParameter(os,pF_DevsTCF,toR,willBeActive,"fisheries, TCF mortality/capture rate devs, , ln-scale devs [1965+]");    
    wts::writeParameter(os,pF_DevsSCF,toR,willBeActive,"fisheries, SCF mortality/capture rate devs, , ln-scale devs [1992+]");    
    wts::writeParameter(os,pF_DevsRKF,toR,willBeActive,"fisheries, RKF mortality/capture rate devs, , ln-scale devs [1992+]");    
    wts::writeParameter(os,pF_DevsGTF,toR,willBeActive,"fisheries, GTF mortality/capture rate devs, , ln-scale devs [1973+]");    
    
    wts::writeParameter(os,pRetTCFM_slpA1,toR,willBeActive,"fisheries, TCF retention, -1990\\n slope, slope [-1990]");   
    wts::writeParameter(os,pRetTCFM_z50A1,toR,willBeActive,"fisheries, TCF retention, -1990\\n z50, size at 50%-selected [-1990]");    
    wts::writeParameter(os,pRetTCFM_slpA2,toR,willBeActive,"fisheries, TCF retention, 1991+\\n slope, slope [1991+]");   
    wts::writeParameter(os,pRetTCFM_z50A2,toR,willBeActive,"fisheries, TCF retention, 1991+\\n z50,size at 50%-selected [1991+]");    
    
    wts::writeParameter(os,pSelTCFM_slpA1,toR,willBeActive,    "fisheries, TCF selectivity, males (-1996)\\n slope, male slope [-1996]");   
    wts::writeParameter(os,pSelTCFM_slpA2,toR,willBeActive,    "fisheries, TCF selectivity, males (1997+)\\n slope, male slope [1997+]");   
    wts::writeParameter(os,pSelTCFM_mnLnZ50A2,toR,willBeActive,"fisheries, TCF selectivity, males \\n mean(ln(z50)), male ln-scale mean size at 50%-selected");    
    wts::writeParameter(os,pSelTCFM_devsZ50,toR,willBeActive,  "fisheries, TCF selectivity, males z50\\n dev, male ln-scale devs in size at 50%-selected [1991+]");    
    
    wts::writeParameter(os,pSelTCFF_slp,toR,willBeActive,"fisheries, TCF selectivity, female slope\\n all years, female slope [all years]");   
    wts::writeParameter(os,pSelTCFF_z50,toR,willBeActive,"fisheries, TCF selectivity, female z50\\n all years, female size at 50%-selected [all years]");    
    
    wts::writeParameter(os,pSelSCFF_slpA1,toR,willBeActive,"fisheries, SCF selectivity, females (-1996)\\n  slope, female slope [-1996]");   
    wts::writeParameter(os,pSelSCFF_z50A1,toR,willBeActive,"fisheries, SCF selectivity, females (-1996)\\n  z50, female size at 50%-selected [-1996]");    
    wts::writeParameter(os,pSelSCFF_slpA2,toR,willBeActive,"fisheries, SCF selectivity, females (1997-2004)\\n  slope, female slope [1997-2004]");   
    wts::writeParameter(os,pSelSCFF_z50A2,toR,willBeActive,"fisheries, SCF selectivity, females (1997-2004)\\n  z50, female size at 50%-selected [1997-2004]");    
    wts::writeParameter(os,pSelSCFF_slpA3,toR,willBeActive,"fisheries, SCF selectivity, females (2005+)\\n  slope, female slope [2005+]");   
    wts::writeParameter(os,pSelSCFF_z50A3,toR,willBeActive,"fisheries, SCF selectivity, females (2005+)\\n  z50, female size at 50%-selected [2005+]");    
    
    wts::writeParameter(os,pSelSCFM_slpA1,toR,willBeActive,  "fisheries, SCF selectivity, males (-1996)\\n  asc. slope, male ascending slope [-1996]");   
    wts::writeParameter(os,pSelSCFM_z50A1,toR,willBeActive,  "fisheries, SCF selectivity, males (-1996)\\n  asc. z50, male ascending size at 50%-selected [-1996]");    
    wts::writeParameter(os,pSelSCFM_slpD1,toR,willBeActive,  "fisheries, SCF selectivity, males (-1996)\\n  dsc. slope, male descending slope [-1996]");   
    wts::writeParameter(os,pSelSCFM_lnZ50D1,toR,willBeActive,"fisheries, SCF selectivity, males (-1996)\\n  dsc. z50 offset, male descending ln-scale offset to size at 50%-selected [-1996]");    
    
    wts::writeParameter(os,pSelSCFM_slpA2,toR,willBeActive,  "fisheries, SCF selectivity, males (1997-2004)\\n  asc. slope, male ascending slope [1997-2004]");   
    wts::writeParameter(os,pSelSCFM_z50A2,toR,willBeActive,  "fisheries, SCF selectivity, males (1997-2004)\\n  asc. z50, male ascending size at 50%-selected [1997-2004]");    
    wts::writeParameter(os,pSelSCFM_slpD2,toR,willBeActive,  "fisheries, SCF selectivity, males (1997-2004)\\n  dsc. slope, male descending slope [1997-2004]");   
    wts::writeParameter(os,pSelSCFM_lnZ50D2,toR,willBeActive,"fisheries, SCF selectivity, males (1997-2004)\\n  dsc. z50 offset, male descending ln-scale offset to size at 50%-selected [1997-2004]");    
    
    wts::writeParameter(os,pSelSCFM_slpA3,toR,willBeActive,  "fisheries, SCF selectivity, males (2005+)\\n  asc. slope, male ascending slope [2005+]");   
    wts::writeParameter(os,pSelSCFM_z50A3,toR,willBeActive,  "fisheries, SCF selectivity, males (2005+)\\n  asc. z50, male ascending size at 50%-selected [2005+]");    
    wts::writeParameter(os,pSelSCFM_slpD3,toR,willBeActive,  "fisheries, SCF selectivity, males (2005+)\\n  dsc. slope, male descending slope [2005+]");   
    wts::writeParameter(os,pSelSCFM_lnZ50D3,toR,willBeActive,"fisheries, SCF selectivity, males (2005+)\\n  dsc. z50 offset, male descending ln-scale offset to size at 50%-selected [2005+]");    
    
    wts::writeParameter(os,pSelRKFF_slpA1,toR,willBeActive,"fisheries, RKF selectivity, females (-1996)\\n  slope, female slope [-1996]");   
    wts::writeParameter(os,pSelRKFF_z50A1,toR,willBeActive,"fisheries, RKF selectivity, females (-1996)\\n  z50, female size at 50%-selected [-1996]");    
    wts::writeParameter(os,pSelRKFF_slpA2,toR,willBeActive,"fisheries, RKF selectivity, females (1997-2004)\\n  slope, female slope [1997-2004]");   
    wts::writeParameter(os,pSelRKFF_z50A2,toR,willBeActive,"fisheries, RKF selectivity, females (1997-2004)\\n  z50, female size at 50%-selected [1997-2004]");    
    wts::writeParameter(os,pSelRKFF_slpA3,toR,willBeActive,"fisheries, RKF selectivity, females (2005+)\\n  slope, female slope [2005+]");   
    wts::writeParameter(os,pSelRKFF_z50A3,toR,willBeActive,"fisheries, RKF selectivity, females (2005+)\\n  z50, female size at 50%-selected [2005+]");    
    
    wts::writeParameter(os,pSelRKFM_slpA1,toR,willBeActive,"fisheries, RKF selectivity, males (-1996)\\n  slope, male slope [-1996]");   
    wts::writeParameter(os,pSelRKFM_z50A1,toR,willBeActive,"fisheries, RKF selectivity, males (-1996)\\n  z50, male size at 50%-selected [-1996]");    
    wts::writeParameter(os,pSelRKFM_slpA2,toR,willBeActive,"fisheries, RKF selectivity, males (1997-2004)\\n  slope, male slope [1997-2004]");   
    wts::writeParameter(os,pSelRKFM_z50A2,toR,willBeActive,"fisheries, RKF selectivity, males (1997-2004)\\n  z50, male size at 50%-selected [1997-2004]");    
    wts::writeParameter(os,pSelRKFM_slpA3,toR,willBeActive,"fisheries, RKF selectivity, males (2005+)\\n  slope, male slope [2005+]");   
    wts::writeParameter(os,pSelRKFM_z50A3,toR,willBeActive,"fisheries, RKF selectivity, males (2005+)\\n  z50, male size at 50%-selected [2005+]");    
    
    wts::writeParameter(os,pSelGTFF_slpA1,toR,willBeActive,"fisheries, GTF selectivity, females (-1987)\\n  slope, female slope [-1987]");   
    wts::writeParameter(os,pSelGTFF_z50A1,toR,willBeActive,"fisheries, GTF selectivity, females (-1987)\\n  z50, female size at 50%-selected [-1987]");   
    wts::writeParameter(os,pSelGTFF_slpA2,toR,willBeActive,"fisheries, GTF selectivity, females (1988-1996)\\n  slope, female slope [1988-1996]");   
    wts::writeParameter(os,pSelGTFF_z50A2,toR,willBeActive,"fisheries, GTF selectivity, females (1988-1996)\\n  z50, female size at 50%-selected [1988-1996]");   
    wts::writeParameter(os,pSelGTFF_slpA3,toR,willBeActive,"fisheries, GTF selectivity, females (1997+)\\n  slope, female slope [1997+]");   
    wts::writeParameter(os,pSelGTFF_z50A3,toR,willBeActive,"fisheries, GTF selectivity, females (1997+)\\n  z50, female size at 50%-selected [1997+]");   
    
    wts::writeParameter(os,pSelGTFM_slpA1,toR,willBeActive,"fisheries, GTF selectivity, males (-1987)\\n  slope, male slope [-1987]");   
    wts::writeParameter(os,pSelGTFM_z50A1,toR,willBeActive,"fisheries, GTF selectivity, males (-1987)\\n  z50, male size at 50%-selected [-1987]");   
    wts::writeParameter(os,pSelGTFM_slpA2,toR,willBeActive,"fisheries, GTF selectivity, males (1988-1996)\\n  slope, male slope [1988-1996]");   
    wts::writeParameter(os,pSelGTFM_z50A2,toR,willBeActive,"fisheries, GTF selectivity, males (1988-1996)\\n  z50, male size at 50%-selected [1988-1996]");   
    wts::writeParameter(os,pSelGTFM_slpA3,toR,willBeActive,"fisheries, GTF selectivity, males (1997+)\\n  slope, male slope [1997+]");   
    wts::writeParameter(os,pSelGTFM_z50A3,toR,willBeActive,"fisheries, GTF selectivity, males (1997+)\\n  z50, male size at 50%-selected [1997+]");   
    
    wts::writeParameter(os,pLnEffXtr_TCF,toR,willBeActive,"fisheries, mortality/capture rate, TCF lnQ, TCF effort extrapolation");      
    wts::writeParameter(os,pLnEffXtr_SCF,toR,willBeActive,"fisheries, mortality/capture rate, SCF lnQ, SCF effort extrapolation");      
    wts::writeParameter(os,pLnEffXtr_RKF,toR,willBeActive,"fisheries, mortality/capture rate, RKF lnQ, RKF effort extrapolation");      
    wts::writeParameter(os,pLnEffXtr_GTF,toR,willBeActive,"fisheries, mortality/capture rate, GTF lnQ, GTF effort extrapolation");      
    
// ----------------------------------------------------------------------
FUNCTION void jitterParameters(double fac)   //wts: new 2014-05-10
    cout<<"starting jitterParameters"<<endl;
    
    pGrAF1 = wts::jitterParameter(pGrAF1,fac,rng);     // Female growth-increment
    pGrBF1 = wts::jitterParameter(pGrBF1,fac,rng);     // Female growth-increment
    pGrAM1 = wts::jitterParameter(pGrAM1,fac,rng);     // Male growth-increment
    pGrBM1 = wts::jitterParameter(pGrBM1,fac,rng);     // Male growth-increment
    
    pGrBeta_x = wts::jitterParameter(pGrBeta_x,fac,rng); // Growth beta                                //this is NOT estimated (why?)
    pMfac_Imm  = wts::jitterParameter(pMfac_Imm,fac,rng);  // natural mortality multiplier for immature females and males
    pMfac_MatM      = wts::jitterParameter(pMfac_MatM,fac,rng);      // natural mortality multiplier for mature new and old shell male
    pMfac_MatF      = wts::jitterParameter(pMfac_MatF,fac,rng);      // natural mortality multiplier for mature new and old shell female
    pMfac_Big     = wts::jitterParameter(pMfac_Big,fac,rng);     // mult. on 1980-1984 M for mature males and females                     
    pRecAlpha  = wts::jitterParameter(pRecAlpha,fac,rng);  // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    pRecBeta    = wts::jitterParameter(pRecBeta,fac,rng);    // Parameters related to fraction recruiting  //this is NOT estimated (why?)
    
    pMnLnRec = wts::jitterParameter(pMnLnRec,fac,rng);         // Mean log-scale recruitment mnYrRecCurr+ (males, females are equal)
    pRecDevs = wts::jitterParameter(pRecDevs,0.1*fac,rng);     // Deviations about mean recruitment mnYrRecCurr+ (IMPORTANT CHANGE: used to be "endyr-1")
    pMnLnRecInit = wts::jitterParameter(pMnLnRecInit,fac,rng);    // Mean log-scale recruitment in early phase (pre-mnYrRecCurr)
    pRecDevsHist = wts::jitterParameter(pRecDevsHist,0.1*fac,rng);// Deviations about logscale mean recruitment in early phase (pre-mnYrRecCurr)
    
    pAvgLnF_TCF = wts::jitterParameter(pAvgLnF_TCF,fac,rng);           //log-scale mean directed fishing mortality
    pF_DevsTCF  = wts::jitterParameter(pF_DevsTCF,0.1*fac,rng);//log-scale directed fishing mortality devs IMPORTANT CHANGE: USED TO BE "1966,endyr-12"
    pAvgLnF_GTF = wts::jitterParameter(pAvgLnF_GTF,fac,rng);   // fishing mortality (trawl)
    pF_DevsGTF  = wts::jitterParameter(pF_DevsGTF,0.1*fac,rng);// trawl fishery f-devs       (IMPORTANT CHANGE: used to be "endyr") 1973 seems OK
    pAvgLnF_SCF = wts::jitterParameter(pAvgLnF_SCF,fac,rng);   // fishing mortality snow crab fishery discards
    pF_DevsSCF  = wts::jitterParameter(pF_DevsSCF,0.1*fac,rng);// snow crab fishery f-devs   (IMPORTANT CHANGE: used to be "endyr")  1992 is OK
    pAvgLnF_RKF = wts::jitterParameter(pAvgLnF_RKF,fac,rng);   // fishing mortality red king crab fishery discards //this is NOT estimated (why?)
    pF_DevsRKF  = wts::jitterParameter(pF_DevsRKF,0.1*fac,rng);//this is NOT estimated (why?)  IMPORTANT CHANGEA: was nObsDscRKF-1.  why -1 in "nObsDscRKF-1"
    
    // Retention function
    // 1981 - 1992
    pRetTCFM_slpA1 = wts::jitterParameter(pRetTCFM_slpA1,fac,rng);
    pRetTCFM_z50A1 = wts::jitterParameter(pRetTCFM_z50A1,fac,rng);
    // 2005-endyr  
    pRetTCFM_slpA2 = wts::jitterParameter(pRetTCFM_slpA2,fac,rng);
    pRetTCFM_z50A2 = wts::jitterParameter(pRetTCFM_z50A2,fac,rng);
    
    // Directed fishery selectivity pattern for period-1: 1993-1996
    pSelTCFM_slpA1 = wts::jitterParameter(pSelTCFM_slpA1,fac,rng);      
    
    // Directed fishery selectivity pattern changing by year for period-3: 2005-P
    pSelTCFM_slpA2 = wts::jitterParameter(pSelTCFM_slpA2,fac,rng);      
    pSelTCFM_mnLnZ50A2 = wts::jitterParameter(pSelTCFM_mnLnZ50A2,fac,rng);
    pSelTCFM_devsZ50 = wts::jitterParameter(pSelTCFM_devsZ50,0.1*fac,rng);
    
    // Female discards
    pSelTCFF_slp = wts::jitterParameter(pSelTCFF_slp,fac,rng);
    pSelTCFF_z50 = wts::jitterParameter(pSelTCFF_z50,fac,rng);
    
    // snow fishery female discards for period-1: 1989-1996
    pSelSCFF_slpA1 = wts::jitterParameter(pSelSCFF_slpA1,fac,rng);
    pSelSCFF_z50A1 = wts::jitterParameter(pSelSCFF_z50A1,fac,rng);
    
    // snow fishery female discards for period-2: 1997-2004
    pSelSCFF_slpA2 = wts::jitterParameter(pSelSCFF_slpA2,fac,rng);
    pSelSCFF_z50A2 = wts::jitterParameter(pSelSCFF_z50A2,fac,rng);
    
    // snow fishery female discards for period-3: 2005-P
    pSelSCFF_slpA3 = wts::jitterParameter(pSelSCFF_slpA3,fac,rng);
    pSelSCFF_z50A3 = wts::jitterParameter(pSelSCFF_z50A3,fac,rng);
    
    // snow fishery male discards for period-1: 1989-1996
    pSelSCFM_slpA1   = wts::jitterParameter(pSelSCFM_slpA1,fac,rng);
    pSelSCFM_z50A1   = wts::jitterParameter(pSelSCFM_z50A1,fac,rng);
    pSelSCFM_slpD1   = wts::jitterParameter(pSelSCFM_slpD1,fac,rng);
    pSelSCFM_lnZ50D1 = wts::jitterParameter(pSelSCFM_lnZ50D1,fac,rng);
    
    // snow fishery male discards for period-2: 1997-2004
    pSelSCFM_slpA2   = wts::jitterParameter(pSelSCFM_slpA2,fac,rng);
    pSelSCFM_z50A2   = wts::jitterParameter(pSelSCFM_z50A2,fac,rng);
    pSelSCFM_slpD2   = wts::jitterParameter(pSelSCFM_slpD2,fac,rng);
    pSelSCFM_lnZ50D2 = wts::jitterParameter(pSelSCFM_lnZ50D2,fac,rng);
    
    // snow fishery male discards for period-3: 2005-P
    pSelSCFM_slpA3   = wts::jitterParameter(pSelSCFM_slpA3,fac,rng);
    pSelSCFM_z50A3   = wts::jitterParameter(pSelSCFM_z50A3,fac,rng);
    pSelSCFM_slpD3   = wts::jitterParameter(pSelSCFM_slpD3,fac,rng);
    pSelSCFM_lnZ50D3 = wts::jitterParameter(pSelSCFM_lnZ50D3,fac,rng);
    
    // red king fishery female discards

    pSelRKFF_slpA1 = wts::jitterParameter(pSelRKFF_slpA1,fac,rng);
    pSelRKFF_z50A1 = wts::jitterParameter(pSelRKFF_z50A1,fac,rng);
    pSelRKFF_slpA2 = wts::jitterParameter(pSelRKFF_slpA2,fac,rng);
    pSelRKFF_z50A2 = wts::jitterParameter(pSelRKFF_z50A2,fac,rng);
    pSelRKFF_slpA3 = wts::jitterParameter(pSelRKFF_slpA3,fac,rng);
    pSelRKFF_z50A3 = wts::jitterParameter(pSelRKFF_z50A3,fac,rng);
    
    // red king fishery male discards
    pSelRKFM_slpA1 = wts::jitterParameter(pSelRKFM_slpA1,fac,rng);
    pSelRKFM_z50A1 = wts::jitterParameter(pSelRKFM_z50A1,fac,rng);
    pSelRKFM_slpA2 = wts::jitterParameter(pSelRKFM_slpA2,fac,rng);
    pSelRKFM_z50A2 = wts::jitterParameter(pSelRKFM_z50A2,fac,rng);
    pSelRKFM_slpA3 = wts::jitterParameter(pSelRKFM_slpA3,fac,rng);
    pSelRKFM_z50A3 = wts::jitterParameter(pSelRKFM_z50A3,fac,rng);
    
    // Trawl fishery selectivity female, 1973-1987
    pSelGTFF_slpA1 = wts::jitterParameter(pSelGTFF_slpA1,fac,rng);
    pSelGTFF_z50A1 = wts::jitterParameter(pSelGTFF_z50A1,fac,rng);
    // Trawl fishery selectivity female, 1988-1996
    pSelGTFF_slpA2 = wts::jitterParameter(pSelGTFF_slpA2,fac,rng);
    pSelGTFF_z50A2 = wts::jitterParameter(pSelGTFF_z50A2,fac,rng);
    // Trawl fishery selectivity female, 1997-P
    pSelGTFF_slpA3 = wts::jitterParameter(pSelGTFF_slpA3,fac,rng);
    pSelGTFF_z50A3 = wts::jitterParameter(pSelGTFF_z50A3,fac,rng);
    // Trawl fishery selectivity male, 1973-1987
    pSelGTFM_slpA1 = wts::jitterParameter(pSelGTFM_slpA1,fac,rng);
    pSelGTFM_z50A1 = wts::jitterParameter(pSelGTFM_z50A1,fac,rng);
    // Trawl fishery selectivity male, 1988-1996
    pSelGTFM_slpA2 = wts::jitterParameter(pSelGTFM_slpA2,fac,rng);
    pSelGTFM_z50A2 = wts::jitterParameter(pSelGTFM_z50A2,fac,rng);
    // Trawl fishery selectivity male, 1997-P
    pSelGTFM_slpA3 = wts::jitterParameter(pSelGTFM_slpA3,fac,rng);
    pSelGTFM_z50A3 = wts::jitterParameter(pSelGTFM_z50A3,fac,rng);
    //1974 to 1981 
    pSrv1_QM       = wts::jitterParameter(pSrv1_QM,fac,rng);
    pSrv1M_dz5095 = wts::jitterParameter(pSrv1M_dz5095,fac,rng);
    pSrv1M_z50   = wts::jitterParameter(pSrv1M_z50,fac,rng);
    //1982-P
    pSrv2_QM       = wts::jitterParameter(pSrv2_QM,fac,rng);
    pSrv2M_dz5095 = wts::jitterParameter(pSrv2M_dz5095,fac,rng);
    pSrv2M_z50   = wts::jitterParameter(pSrv2M_z50,fac,rng);
    
//    pPrM2MF = wts::jitterParameter(pPrM2MF,fac,rng);
//    pPrM2MM = wts::jitterParameter(pPrM2MM,fac,rng);
    
    pSrv1_QF      = wts::jitterParameter(pSrv1_QF,fac,rng);
    pSrv1F_dz5095 = wts::jitterParameter(pSrv1F_dz5095,fac,rng);    
    pSrv1F_z50   = wts::jitterParameter(pSrv1F_z50,fac,rng);
    
    pSrv2_QF      = wts::jitterParameter(pSrv2_QF,fac,rng);
    pSrv2F_dz5095 = wts::jitterParameter(pSrv2F_dz5095,fac,rng);
    pSrv2F_z50   = wts::jitterParameter(pSrv2F_z50,fac,rng);
    
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
FUNCTION closeMCMCFile
    mcmc<<"dummy=0)"<<endl;
    mcmc.close();
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION WriteMCMC
    post<<
        pSelTCFF_slp <<","<<
        pSelTCFF_z50 <<","<<
    endl;

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
FUNCTION get_maturity
    //prMoltToMaturity(female|size)
    if(active(pPrM2MF)){
        //prM2M(females> zBs(16)) assumed = 1
        if (optPrM2M==0){
            modPrM2M(FEMALE)(1,16) = mfexp(pPrM2MF);
        } else {
            modPrM2M(FEMALE)(1,16) = 1.0/(1.0+mfexp(-pPrM2MF));
        }
    }

    //prMoltToMaturity(female|size)
    if(active(pPrM2MM)){
        if (optPrM2M==0){
            modPrM2M(MALE) = mfexp(pPrM2MM);
        } else {
            modPrM2M(MALE) = 1.0/(1.0+mfexp(-pPrM2MM));
        }
    }
//     CheckFile<<"modPrM2M"<<endl;
//     CheckFile<<modPrM2M(1)<<endl;
//     CheckFile<<modPrM2M(2)<<endl;
    
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// NOTE: function renamed from get_growth() and growth_switch flag check removed
FUNCTION get_growth1                                   //wts: revised
    
    meanPostMoltSize(FEMALE)= mfexp(pGrAF1)* pow(zBs,pGrBF1);
    meanPostMoltSize(MALE)  = mfexp(pGrAM1)* pow(zBs,pGrBM1);
    
    //  CheckFile << af<<" " <<bf<<" "<<am<<" "<<bm<<endl;
//     CheckFile<<"meanPostMoltSize"<<endl;
//     CheckFile<<meanPostMoltSize(1)<<endl;
//     CheckFile<<meanPostMoltSize(2)<<endl;
    
    // using Gamma function for transition matrix
    // devia is the bounds of growth bins to evaluate
    // the gamma function (x) in prop = integral(i1 to i2) g(x|alpha,beta) dx
    // alpha and pGrBeta_x are parameters 
    // alpha is the mean growth increment per molt for some premolt length class
    // alpha = mean growth increment per molt divided by beta
    // beta is the shape parameter - larger beta - more variance 

    double devia;
    dvariable alpha;    
    prGr_xzz.initialize();
    for (int sex=FEMALE;sex<=MALE;sex++) {
        for (int ilen=1;ilen<=nZBs;ilen++){    
            // subract the 2.5 from the midpoint of the length bin to get the lower bound
            alpha = (meanPostMoltSize(sex,ilen)-(zBs(ilen)-2.5))/pGrBeta_x(sex);
            //    cout<<"alpha = "<<alpha<<endl;
            //    cout<<"pGrBeta_x = "<<pGrBeta_x<<endl;
            //    cout<<"meanPostMoltSize = "<<meanPostMoltSize(sex,ilen)<<endl;
             //    truncate growth transition to max=10 bins
            for (int il2=ilen;il2<=ilen+min(10,nZBs-ilen);il2++) {
                devia = zBs(il2)+2.5-zBs(ilen);
                prGr_xzz(sex,ilen,il2) = pow(devia,(alpha-1.))*exp(-devia/pGrBeta_x(sex));
            }  
            //standardize so each row sums to 1.0
            prGr_xzz(sex,ilen) /= sum(prGr_xzz(sex,ilen));
        }
    }
//     CheckFile<<"prGr_xzz"<<endl;
//     CheckFile<<prGr_xzz(1)<<endl;
//     CheckFile<<prGr_xzz(2)<<endl;
    
    // Fraction recruiting
    dvariable alpha_rec = pRecAlpha/pRecBeta;
    for (int ilen=1;ilen<=nZBs;ilen++) {
        devia = zBs(ilen)+2.5-zBs(1);     //wts: think adding 2.5 here is wrong! just want offset from 1st bin (unless these are cutputs)
        prRec_z(ilen) = pow(devia,alpha_rec-1.)*exp(-devia/pRecBeta);//mfexp(...) doesn't work here?
    }
    prRec_z /= sum(prRec_z); //standardize so each row sums to 1.0
    
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
FUNCTION get_moltingp                         //wts: revised
    //assume all immature crab molt
    prMoltImm_xz(FEMALE)=1.0;
    prMoltImm_xz(MALE)  =1.0;
    
    //assume mature crab do not molt
    prMoltMat_xz(FEMALE)=0.0;
    prMoltMat_xz(  MALE)=0.0;

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
FUNCTION get_selectivity                  //wts: revised
//    cout<<"get_selectivity"<<endl;
    int ii = 1;
    
    selTCFM_syz.initialize();
    retFcn_syz.initialize();
    dvariable tmpSel50 = mean(mfexp(pSelTCFM_mnLnZ50A2+pSelTCFM_devsZ50(1,6)));
    for(int iy=styr;iy<=1990;iy++){ 
        selTCFM_syz(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*pSelTCFM_slpA1*(zBs-tmpSel50)));    
        retFcn_syz(NEW_SHELL, iy) = 1./(1.+mfexp(-1.*pRetTCFM_slpA1*(zBs-pRetTCFM_z50A1)));
    }
//    cout<<"get_sel: 1a"<<endl;
    int ctr = 1;
//    cout<<"max index of pSelTCFM_devsZ50: "<<pSelTCFM_devsZ50.indexmax()<<endl;
    for(int iy=1991;iy<=1996;iy++){ 
        if (hasDirectedFishery_y(iy)) {
//            cout<<"yr = "<<iy<<".  ctr = "<<ctr<<endl;
            selTCFM_syz(NEW_SHELL,iy)=1./(1.+mfexp(-1.*pSelTCFM_slpA1*(zBs-mfexp(pSelTCFM_mnLnZ50A2+pSelTCFM_devsZ50(ctr++)))));//ctr was iy-1990
        } else {
//            cout<<"yr = "<<iy<<".  no fishery."<<endl;
        }
    }
//    cout<<"get_sel: 1b"<<endl;
    //no directed fishery, set 50% selectivity to mean
    for(int iy=1997;iy<endyr;iy++){ 
        if (hasDirectedFishery_y(iy)) {
//            cout<<"yr = "<<iy<<".  ctr = "<<ctr<<endl;
            selTCFM_syz(NEW_SHELL,iy)=1./(1.+mfexp(-1.*pSelTCFM_slpA2*(zBs-mfexp(pSelTCFM_mnLnZ50A2+pSelTCFM_devsZ50(ctr++)))));//ctr was iy-1998
        } else {
//            cout<<"yr = "<<iy<<".  no fishery."<<endl;
            selTCFM_syz(NEW_SHELL,iy)=1./(1.+mfexp(-1.*pSelTCFM_slpA2*(zBs-mfexp(pSelTCFM_mnLnZ50A2))));
        }
    }
//    cout<<"get_sel: 1d"<<endl;
    
    for(int iy=styr;iy<=1990;iy++) retFcn_syz(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*pRetTCFM_slpA1*(zBs-pRetTCFM_z50A1)));
//    cout<<"get_sel: 1f"<<endl;
    for(int iy=1991;iy<endyr;iy++) retFcn_syz(NEW_SHELL,iy) = 1./(1.+mfexp(-1.*pRetTCFM_slpA2*(zBs-pRetTCFM_z50A2)));
//    cout<<"get_sel: 1g"<<endl;
    
    for(int iy=styr;iy<endyr;iy++){      //used to be iy<=endyr            
        // set new and old selTCFM_syz same
        selTCFM_syz(OLD_SHELL,iy) = selTCFM_syz(NEW_SHELL,iy);
        retFcn_syz(OLD_SHELL,iy)  = retFcn_syz(NEW_SHELL,iy);
    }//year loop
//    cout<<"get_sel: 2"<<endl;
    
    // female discards ascending logistic curve 
    selTCFF_z=1./(1.+mfexp(-1.*pSelTCFF_slp*(zBs-pSelTCFF_z50)));
    
    //  snow fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selSCF_cxz(1,FEMALE)=1./(1.+mfexp(-1.*pSelSCFF_slpA1*(zBs-pSelSCFF_z50A1))); 
    selSCF_cxz(2,FEMALE)=1./(1.+mfexp(-1.*pSelSCFF_slpA2*(zBs-pSelSCFF_z50A2))); 
    selSCF_cxz(3,FEMALE)=1./(1.+mfexp(-1.*pSelSCFF_slpA3*(zBs-pSelSCFF_z50A3))); 
//    cout<<"get_sel: 2a"<<endl;
        
    //  snow fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selSCF_cxz(1,MALE)=elem_prod(1./(1.+mfexp(-1.*pSelSCFM_slpA1*(zBs-pSelSCFM_z50A1))),
                                 1./(1.+mfexp(pSelSCFM_slpD1*(zBs-(pSelSCFM_z50A1+mfexp(pSelSCFM_lnZ50D1))))));
    selSCF_cxz(2,MALE)=elem_prod(1./(1.+mfexp(-1.*pSelSCFM_slpA2*(zBs-pSelSCFM_z50A2))),
                                 1./(1.+mfexp(pSelSCFM_slpD2*(zBs-(pSelSCFM_z50A2+mfexp(pSelSCFM_lnZ50D2))))));
    selSCF_cxz(3,MALE)=elem_prod(1./(1.+mfexp(-1.*pSelSCFM_slpA3*(zBs-pSelSCFM_z50A3))),
                                 1./(1.+mfexp(pSelSCFM_slpD3*(zBs-(pSelSCFM_z50A3+mfexp(pSelSCFM_lnZ50D3))))));
//    cout<<"get_sel: 2b"<<endl;
    
    //  red fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selRKF_cxz(1,FEMALE)=1./(1.+mfexp(-1.*pSelRKFF_slpA1*(zBs-pSelRKFF_z50A1))); 
    selRKF_cxz(2,FEMALE)=1./(1.+mfexp(-1.*pSelRKFF_slpA2*(zBs-pSelRKFF_z50A2))); 
    selRKF_cxz(3,FEMALE)=1./(1.+mfexp(-1.*pSelRKFF_slpA3*(zBs-pSelRKFF_z50A3))); 
    
    //  red fishery selectivity for 3 time periods, #1 (1989-1996), #2 (1997-2004) and #3 (2005-P)      
    selRKF_cxz(1,MALE)=1./(1.+mfexp(-1.*pSelRKFM_slpA1*(zBs-pSelRKFM_z50A1))); 
    selRKF_cxz(2,MALE)=1./(1.+mfexp(-1.*pSelRKFM_slpA2*(zBs-pSelRKFM_z50A2))); 
    selRKF_cxz(3,MALE)=1./(1.+mfexp(-1.*pSelRKFM_slpA3*(zBs-pSelRKFM_z50A3))); 
//    cout<<"get_sel: 2c"<<endl;
    
    //  trawl fishery selectivity for 3 time periods, #1 (1973-1987), #2 (1988-1996) and #3 (1997-P)
    selGTF_cxz(1,FEMALE)=1./(1.+mfexp(-1.*pSelGTFF_slpA1*(zBs-pSelGTFF_z50A1)));
    selGTF_cxz(2,FEMALE)=1./(1.+mfexp(-1.*pSelGTFF_slpA2*(zBs-pSelGTFF_z50A2)));
    selGTF_cxz(3,FEMALE)=1./(1.+mfexp(-1.*pSelGTFF_slpA3*(zBs-pSelGTFF_z50A3)));
    
    selGTF_cxz(1,MALE)=1./(1.+mfexp(-1.*pSelGTFM_slpA1*(zBs-pSelGTFM_z50A1)));    
    selGTF_cxz(2,MALE)=1./(1.+mfexp(-1.*pSelGTFM_slpA2*(zBs-pSelGTFM_z50A2)));
    selGTF_cxz(3,MALE)=1./(1.+mfexp(-1.*pSelGTFM_slpA3*(zBs-pSelGTFM_z50A3)));
//    cout<<"get_sel: 2d"<<endl;
        
    dvariable maxsel;
    selTCFR_syz.initialize();
    for(int iy=styr;iy<endyr;iy++){          //used to be iy<=endyr
        maxsel = max(selTCFM_syz(NEW_SHELL,iy));
        if(maxsel<max(selTCFM_syz(OLD_SHELL,iy))) maxsel = max(selTCFM_syz(OLD_SHELL,iy));  //wts: is this differentiable??
        for (int shell=NEW_SHELL;shell<=OLD_SHELL;shell++){
            selTCFM_syz(shell,iy) = selTCFM_syz(shell,iy)/maxsel;
            selTCFR_syz(shell,iy) = elem_prod(retFcn_syz(shell,iy),selTCFM_syz(shell,iy));
        }
    }
//    cout<<"get_sel: 4"<<endl;
//    cout<<"done"<<endl;
    
     if (optFshSel==1){//set logistic selectivity = 1 in largest size bin
        //TCFM and retFcn_syz
        for(int iy=styr;iy<endyr;iy++){
            selTCFM_syz(NEW_SHELL,iy) /= selTCFM_syz(NEW_SHELL,iy,nZBs);
            retFcn_syz( NEW_SHELL,iy) /= retFcn_syz( NEW_SHELL,iy,nZBs);
            selTCFM_syz(OLD_SHELL,iy) /= selTCFM_syz(OLD_SHELL,iy,nZBs);
            retFcn_syz( OLD_SHELL,iy) /= retFcn_syz( OLD_SHELL,iy,nZBs);
        }
        //TCFF
        selTCFF_z /= selTCFF_z(nZBs);
        //SCF females (only)
        selSCF_cxz(1,FEMALE) /= selSCF_cxz(1,FEMALE,nZBs);
        selSCF_cxz(2,FEMALE) /= selSCF_cxz(2,FEMALE,nZBs);
        selSCF_cxz(3,FEMALE) /= selSCF_cxz(3,FEMALE,nZBs);
        //RKF
        selRKF_cxz(1,  MALE) /= selRKF_cxz(1,  MALE,nZBs);
        selRKF_cxz(2,  MALE) /= selRKF_cxz(2,  MALE,nZBs);
        selRKF_cxz(3,  MALE) /= selRKF_cxz(3,  MALE,nZBs);
        selRKF_cxz(1,FEMALE) /= selRKF_cxz(1,FEMALE,nZBs);
        selRKF_cxz(2,FEMALE) /= selRKF_cxz(2,FEMALE,nZBs);
        selRKF_cxz(3,FEMALE) /= selRKF_cxz(3,FEMALE,nZBs);
        //GTF
        selGTF_cxz(1,  MALE) /= selGTF_cxz(1,  MALE,nZBs);
        selGTF_cxz(2,  MALE) /= selGTF_cxz(2,  MALE,nZBs);
        selGTF_cxz(3,  MALE) /= selGTF_cxz(3,  MALE,nZBs);
        selGTF_cxz(1,FEMALE) /= selGTF_cxz(1,FEMALE,nZBs);
        selGTF_cxz(2,FEMALE) /= selGTF_cxz(2,FEMALE,nZBs);
        selGTF_cxz(3,FEMALE) /= selGTF_cxz(3,FEMALE,nZBs);
     }
    
    //calculate survey selectivities
    selSrv1_xz.initialize();
    selSrv2_xz.initialize();
    selSrv3_xz.initialize();
    selSrv1_xz(MALE) = 1./(1.+mfexp(-1.*log(19.)*(zBs-pSrv1M_z50)/(pSrv1M_dz5095)));
    selSrv2_xz(MALE) = 1./(1.+mfexp(-1.*log(19.)*(zBs-pSrv2M_z50)/(pSrv2M_dz5095)));
    selSrv3_xz(MALE) = 1./(1.+mfexp(-1.*log(19.)*(zBs-pSrv2M_z50)/(pSrv2M_dz5095)));
        
    selSrv1_xz(FEMALE) = 1./(1.+mfexp(-1.*log(19.)*(zBs-pSrv1F_z50)/(pSrv1F_dz5095)));
    selSrv2_xz(FEMALE) = 1./(1.+mfexp(-1.*log(19.)*(zBs-pSrv2F_z50)/(pSrv2F_dz5095)));
    selSrv3_xz(FEMALE) = 1./(1.+mfexp(-1.*log(19.)*(zBs-pSrv2F_z50)/(pSrv2F_dz5095)));
//    cout<<"get_sel: 3"<<endl;
    
    //survey selectivities
     if (optSrvSel==1){//set logistic selectivity = 1 in largest size bin
        selSrv1_xz(  MALE) /= selSrv1_xz(  MALE,nZBs);
        selSrv2_xz(  MALE) /= selSrv2_xz(  MALE,nZBs);
        selSrv3_xz(  MALE) /= selSrv3_xz(  MALE,nZBs);
        selSrv1_xz(FEMALE) /= selSrv1_xz(FEMALE,nZBs);
        selSrv2_xz(FEMALE) /= selSrv2_xz(FEMALE,nZBs);
        selSrv3_xz(FEMALE) /= selSrv3_xz(FEMALE,nZBs);
    }
    
    selSrv1_xz(MALE) *= pSrv1_QM;
    selSrv2_xz(MALE) *= pSrv2_QM;
    selSrv3_xz(MALE) *= pSrv2_QM;
    
    // use somerton and otto curve for survey selectivities, as required
    if (useSomertonOtto1) selSrv1_xz(MALE) = selSO_z;
    if (useSomertonOtto2) selSrv2_xz(MALE) = selSO_z;
    if (useSomertonOtto3) selSrv3_xz(MALE) = selSO_z;
        
    //scale survey selectivity by survey catchability
    selSrv1_xz(FEMALE)  *= pSrv1_QF;
    selSrv2_xz(FEMALE) *= pSrv2_QF;
    selSrv3_xz(FEMALE)  *= pSrv2_QF;
//    cout<<"get_sel: 3"<<endl;
    //<--
    
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
FUNCTION get_mortality
    int debug = 0;
    if (debug) cout<<"get_mortality"<<endl;
    int ii;
    int inc;
    
    M_msx(IMMATURE,NEW_SHELL,FEMALE) = baseM_msx(IMMATURE,NEW_SHELL,FEMALE)*pMfac_Imm;
    M_msx(IMMATURE,NEW_SHELL,  MALE) = baseM_msx(IMMATURE,NEW_SHELL,  MALE)*pMfac_Imm;
    M_msx(IMMATURE,OLD_SHELL,FEMALE) = baseM_msx(IMMATURE,OLD_SHELL,FEMALE)*pMfac_Imm;
    M_msx(IMMATURE,OLD_SHELL,  MALE) = baseM_msx(IMMATURE,OLD_SHELL,  MALE)*pMfac_Imm;
    M_msx(  MATURE,NEW_SHELL,FEMALE) = baseM_msx(  MATURE,NEW_SHELL,FEMALE)*pMfac_MatF;
    M_msx(  MATURE,NEW_SHELL,  MALE) = baseM_msx(  MATURE,NEW_SHELL,  MALE)*pMfac_MatM;
    M_msx(  MATURE,OLD_SHELL,FEMALE) = baseM_msx(  MATURE,OLD_SHELL,FEMALE)*pMfac_MatF;
    M_msx(  MATURE,OLD_SHELL,  MALE) = baseM_msx(  MATURE,OLD_SHELL,  MALE)*pMfac_MatM;
    if (debug) cout<<"0"<<endl;
    
    //first year retained catch 1965(1966 fishery) no fishery 1985, 1986 or 1997-2004 or 2010-2012
    fTCF_xy.initialize();
    fTCF_xy(MALE)(styr,1964) = 0.05;//was 1965!!
//    cout<<"0a"<<endl;
    int idx = 1;
    for(int iy =1965;iy<endyr;iy++){
        if(hasDirectedFishery_y(iy)) fTCF_xy(MALE,iy) = mfexp(pAvgLnF_TCF+pF_DevsTCF(idx++));
    }
    fTCF_xy(FEMALE) = fTCF_xy(MALE)*mfexp(pAvgLnF_TCFF);
    if (debug) cout<<"1"<<endl;
    
    // fmortdf=mfexp(log_avg_fmortdf+fmortdf_dev); using overall fmTCFM_syz for females as well as males in directed fishery
    //Fs in snow and BBRKC fishery are scalars need to multiply in projections by retained snow crab/average snow catch * fmTCFM_syz to get fmTCFM_syz.
    //20150601: ratio is now either mortality rate/effort OR fishing capture rate/effort
    //20160325: scaling factor can now be related to a model parameter
    if (active(pLnEffXtr_SCF)){
        qSCF = mfexp(pLnEffXtr_SCF);//parameter being estimated
    } else {
        dvar_vector f_SCF1 = mfexp(pAvgLnF_SCF+pF_DevsSCF);
        qSCF = mean(f_SCF1);
    }
    if (debug) cout<<" qSCF = "<<qSCF<<endl; 
    if (debug) cout<<"2"<<endl;
    
    fSCF_xy.initialize();
    fSCF_xy(MALE)(styr,1977)= 0.01;
//    for(int iy=1978;iy<=1991;iy++) fSCF_xy(MALE)(iy) = qSCF*effSCF_y(iy);
    fSCF_xy(MALE)(1978,endyr-1) = qSCF*effSCF_y(1978,endyr-1)/mnEff_SCF;
    fSCF_xy(MALE)(1992,endyr-1) = mfexp(pAvgLnF_SCF+pF_DevsSCF);
    fSCF_xy(FEMALE) = fSCF_xy(MALE)*mfexp(pAvgLnF_SCFF);
    if (debug) cout<<"3"<<endl;
    
    // need to have the devs 1992 to present
    //20150601: ratio is now either mortality rate/effort OR fishing capture rate/effort
//    cout<<"4a"<<endl;
    fRKF_xy.initialize();
    fRKF_xy(MALE)(styr,1952)= 0.02; //qRKF*mean(rkccatch(1969,1973));
//    cout<<"4a"<<endl;
    //20160325: scaling factor can now be related to a model parameter
    if (active(pLnEffXtr_RKF)){
        qRKF = mfexp(pLnEffXtr_RKF);//parameter being estimated
    } else {
        dvar_vector f_RKF1(1,nObsDscRKF);   //was nObsDscRKF-1
        f_RKF1 = mfexp(pAvgLnF_RKF+pF_DevsRKF);
        if (optEffXtr_RKF==1){
            qRKF = mean(f_RKF1);
        } else if (optEffXtr_RKF==2) {
            qRKF = mean(1-mfexp(-f_RKF1));//old style
        }
    }
    if (optEffXtr_RKF==1){
        fRKF_xy(MALE)(1953,endyr-1) = qRKF*effRKF_y(1953,endyr-1)/mnEff_RKF;
    } else if (optEffXtr_RKF==2) {
        fRKF_xy(MALE)(1953,endyr-1) = -log(1-qRKF*effRKF_y(1953,endyr-1)/mnEff_RKF);//old style
    }
    if (optMinFs>0) {for (int iy=1953;iy<=1991;iy++) if(fRKF_xy(MALE)(iy)< 0.01) fRKF_xy(MALE)(iy) = 0.01;}
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
            tmpSelSCF(FEMALE)=selSCF_cxz(1,FEMALE);
            tmpSelSCF(MALE)  =selSCF_cxz(1,MALE);
        }
        if (1997<=iy && iy<=2004) {
            tmpSelSCF(FEMALE)=selSCF_cxz(2,FEMALE);
            tmpSelSCF(MALE)  =selSCF_cxz(2,MALE);
        }
        if (2005<=iy) {
            tmpSelSCF(FEMALE)=selSCF_cxz(3,FEMALE);
            tmpSelSCF(MALE)  =selSCF_cxz(3,MALE);
        }
        
        // test on year for 3 red selectivity periods
        if (iy<=1996) {
            tmpSelRKF(FEMALE)=selRKF_cxz(1,FEMALE);
            tmpSelRKF(MALE)  =selRKF_cxz(1,MALE);
        }
        if (1997<=iy && iy<=2004) {
            tmpSelRKF(FEMALE)=selRKF_cxz(2,FEMALE);
            tmpSelRKF(MALE)  =selRKF_cxz(2,MALE);
        }
        if (2005<=iy) {
            tmpSelRKF(FEMALE)=selRKF_cxz(3,FEMALE);
            tmpSelRKF(MALE)  =selRKF_cxz(3,MALE);
        }
        
        // test on year for 3 trawl selectivity periods
        if (iy<=1986) {
            tmpSelGTF(FEMALE)=selGTF_cxz(1,FEMALE);
            tmpSelGTF(MALE)  =selGTF_cxz(1,MALE);
        }
        if (1987<=iy && iy<=1996) {
            tmpSelGTF(FEMALE)=selGTF_cxz(2,FEMALE);
            tmpSelGTF(MALE)  =selGTF_cxz(2,MALE);
        }
        if (1997<=iy) {
            tmpSelGTF(FEMALE)=selGTF_cxz(3,FEMALE);
            tmpSelGTF(MALE)  =selGTF_cxz(3,MALE);
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
            fmTCFF_yz(iy)= selTCFF_z*fTCF_xy(FEMALE,iy);
            for(int s=NEW_SHELL;s<=OLD_SHELL;s++) {
                fmTCFM_syz(s,iy)     = selTCFM_syz(s,iy)*fTCF_xy(MALE,iy);//total fishing mortality on males in directed fishery       
                fmTCFR_syz(s,iy)     = selTCFR_syz(s,iy)*fTCF_xy(MALE,iy);//retained fishing mortality on males in directed fishery
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
            fcTCFF_yz(iy)= selTCFF_z*fTCF_xy(FEMALE,iy);
            for(int s=NEW_SHELL;s<=OLD_SHELL;s++) fcTCFM_syz(s,iy) = selTCFM_syz(s,iy)*fTCF_xy(MALE,iy);       
        
            //mortality rates
            fmSCF_xyz(FEMALE,iy) = hm_pot*fcSCF_xyz(FEMALE,iy);   
            fmSCF_xyz(  MALE,iy) = hm_pot*fcSCF_xyz(  MALE,iy);   
            fmRKF_xyz(FEMALE,iy) = hm_pot*fcRKF_xyz(FEMALE,iy);   
            fmRKF_xyz(  MALE,iy) = hm_pot*fcRKF_xyz(  MALE,iy);  
            fmGTF_xyz(FEMALE,iy) = hm_trawl*fcGTF_xyz(FEMALE,iy);  
            fmGTF_xyz(  MALE,iy) = hm_trawl*fcGTF_xyz(  MALE,iy);   
            fmTCFF_yz(iy)        = hm_pot*fcTCFF_yz(iy);
            for (int s=NEW_SHELL;s<=OLD_SHELL;s++) { //over new (shell=1) and old (shell=2) shell...
                fmTCFR_syz(s,iy) = elem_prod(retFcn_syz(s,iy),             fcTCFM_syz(s,iy)); //retention rate
                fmTCFD_syz(s,iy) = elem_prod(hm_pot*(1.0-retFcn_syz(s,iy)),fcTCFM_syz(s,iy));  //discard mortality rate
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
    modNum_xyz.initialize();
    modNum_yxmsz.initialize();
    natlength_new.initialize();
    natlength_old.initialize();
    natlength_imm.initialize();
    natlength_mat.initialize();
    
    //nAtZ_msxy.initialize();
    rec_y(styr,mnYrRecCurr-1)            = mfexp(pMnLnRecInit);
    rec_y(mnYrRecDevsHist,mnYrRecCurr-1) = mfexp(pMnLnRecInit+pRecDevsHist);
    rec_y(mnYrRecCurr,endyr)             = mfexp(pMnLnRec+pRecDevs);
    
    for (int x=1;x<=nSXs;x++) {  
        //initialize modNum_yxmsz in styr with recruitment in new shell, immature
        modNum_yxmsz(styr,x,IMMATURE,NEW_SHELL) += 0.5*rec_y(styr)*prRec_z;  
        modNum_xyz(x,styr) = modNum_yxmsz(styr,x,IMMATURE,NEW_SHELL);
        for (int yr=styr;yr<endyr;yr++) {
            Surv1 = mfexp(-M_msx(IMMATURE,NEW_SHELL,x));
            Surv2 = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,x));
            if((lyr_mort<=yr) && (yr<=uyr_mort) && (mort_switch==1)) {
//                 cout<<yr<<": applying big_mort"<<endl;
                Surv3 = mfexp(-M_msx(MATURE,NEW_SHELL,x)*pMfac_Big(x));
                Surv4 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,x)*pMfac_Big(x));
                Surv5 = mfexp(-M_msx(MATURE,OLD_SHELL, x)*pMfac_Big(x));
                Surv6 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL, x)*pMfac_Big(x));
            } else {
                Surv3 = mfexp(-M_msx(MATURE,NEW_SHELL,x));
                Surv4 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,x));
                Surv5 = mfexp(-M_msx(MATURE,OLD_SHELL,x));
                Surv6 = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,x));
            }
            // Numbers advancing to new shell...
            dvar_vector tmpn = Surv1*elem_prod(prMoltImm_xz(x),elem_prod(S_xsyz(x,NEW_SHELL,yr),modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL)));
            dvar_vector tmpo = Surv1*elem_prod(prMoltImm_xz(x),elem_prod(S_xsyz(x,OLD_SHELL,yr),modNum_yxmsz(yr,x,IMMATURE,OLD_SHELL)));
            natlength_new(x,yr+1) =  (tmpn+tmpo) * prGr_xzz(x);
            modNum_yxmsz(yr+1,x,IMMATURE,OLD_SHELL) = Surv1*(elem_prod(S_xsyz(x,NEW_SHELL,yr),modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL))+
                                                             elem_prod(S_xsyz(x,OLD_SHELL,yr),modNum_yxmsz(yr,x,IMMATURE,OLD_SHELL))) - tmpn-tmpo;
            
            dvar_vector tmpmn = Surv3*elem_prod(prMoltMat_xz(x),elem_prod(S_xsyz(x,NEW_SHELL,yr),modNum_yxmsz(yr,x,  MATURE,NEW_SHELL)));
            dvar_vector tmpmo = Surv5*elem_prod(prMoltMat_xz(x),elem_prod(S_xsyz(x,OLD_SHELL,yr),modNum_yxmsz(yr,x,  MATURE,OLD_SHELL)));
            modNum_yxmsz(yr+1,x,  MATURE,NEW_SHELL) = (tmpmn+tmpmo) * prGr_xzz(x);
            modNum_yxmsz(yr+1,x,  MATURE,OLD_SHELL) = Surv3 * elem_prod(S_xsyz(x,NEW_SHELL,yr),modNum_yxmsz(yr,x,  MATURE,NEW_SHELL)) + 
                                                      Surv5 * elem_prod(S_xsyz(x,OLD_SHELL,yr),modNum_yxmsz(yr,x,  MATURE,OLD_SHELL)) - tmpmn-tmpmo;
            
            // this is for estimating the fraction of new shell that move to old shell to fit
            // the survey data that is split by immature and mature
            modNum_yxmsz(yr+1,x,  MATURE,NEW_SHELL) += elem_prod(    modPrM2M(x),natlength_new(x,yr+1));//add new shell that become mature
            modNum_yxmsz(yr+1,x,IMMATURE,NEW_SHELL)  = elem_prod(1.0-modPrM2M(x),natlength_new(x,yr+1));//new shell that stay immature
            //     cout<<" to 2 "<<endl;
            // add in recruits for next year
            // put all recruits in new shell immature
            modNum_yxmsz(yr+1,x,IMMATURE,NEW_SHELL) += 0.5*rec_y(yr+1)*prRec_z;
            natlength_new(x,yr+1)   = modNum_yxmsz(yr+1,x,IMMATURE,NEW_SHELL) + modNum_yxmsz(yr+1,x,  MATURE,NEW_SHELL);//no real reason to recalc this!
            natlength_old(x,yr+1)   = modNum_yxmsz(yr+1,x,  MATURE,OLD_SHELL) + modNum_yxmsz(yr+1,x,IMMATURE,OLD_SHELL);
            natlength_mat(x,yr+1)   = modNum_yxmsz(yr+1,x,  MATURE,NEW_SHELL) + modNum_yxmsz(yr+1,x,  MATURE,OLD_SHELL);
            natlength_imm(x,yr+1)   = modNum_yxmsz(yr+1,x,IMMATURE,NEW_SHELL) + modNum_yxmsz(yr+1,x,IMMATURE,OLD_SHELL);
            modNum_xyz(x,yr+1)      = natlength_mat(x,yr+1)  + natlength_imm(x,yr+1);
        }//yr
    }//x
//    cout<<"2"<<endl;
    
    //n-at-length just before fishing occurs
    for (int yr=styr;yr<endyr;yr++){
        for (int x=1;x<=nSXs;x++) {
            if (mort_switch2) {//correct way of applying this; wts: new 2013-08-28
                modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL);
                modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,OLD_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,OLD_SHELL);
                if(lyr_mort<=yr && yr<=uyr_mort && mort_switch==1) {
                    modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,x)*pMfac_Big(x))*modNum_yxmsz(yr,x,MATURE,NEW_SHELL);
                    modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,x)*pMfac_Big(x))*modNum_yxmsz(yr,x,MATURE,OLD_SHELL);
                } else {  
                    modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,MATURE,NEW_SHELL);
                    modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,x))*modNum_yxmsz(yr,x,MATURE,OLD_SHELL);
                }
                natl_new_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL); 
                natl_old_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL);
            } else { //2012 way of doing it               
                if(lyr_mort<=yr && yr<=uyr_mort && mort_switch==1) {
                    modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL);
                    modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,OLD_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,OLD_SHELL);
                    modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,x)*pMfac_Big(x))*modNum_yxmsz(yr,x,MATURE,NEW_SHELL);
                    modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,x)*pMfac_Big(x))*modNum_yxmsz(yr,x,MATURE,OLD_SHELL);
                    natl_new_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL); 
                    natl_old_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL);
                }
                modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL);
                modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(IMMATURE,OLD_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,OLD_SHELL);
                modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,  MATURE,NEW_SHELL);
                modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,x))*modNum_yxmsz(yr,x,  MATURE,OLD_SHELL);
                natl_new_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL); 
                natl_old_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL);
            }//mort_switch2
        }//x
    }//yr
//    cout<<"3"<<endl;
    //assume mdptFshs_y(endyr) = mdptFshs_y(endyr-1)
    {   int yr = endyr;    //don't worry about mort_switch here
        for (int x=1;x<=nSXs;x++) {
            modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr-1)*M_msx(IMMATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL);
            modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr-1)*M_msx(IMMATURE,OLD_SHELL,x))*modNum_yxmsz(yr,x,IMMATURE,OLD_SHELL);
            modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL) = mfexp(-mdptFshs_y(yr-1)*M_msx(MATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,  MATURE,NEW_SHELL);
            modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL) = mfexp(-mdptFshs_y(yr-1)*M_msx(MATURE,OLD_SHELL,x))*modNum_yxmsz(yr,x,  MATURE,OLD_SHELL);
            natl_new_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,NEW_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,NEW_SHELL); 
            natl_old_fishtime(x,yr)  = modFT_PopNum_yxmsz(yr,x,IMMATURE,OLD_SHELL)+modFT_PopNum_yxmsz(yr,x,  MATURE,OLD_SHELL);
        }//x
    }//endyr
//    cout<<"3a"<<endl;
    
    // predicted survey values 
    for (int yr=styr;yr<=endyr;yr++){
        for(int x=1;x<=nSXs;x++) {
            if (yr<1982)             {modSrvNum_xy(x,yr) = (modNum_xyz(x,yr)*selSrv1_xz(x));}  else
            if (1982<=yr && yr<1988) {modSrvNum_xy(x,yr) = (modNum_xyz(x,yr)*selSrv2_xz(x));} else
            if (1988<=yr)            {modSrvNum_xy(x,yr) = (modNum_xyz(x,yr)*selSrv3_xz(x));}
        }
    }
//    cout<<"4"<<endl;
    
    dvariable totSrvNum;
    dvar_matrix useSelSrv(1,nSXs,1,nZBs);
    modPopBio_y.initialize();
    fspbio.initialize();
    mspbio.initialize(); 
    modSrvPrNatZ_msxyz.initialize();
    modSrvPrNatZ_mxyz.initialize();
    for (int yr=styr;yr<=endyr;yr++) {
        fspbio(yr) = natlength_mat(FEMALE,yr)*wt_xmz(FEMALE,MATURE);//dot product sum
        mspbio(yr) = natlength_mat(  MALE,yr)*wt_xmz(  MALE,MATURE);//dot product sum
        
        // Selection pattern
        if (yr<1982)             {useSelSrv = selSrv1_xz;}  else
        if (1982<=yr && yr<1988) {useSelSrv = selSrv2_xz;} else
        if (1988<=yr)            {useSelSrv = selSrv3_xz;}
        
        modSrvImmBio_xy(FEMALE,yr) = multQ*natlength_imm(FEMALE,yr)*elem_prod(wt_xmz(FEMALE,IMMATURE),useSelSrv(FEMALE));//dot product sum
        modSrvImmBio_xy(  MALE,yr) = multQ*natlength_imm(  MALE,yr)*elem_prod(wt_xmz(MALE,  IMMATURE),useSelSrv(  MALE));//dot product sum
        modSrvMatBio_xy(FEMALE,yr) = multQ*natlength_mat(FEMALE,yr)*elem_prod(wt_xmz(FEMALE,  MATURE),useSelSrv(FEMALE));//dot product sum
        modSrvMatBio_xy(  MALE,yr) = multQ*natlength_mat(  MALE,yr)*elem_prod(wt_xmz(MALE,    MATURE),useSelSrv(  MALE));//dot product sum
        
        // this is predicted survey in numbers not biomass-don't adjust by max selectivity 
        for(int x=1;x<=nSXs;x++) modSrvNum_xy(x,yr) = (modNum_xyz(x,yr)*useSelSrv(x));
        totSrvNum = modSrvNum_xy(FEMALE,yr) + modSrvNum_xy(MALE,yr);
        if(totSrvNum<0.001) totSrvNum = 1.0;                     //this is non-differentiable, but PROBABLY just means no survey was done
        for(int x=1;x<=nSXs;x++) {
            for(int m=1;m<=nMSs;m++) {
                for(int s=1;s<=nSCs;s++) {
                    modSrvPrNatZ_msxyz(m,s,x,yr) = elem_prod(useSelSrv(x),modNum_yxmsz(yr,x,m,s))/totSrvNum;
                    modSrvPrNatZ_mxyz(m,x,yr) += modSrvPrNatZ_msxyz(m,s,x,yr);
                }//s
            }//m
        } //x
        modPopBio_y(yr) +=  modNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL)*wt_xmz(FEMALE,IMMATURE)
                          + modNum_yxmsz(yr,  MALE,IMMATURE,NEW_SHELL)*wt_xmz(  MALE,  MATURE)
                          +(modNum_yxmsz(yr,FEMALE,MATURE,NEW_SHELL)+modNum_yxmsz(yr,FEMALE,MATURE,OLD_SHELL))*wt_xmz(FEMALE,MATURE)
                          +(modNum_yxmsz(yr,  MALE,MATURE,NEW_SHELL)+modNum_yxmsz(yr,  MALE,MATURE,OLD_SHELL))*wt_xmz(  MALE,MATURE);
    }//yr  
    //    cout<<" end srv 2"<<endl;
//    cout<<"5"<<endl;
    
//    // Legal males
//    modPopNumLegal_y.initialize();
//    modSrvNumLegal_y.initialize();
//    for (int yr=styr;yr<=endyr;yr++) {
//        // Selection pattern//
//        if (yr<1982)             {useSelSrv = selSrv1_xz;}  else
//        if (1982<=yr && yr<1988) {useSelSrv = selSrv2_xz;} else
//        if (1988<=yr)            {useSelSrv = selSrv3_xz;}        
//        // legal is >=138mm take half the numbers in the 135-139 bin
//        modPopNumLegal_y(yr)    = 0.5*modNum_xyz(MALE,yr,23);
//        modSrvNumLegal_y(yr) = 0.5*modNum_xyz(MALE,yr,23)*useSelSrv(MALE,23);  //fixed indices; need vector of 0's, 0.5 and 1's to mult here
//        for(int j=24;j<=nZBs;j++) {
//            modPopNumLegal_y(yr) += modNum_xyz(MALE,yr,j);
//            modSrvNumLegal_y(yr) += modNum_xyz(MALE,yr,j)*useSelSrv(MALE,j);
//        }
//    }
    //  cout<<"6"<<endl;
    
    //  cout<<" to end of number at len "<<endl;
    //  cout<<"done"<<endl;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
FUNCTION get_catch_at_len                  //wts: new version
//    cout<<"get_catch_at_len"<<endl;
    dvar_vector ratio1(1,nZBs);
    dvar_vector ratio2(1,nZBs);
    
    predRetNumMortTCFM_yz.initialize();
    predTotNumMortTCFM_yz.initialize();
    predDscNumMortTCF_xyz.initialize();
    predDscNumMortSCF_xyz.initialize();
    predDscNumMortRKF_xyz.initialize();
    predDscNumMortGTF_xyz.initialize();
    
    predRetBioMortTCFM_y.initialize();
    predTotBioMortTCFM_y.initialize();
    predDscBioMortTCF_xy.initialize();
    predDscBioMortSCF_xy.initialize();
    predDscBioMortRKF_xy.initialize();
    predDscBioMortGTF_xy.initialize();
    //   cout<<" to get catch at length "<<endl;
    for (int yr=styr;yr<endyr;yr++){                //(IMPORTANT CHANGE: used to be "endyr")        
        //total male directed catch mortality
        ratio1 = elem_prod(elem_div(fmTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iTCF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iTCF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
        predTotNumMortTCFM_syz(NEW_SHELL,yr) = tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL);
        predTotNumMortTCFM_syz(OLD_SHELL,yr) = tmN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL);
        predTotNumMortTCFM_yz(yr) = predTotNumMortTCFM_syz(NEW_SHELL,yr)+predTotNumMortTCFM_syz(OLD_SHELL,yr);
        predTotBioMortTCFM_y(yr) = predTotNumMortTCFM_yz(yr)*wt_xmz(MALE,  MATURE);//note dot product sum over size bins here
        //         cout<<ratio1<<endl;
        //         cout<<ratio2<<endl;
        
        //retained male directed catch mortality
        ratio1 = elem_prod(elem_div(fmTCFR_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmTCFR_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        rmN_ymsz(yr,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
        rmN_ymsz(yr,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
        rmN_ymsz(yr,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
        rmN_ymsz(yr,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
        predRetNumMortTCFM_syz(NEW_SHELL,yr) = rmN_ymsz(yr,IMMATURE,NEW_SHELL) + rmN_ymsz(yr,MATURE,NEW_SHELL); 
        predRetNumMortTCFM_syz(OLD_SHELL,yr) = rmN_ymsz(yr,IMMATURE,OLD_SHELL) + rmN_ymsz(yr,MATURE,OLD_SHELL);
        predRetNumMortTCFM_yz(yr) = predRetNumMortTCFM_syz(NEW_SHELL,yr)+predRetNumMortTCFM_syz(OLD_SHELL,yr);
        predRetBioMortTCFM_y(yr)  = predRetNumMortTCFM_yz(yr)*wt_xmz(MALE,  MATURE);
        
        predDscNumMortTCF_xyz(MALE,yr) = predTotNumMortTCFM_yz(yr) - predRetNumMortTCFM_yz(yr);
        predDscBioMortTCF_xy(MALE,yr)  = predTotBioMortTCFM_y(yr)  - predRetBioMortTCFM_y(yr);
        
        //snow crab discard catch male mortality
        ratio1 = elem_prod(elem_div(fmSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iSCF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iSCF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
        predDscNumMortSCF_xyz(MALE,yr) = tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iSCF,yr,MALE,MATURE,NEW_SHELL)
                                        +tmN_fyxmsz(iSCF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iSCF,yr,MALE,MATURE,OLD_SHELL);
        predDscBioMortSCF_xy(MALE,yr) = predDscNumMortSCF_xyz(MALE)(yr)*wt_xmz(MALE,  MATURE);
        
        //red king crab discard catch male     
        ratio1 = elem_prod(elem_div(fmRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iRKF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iRKF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
        predDscNumMortRKF_xyz(MALE,yr) = tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iRKF,yr,MALE,MATURE,NEW_SHELL)
                                        +tmN_fyxmsz(iRKF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iRKF,yr,MALE,MATURE,OLD_SHELL);
        predDscBioMortRKF_xy(MALE,yr) = predDscNumMortRKF_xyz(MALE)(yr)*wt_xmz(MALE,  MATURE);
        
        //trawl bycatch male (20150601: is this correct?!->rates are indep of shell condition)
        ratio1 = elem_prod(elem_div(fmGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
        tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iGTF,yr,MALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iGTF,yr,MALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
        predDscNumMortGTF_xyz(MALE,yr) = tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iGTF,yr,MALE,MATURE,NEW_SHELL) 
                                        +tmN_fyxmsz(iGTF,yr,MALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iGTF,yr,MALE,MATURE,OLD_SHELL); 
        predDscBioMortGTF_xy(MALE,yr) = predDscNumMortGTF_xyz(MALE,yr)*wt_xmz(MALE,  MATURE);
        
        //directed tanner discard catch female
        ratio1 = elem_prod(elem_div(fmTCFF_yz(yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmTCFF_yz(yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iTCF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iTCF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
        predDscNumMortTCF_xyz(FEMALE,yr) = tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iTCF,yr,FEMALE,MATURE,NEW_SHELL)
                                          +tmN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iTCF,yr,FEMALE,MATURE,OLD_SHELL);
        predDscBioMortTCF_xy(FEMALE,yr)  = predDscNumMortTCF_xyz(FEMALE,yr)*wt_xmz(FEMALE,MATURE);//all crab assumed mature
        
        //snow crab discard catch female
        ratio1 = elem_prod(elem_div(fmSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iSCF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iSCF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
        predDscNumMortSCF_xyz(FEMALE,yr) = tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iSCF,yr,FEMALE,MATURE,NEW_SHELL)
                                          +tmN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iSCF,yr,FEMALE,MATURE,OLD_SHELL);
        predDscBioMortSCF_xy(FEMALE,yr) = predDscNumMortSCF_xyz(FEMALE)(yr)*wt_xmz(FEMALE,MATURE);
        
        //red king crab discard catch female
        ratio1 = elem_prod(elem_div(fmRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iRKF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iRKF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
        predDscNumMortRKF_xyz(FEMALE,yr) = tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iRKF,yr,FEMALE,MATURE,NEW_SHELL)
                                          +tmN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iRKF,yr,FEMALE,MATURE,OLD_SHELL);
        predDscBioMortRKF_xy(FEMALE,yr) = predDscNumMortRKF_xyz(FEMALE,yr)*wt_xmz(FEMALE,MATURE);
        
        //trawl bycatch female (20150601: explicitly calculating ratio 2 rather than assuming ratio2=ratio1)
        ratio1 = elem_prod(elem_div(fmGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
        ratio2 = elem_prod(elem_div(fmGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
        tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
        tmN_fyxmsz(iGTF,yr,FEMALE,  MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
        tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
        tmN_fyxmsz(iGTF,yr,FEMALE,  MATURE,OLD_SHELL) = elem_prod(ratio2,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
        predDscNumMortGTF_xyz(FEMALE,yr) = tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,NEW_SHELL) + tmN_fyxmsz(iGTF,yr,FEMALE,MATURE,NEW_SHELL) 
                                          +tmN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,OLD_SHELL) + tmN_fyxmsz(iGTF,yr,FEMALE,MATURE,OLD_SHELL); 
        predDscBioMortGTF_xy(FEMALE,yr) = predDscNumMortGTF_xyz(FEMALE,yr)*wt_xmz(FEMALE,MATURE);//add in females to males for trawl catches
    }    //end of year loop
    
    modPrNatZ_TCFR_syz.initialize();
    modPrNatZ_TCFM_syz.initialize();
    modPrNatZ_TCFF_yz.initialize();
    modPrNatZ_SCF_xyz.initialize();
    modPrNatZ_RKF_xyz.initialize();
    modPrNatZ_GTF_xyz.initialize();
    for (int yr=styr;yr<endyr;yr++) {          //(IMPORTANT CHANGE: used to be "endyr")        
        // Retained catch (males)
        if(sum(predRetNumMortTCFM_yz(yr))>0.000001){                         //non-differentiable--does it matter
            dvariable tot = sum(predRetNumMortTCFM_yz(yr));
            modPrNatZ_TCFR_syz(NEW_SHELL,yr) = predRetNumMortTCFM_syz(NEW_SHELL,yr)/tot;
            modPrNatZ_TCFR_syz(OLD_SHELL,yr) = predRetNumMortTCFM_syz(OLD_SHELL,yr)/tot;
        }
        
        // Total catch (males)
        if(sum(predTotNumMortTCFM_yz(yr))>0.0000001){                         //non-differentiable--does it matter
            dvariable tot = sum(predTotNumMortTCFM_yz(yr));
            modPrNatZ_TCFM_syz(NEW_SHELL,yr) = predTotNumMortTCFM_syz(NEW_SHELL,yr)/tot;
            modPrNatZ_TCFM_syz(OLD_SHELL,yr) = predTotNumMortTCFM_syz(OLD_SHELL,yr)/tot;
        }
         
        
        // female discards
        if(sum(predDscNumMortTCF_xyz(FEMALE,yr))>0.0000001){                         //non-differentiable--does it matter
            modPrNatZ_TCFF_yz(yr) = predDscNumMortTCF_xyz(FEMALE,yr)/sum(predDscNumMortTCF_xyz(FEMALE,yr));
        }
        
        // snow crab discards female male
        if(sum(predDscNumMortSCF_xyz(FEMALE,yr))>0.00000001){                    //non-differentiable--does it matter
            modPrNatZ_SCF_xyz(FEMALE,yr) = predDscNumMortSCF_xyz(FEMALE,yr)/sum(predDscNumMortSCF_xyz(FEMALE,yr));
        }
        if(sum(predDscNumMortSCF_xyz(MALE,yr))>0.00000001){                       //non-differentiable--does it matter
            modPrNatZ_SCF_xyz(MALE,yr) = predDscNumMortSCF_xyz(MALE,yr)/sum(predDscNumMortSCF_xyz(MALE,yr));
        }
        
        // red king crab discards female male
        if(sum(predDscNumMortRKF_xyz(FEMALE,yr))>0.00000001){                       //non-differentiable--does it matter
            modPrNatZ_RKF_xyz(FEMALE,yr) = predDscNumMortRKF_xyz(FEMALE,yr)/sum(predDscNumMortRKF_xyz(FEMALE,yr));
        }
        if(sum(predDscNumMortRKF_xyz(MALE,yr))>0.00000001){                         //non-differentiable--does it matter
            modPrNatZ_RKF_xyz(MALE,yr) = predDscNumMortRKF_xyz(MALE,yr)/sum(predDscNumMortRKF_xyz(MALE,yr));
        }
        
        // total trawl selected numbers
        {
            dvariable tot = sum(predDscNumMortGTF_xyz(MALE,yr))+sum(predDscNumMortGTF_xyz(FEMALE,yr));
            // Trawl proportions (adds to 1 over sex, shell and length)
            modPrNatZ_GTF_xyz(FEMALE,yr) = predDscNumMortGTF_xyz(FEMALE,yr)/tot;
            modPrNatZ_GTF_xyz(  MALE,yr) = predDscNumMortGTF_xyz(  MALE,yr)/tot;
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
//    if (active(pRecDevs)) {        
        //recruitment likelihood - norm2 is sum of square values   
        penal_rec = 1.0*like_wght_rec*norm2(pRecDevs);
        //   penalty on devs in historic period     
        //penal_rec += 1.0*norm2(pRecDevsHist);
        //   penalty on dev first differences in historic period     
        penal_rec += 1.0*norm2(first_difference(pRecDevsHist));
        
        f += penal_rec; objfOut(1) = penal_rec; likeOut(1) = penal_rec; wgtsOut(1) = 1;
//    } 
    
    //nat Mort. penalty
//    if(active(pMfac_Imm)) {
        nat_penalty = 0.5 * square((pMfac_Imm - 1.0) / 0.05);                 //hard-wired
        f += nat_penalty; objfOut(3) = nat_penalty; likeOut(3) = nat_penalty; wgtsOut(3) = 1;
//    }
//    if(active(pMfac_MatM)) {  
        nat_penalty = 0.5 * square((pMfac_MatM - 1.0) / 0.05);                     //hard-wired
        f += nat_penalty; objfOut(4) = nat_penalty; likeOut(4) = nat_penalty; wgtsOut(4) = 1;
//    }
//    if(active(pMfac_MatF)) {  
        nat_penalty = 0.5 * square((pMfac_MatF - 1.0) / 0.05);                     //hard-wired
        f += nat_penalty; objfOut(5) = nat_penalty; likeOut(5) = nat_penalty; wgtsOut(5) = 1;
//    }
    
    //penalty on survey Q
//    if(active(pSrv2_QM) && srv3_qPriorWgt>0) {  
//        //max of underbag at 182.5 mm is 0.873   
//        srv3q_penalty = 0.5 * square((pSrv2_QM - 0.88) / 0.05);                  //hard-wired
//        //    srv3q_penalty = 0.0 * square((pSrv2_QM - 0.88) / 0.05);
//        f += srv3q_penalty; objfOut(6) = srv3q_penalty; likeOut(6) = srv3q_penalty; wgtsOut(6) = 1;
//    }
//    if(active(pSrv2_QF) && srv3_qPriorWgt>0) {  
//        //peak of females is at about 80mm underbag is 0.75 at this size - less uncertainty  
//        srv3q_penalty = 0.5 * square((pSrv2_QF - 0.88) / 0.05);                //hard-wired
//        //    srv3q_penalty = 0.0 * square((pSrv2_QF - 0.88) / 0.05);
//        f += srv3q_penalty; objfOut(7) = srv3q_penalty; likeOut(7) = srv3q_penalty; wgtsOut(7) = 1;
//    }
//    if(active(pSrv2_QM) && srv3_qPriorWgt>=0) {  
        //max of underbag at 182.5 mm is 0.873   
        srv3q_penalty = 0.5 * square((pSrv2_QM - srv3_qPriorMean) / srv3_qPriorStD);
        f += srv3_qPriorWgt*srv3q_penalty; objfOut(6) = srv3_qPriorWgt*srv3q_penalty; likeOut(6) = srv3q_penalty; wgtsOut(6) = srv3_qPriorWgt;
//    }
//    if(active(pSrv2_QF) && srv3_qFemPriorWgt>=0) {  
        //peak of females is at about 80mm underbag is 0.75 at this size - less uncertainty  
        srv3q_penalty = 0.5 * square((pSrv2_QF - srv3_qFemPriorMean) / srv3_qFemPriorStD);
        f += srv3_qFemPriorWgt*srv3q_penalty; objfOut(7) = srv3_qFemPriorWgt*srv3q_penalty; likeOut(7) = srv3q_penalty; wgtsOut(7) = srv3_qFemPriorWgt;
//    }
    
    // bayesian part - likelihood on growth parameters af,am,bf,bm
    // not used in this case
    af_penal = 0; bf_penal = 0; am_penal = 0; bm_penal = 0;
//    if(active(pGrAF1)) {  
        af_penal = 0.5 * square((pGrAF1 - 0.56560241)/0.1);                  //hard-wired
        f += af_penal; objfOut(8) = af_penal; likeOut(8) = af_penal; wgtsOut(8) = 1;
//    }
//    if(active(pGrBF1)) {  
        bf_penal = 0.5 * square((pGrBF1 - 0.9132661)/0.025);                    //hard-wired
        f += bf_penal; objfOut(9) = bf_penal; likeOut(9) = bf_penal; wgtsOut(9) = 1;
//    }
//    if(active(pGrAM1)) {
        am_penal   = 0.5 * square((pGrAM1 - 0.437941)/0.025);                     //hard-wired
        f += am_penal; objfOut(10) = am_penal; likeOut(10) = am_penal; wgtsOut(10) = 1;
//    }
//    if(active(pGrBM1)) {
        bm_penal = 0.5 * square((pGrBM1 - 0.9487)/0.1);                          //hard-wired
        f += bm_penal; objfOut(11) = bm_penal; likeOut(11) = bm_penal; wgtsOut(11) = 1;
//    }
    
//    if(active(pPrM2MF)) {
        lkPrM2M = norm2(first_difference(first_difference(pPrM2MF)));
        f += 1.0*lkPrM2M; objfOut(12)= 1.0*lkPrM2M; likeOut(12) = lkPrM2M; wgtsOut(12) = 1;
//    }
        
//    if(active(pPrM2MM)) {
        lkPrM2M = norm2(first_difference(first_difference(pPrM2MM)));//wts: why different weights? 1.0 vs. 0.5?
        f += 0.5*lkPrM2M; objfOut(13)= 0.5*lkPrM2M; likeOut(13) = lkPrM2M; wgtsOut(13) = 0.5;
//    }
    
    // various penalties
    // =================
    double llw = 0.0;
    double red = 0.01;
    fpen.initialize();
//    if (active(pSelTCFM_devsZ50)) 
    { 
        int phs = pSelTCFM_devsZ50.get_phase_start();
        llw = llwSelTCFM_devsZ50; 
//        if (doPenRed) llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))*llw;
        nextf = norm2(first_difference(pSelTCFM_devsZ50));
        fpen += llw*nextf; objfOut(14) = llw*nextf; likeOut(14) = nextf; wgtsOut(14) = llw;           
        
        phs = pSelTCFM_devsZ50.get_phase_start();
        llw = llwSelTCFM_devsZ50; 
//        if (doPenRed) llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))*llw;
        nextf = norm2(pSelTCFM_devsZ50);
        fpen += llw*nextf; objfOut(38) = llw*nextf; likeOut(38) = nextf; wgtsOut(38) = llw;   //wts: need to turn this off in last phase?        
    }
//    if (active(pF_DevsTCF)) 
    { 
        int phs = pF_DevsTCF.get_phase_start();
        llw = 1.0; 
        if (doPenRed) {
            llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))*llw;
            if (current_phase()>=maxPhase) llw = 0.0;
        }
        nextf = norm2(pF_DevsTCF);
        fpen += llw*nextf; objfOut(15) = llw*nextf; likeOut(15) = nextf; wgtsOut(15) = llw;   //wts: need to turn this off in last phase?        
    }
//    if(active(pF_DevsSCF)) 
    {
        int phs = pF_DevsSCF.get_phase_start();
        llw = 0.5; 
        if (doPenRed) {
            llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))*llw;
            if (current_phase()>=maxPhase) llw = 0.0;
        }
        nextf = norm2(pF_DevsSCF);
        fpen += llw*nextf; objfOut(16) = llw*nextf; likeOut(16) = nextf; wgtsOut(16) = llw; //wts: need to turn this off in last phase? note that relative weights are hard-wired
        
    }
//    if(active(pF_DevsRKF)) 
    {
        int phs = pF_DevsRKF.get_phase_start();
        llw = 3.0; 
        if (doPenRed) {
            llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))*llw;
            if (current_phase()>=maxPhase) llw = 0.0;
        }
        nextf = norm2(pF_DevsRKF);
        fpen += llw*nextf; objfOut(17) = llw*nextf; likeOut(17) = nextf; wgtsOut(17) = llw; //wts: need to turn this off in last phase?
        
    }
//    if(active(pF_DevsGTF)) 
    {
        int phs = pF_DevsGTF.get_phase_start();
        llw = 0.5; 
        if (doPenRed) {
            llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))*llw;
            if (current_phase()>=maxPhase) llw = 0.0;
        }
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
        yr=yrsObsZCsGTF_n(n);
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
            lkZCs( 9) -= ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,x,n)
                          *obsSrvPrNatZ_mxnz(IMMATURE,x,n)*log(modSrvPrNatZ_mxyz(IMMATURE,x,yr)+p_const);
            // this is for mature new and old shell together
            lkZCs(10) -= ssObsZCsSrv_msxn(  MATURE,NEW_SHELL,x,n)
                          *obsSrvPrNatZ_mxnz(  MATURE,x,n)*log(modSrvPrNatZ_mxyz(  MATURE,x,yr)+p_const);       

            x = FEMALE;   
            // obs(maturity, SC, x, year), pred(maturity,x, year)
            // immature new and old together
            lkZCs(11) -= ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,x,n)
                          *obsSrvPrNatZ_mxnz(IMMATURE,x,n)*log(modSrvPrNatZ_mxyz(IMMATURE,x,yr)+p_const);
            // this is for mature new and old shell together
            lkZCs(12) -= ssObsZCsSrv_msxn(  MATURE,NEW_SHELL,x,n)
                          *obsSrvPrNatZ_mxnz(  MATURE,x,n)*log(modSrvPrNatZ_mxyz(  MATURE,x,yr)+p_const);        
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
        zsRetMortBio_TCFR_y(1965,endyr-1) = (obsRetCatchBio_y(1965,endyr-1)-predRetBioMortTCFM_y(1965,endyr-1))/sqrt(0.5);
        lkRetMortBio_TCFR                  = 0.5*norm2(zsRetMortBio_TCFR_y);
        nextf = like_wght_CatchBio*lkRetMortBio_TCFR; objfOut(32) = nextf; f += nextf; likeOut(32) = lkRetMortBio_TCFR; wgtsOut(32) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrTCFR)));
        zsRetMortBio_TCFR_y(1965,endyr-1) = (log(obsRetCatchBio_y(1965,endyr-1)+smlValFsh)-log(predRetBioMortTCFM_y(1965,endyr-1)+smlValFsh))/stdv;
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
            zsTotMortBio_TCFM_n = (obsTotBioMortTCFM_n - predTotBioMortTCFM_y(yrsObsDscTCF_n))/sqrt(0.5);
            lkTotMortBio_TCFM = 0.5*norm2(zsTotMortBio_TCFM_n);
            nextf = like_wght_CatchBio*lkTotMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkTotMortBio_TCFM; wgtsOut(33) = like_wght_CatchBio;
        } else {
            //lognormal error assumption
            double stdv = sqrt(log(1.0+square(obsErrTCFD)));
            zsTotMortBio_TCFM_n = (log(obsTotBioMortTCFM_n+smlValFsh)-log(predTotBioMortTCFM_y(yrsObsDscTCF_n)+smlValFsh))/stdv;
            lkTotMortBio_TCFM = 0.5*norm2(zsTotMortBio_TCFM_n);
            nextf = lkTotMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkTotMortBio_TCFM; wgtsOut(33) = 1.0;
        }
    } else if (optTCFMfit==1) {
        //fit to discard male mortality
        if (optFshNLLs==0){
            //normal error assumption, sd=sqrt(2) [fixed]
            zsDscMortBio_TCFM_n = (obsDscBioMortTCF_xn(MALE) - predDscBioMortTCF_xy(MALE)(yrsObsDscTCF_n))/sqrt(0.5);
            lkDscMortBio_TCFM = 0.5*norm2(zsDscMortBio_TCFM_n);
            nextf = like_wght_CatchBio*lkDscMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkDscMortBio_TCFM; wgtsOut(33) = like_wght_CatchBio;
        } else {
            //lognormal error assumption
            double stdv = sqrt(log(1.0+square(obsErrTCFD)));
            zsDscMortBio_TCFM_n = (log(obsDscBioMortTCF_xn(MALE)+smlValFsh)-log(predDscBioMortTCF_xy(MALE)(yrsObsDscTCF_n)+smlValFsh))/stdv;
            lkDscMortBio_TCFM = 0.5*norm2(zsDscMortBio_TCFM_n);
            nextf = lkDscMortBio_TCFM; objfOut(33) = nextf; f += nextf; likeOut(33) = lkDscMortBio_TCFM; wgtsOut(33) = 1.0;
        }
    }
    
    //  female catch in directed fishery
    zsDscMortBio_TCFF_n.initialize();
    lkDscMortBio_TCFF.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(2) [fixed]
        zsDscMortBio_TCFF_n = (obsDscBioMortTCF_xn(FEMALE)-predDscBioMortTCF_xy(FEMALE)(yrsObsDscTCF_n))/sqrt(0.5);
        lkDscMortBio_TCFF = 0.5*norm2(zsDscMortBio_TCFF_n);
        nextf = like_wght_CatchBio*lkDscMortBio_TCFF; objfOut(34) = nextf; f+= nextf; likeOut(34) = lkDscMortBio_TCFF; wgtsOut(34) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrTCFD)));
        zsDscMortBio_TCFF_n = (log(obsDscBioMortTCF_xn(FEMALE)+smlValFsh)-log(predDscBioMortTCF_xy(FEMALE)(yrsObsDscTCF_n)+smlValFsh))/stdv;
        lkDscMortBio_TCFF = 0.5*norm2(zsDscMortBio_TCFF_n);
        nextf = lkDscMortBio_TCFF; objfOut(34) = nextf; f+= nextf; likeOut(34) = lkDscMortBio_TCFF; wgtsOut(34) = 1.0;
    }
    
    //snow crab fishery
    zsDscMortBio_SCF_xn.initialize();
    lkDscMortBio_SCF_x.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(2) [fixed]
        for (int x=1;x<=nSXs;x++){
            zsDscMortBio_SCF_xn(x) = (obsDscBioMortSCF_xn(x) - predDscBioMortSCF_xy(x)(yrsObsDscSCF))/sqrt(0.5);
            lkDscMortBio_SCF_x(x)  = 0.5*norm2(zsDscMortBio_SCF_xn(x));
        }
        nextf = like_wght_CatchBio*sum(lkDscMortBio_SCF_x); objfOut(35) = nextf; f += nextf; likeOut(35) = sum(lkDscMortBio_SCF_x); wgtsOut(35) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrSCF)));
        for (int x=1;x<=nSXs;x++){
            zsDscMortBio_SCF_xn(x) = (log(obsDscBioMortSCF_xn(x)+smlValFsh)-log(predDscBioMortSCF_xy(x)(yrsObsDscSCF)+smlValFsh))/stdv;
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
            zsDscMortBio_RKF_xn(x) = (obsDscBioMortRKF_xn(x) - predDscBioMortRKF_xy(x)(yrsObsDscRKF))/sqrt(0.5);
            lkDscMortBio_RKF_x(x)  = 0.5*norm2(zsDscMortBio_RKF_xn(x));
        }
        nextf = like_wght_CatchBio*sum(lkDscMortBio_RKF_x); objfOut(36) = nextf; f += nextf; likeOut(36) = sum(lkDscMortBio_RKF_x); wgtsOut(36) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrRKF)));
        for (int x=1;x<=nSXs;x++){
            zsDscMortBio_RKF_xn(x) = (log(obsDscBioMortRKF_xn(x)+smlValFsh)-log(predDscBioMortRKF_xy(x)(yrsObsDscRKF)+smlValFsh))/stdv;
            lkDscMortBio_RKF_x(x)  = 0.5*norm2(zsDscMortBio_RKF_xn(x));
        }
        nextf = sum(lkDscMortBio_RKF_x); objfOut(36) = nextf; f += nextf; likeOut(36) = sum(lkDscMortBio_RKF_x); wgtsOut(36) = 1.0;
    }
    
    //groundfish trawl fishery
    zsDscMortBio_GTF_n.initialize();
    lkDscMortBio_GTF.initialize();
    if (optFshNLLs==0){
        //normal error assumption, sd=sqrt(1/2) [fixed]
        zsDscMortBio_GTF_n = (obsDscBioMortGTF_n - (predDscBioMortGTF_xy(FEMALE)+predDscBioMortGTF_xy(MALE))(yrsObsDscGTF))/sqrt(0.5);
        lkDscMortBio_GTF = 0.5*norm2(zsDscMortBio_GTF_n);
        nextf = like_wght_CatchBio*lkDscMortBio_GTF;  objfOut(37) = nextf; f += nextf; likeOut(37) = lkDscMortBio_GTF; wgtsOut(37) = like_wght_CatchBio;
    } else {
        //lognormal error assumption
        double stdv = sqrt(log(1.0+square(obsErrGTF)));
        zsDscMortBio_GTF_n = (log(obsDscBioMortGTF_n+smlValFsh) - log((predDscBioMortGTF_xy(FEMALE)+predDscBioMortGTF_xy(MALE))(yrsObsDscGTF)+smlValFsh))/stdv;
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
//Effective N for size comps using the McAllister-Ianelli approach
FUNCTION double effN_McIan(dvector& obsPr, dvar_vector& modPr)
    dvariable num = modPr*(1-modPr);//note dot product
    dvariable den = norm2(obsPr-modPr);
    double effN = value(num/den);
    return effN;
    
// ==========================================================================
// ==========================================================================
//Effective z-score for size comps using the Francis mean size approach
FUNCTION double effZ_Francis(dvector& obsPr, dvar_vector& modPr, double inpSS)
    double obsMnZ = zBs*obsPr;
    double modMnZ = value(zBs*modPr);
    double modSeZ = sqrt(value(modPr*elem_prod(zBs-modMnZ,zBs-modMnZ))/inpSS);
    double effZ = (obsMnZ-modMnZ)/modSeZ;
    return effZ;
    
// ==========================================================================
// ==========================================================================
FUNCTION Misc_output
    
    int debug = 0;
    if (debug) cout<<" to Misc_output()"<<endl;

    dvar_matrix useSelSrv(1,nSXs,1,nZBs);
    dvar_matrix cv_srv1_nowt(1,nSXs,styr,endyr);
    
    //legal males at time of fishery
    modFT_PopNumLegal_y.initialize();
    modFT_PopBioLegal_y.initialize();
    modTotBioMortLegal_TCFM_y.initialize();
    modTotNumMortLegal_TCFM_y.initialize();
    for (int yr=styr;yr<endyr;yr++) {      //(IMPORTANT CHANGE: used to be "endyr")
        if (yr<1982)             {useSelSrv = selSrv1_xz;}  else
        if (1982<=yr && yr<1988) {useSelSrv = selSrv2_xz;} else
        if (1988<=yr)            {useSelSrv = selSrv3_xz;}        
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                modFT_PopNumLegal_y(yr) += sum(modFT_PopNum_yxmsz(yr,MALE,m,s)(iZLegal,nZBs));
                modFT_PopBioLegal_y(yr) +=     modFT_PopNum_yxmsz(yr,MALE,m,s)(iZLegal,nZBs)*wt_xmz(MALE,m)(iZLegal,nZBs);
                modTotNumMortLegal_TCFM_y(yr) += sum(tmN_fyxmsz(iTCF,yr,MALE,m,s)(iZLegal,nZBs)); 
                modTotBioMortLegal_TCFM_y(yr) +=     tmN_fyxmsz(iTCF,yr,MALE,m,s)(iZLegal,nZBs)*wt_xmz(MALE,m)(iZLegal,nZBs);
           }//s
       }//m
    }//yr loop
    
    if (debug) cout<<" to eff N "<<endl;
    // Effective N's                    //can these be vectorized?
    dvariable num; dvariable den;
    //retained male size comps in directed fishery
    effnTCF_ret_y.initialize();
    for (int i=1;i<=nObsRetZCsTCF;i++){
        int yr=yrsObsRetZCsTCF_n(i);
        num = (modPrNatZ_TCFR_syz(NEW_SHELL,yr)+modPrNatZ_TCFR_syz(OLD_SHELL,yr))
                *(1-(modPrNatZ_TCFR_syz(NEW_SHELL,yr)+modPrNatZ_TCFR_syz(OLD_SHELL,yr)));
        den = norm2((obsPrNatZ_TCFR_sn(NEW_SHELL,i)+obsPrNatZ_TCFR_sn(OLD_SHELL,i))
                      -(modPrNatZ_TCFR_syz(NEW_SHELL,yr)+modPrNatZ_TCFR_syz(OLD_SHELL,yr)));
        effnTCF_ret_y(yr) = num/den;
    }//i loop
    //total catch size comps in directed fishery
    if (debug) cout<<" to eff N: TCF tot "<<endl;
    effnTCF_tot_xy.initialize();
    for (int i=1;i<=nObsZCsTCFM;i++) {
        int yr=yrsObsZCsTCFM_n(i);
        num = (modPrNatZ_TCFM_syz(NEW_SHELL,yr)+modPrNatZ_TCFM_syz(OLD_SHELL,yr))
                *(1-(modPrNatZ_TCFM_syz(NEW_SHELL,yr)+modPrNatZ_TCFM_syz(OLD_SHELL,yr)));
        den = norm2((obsPrNatZ_TCFM_snz(NEW_SHELL,i)+obsPrNatZ_TCFM_snz(NEW_SHELL,i))
                      -(modPrNatZ_TCFM_syz(NEW_SHELL,yr)+modPrNatZ_TCFM_syz(OLD_SHELL,yr)));
        effnTCF_tot_xy(  MALE,yr) = num/den;
        num = modPrNatZ_TCFF_yz(yr)*(1-modPrNatZ_TCFF_yz(yr));
        den = norm2(obsPrNatZ_TCFF_nz(i)-modPrNatZ_TCFF_yz(yr));
        effnTCF_tot_xy(FEMALE,yr) = num/den;
    }// i loop
    //SCF size comps
    if (debug) cout<<" to eff N: SCF tot "<<endl;
    effnSCF_tot_xy.initialize();
    for (int i=1;i<=nObsZCsSCF;i++) {
        int yr=yrsObsZCsSCF_n(i);
        for (int x=1;x<=nSXs;x++) {
            num = modPrNatZ_SCF_xyz(x,yr)*(1-modPrNatZ_SCF_xyz(x,yr));
            den = norm2(obsPrNatZ_SCF_xnz(x,i)-modPrNatZ_SCF_xyz(x,yr));
            effnSCF_tot_xy(x,yr) = num/den;
        }
    }// i loop
    //RKF size comps
    if (debug) cout<<" to eff N: RKF tot "<<endl;
    effnRKF_tot_xy.initialize();
    for (int i=1;i<=nObsZCsRKF;i++) {
        int yr=yrsObsZCsRKF_n(i);
        for (int x=1;x<=nSXs;x++) {
            num = modPrNatZ_RKF_xyz(x,yr)*(1-modPrNatZ_RKF_xyz(x,yr));
            den = norm2(obsPrNatZ_RKF_xnz(x,i)-modPrNatZ_RKF_xyz(x,yr));
            effnRKF_tot_xy(x,yr) = num/den;
        }
    }// i loop
    //GTF size comps
    if (debug) cout<<" to eff N: GTF tot "<<endl;
    effnGTF_tot_y.initialize();
    for (int i=1;i<=nObsZCsGTF;i++) {
        int yr=yrsObsZCsGTF_n(i);
        num = modPrNatZ_GTF_xyz(FEMALE,yr)*(1-modPrNatZ_GTF_xyz(FEMALE,yr))+modPrNatZ_GTF_xyz(MALE,yr)*(1-modPrNatZ_GTF_xyz(MALE,yr));
        den = norm2(obsPrNatZ_GTF_xnz(FEMALE,i)-modPrNatZ_GTF_xyz(FEMALE,yr))+norm2(obsPrNatZ_GTF_xnz(MALE,i)-modPrNatZ_GTF_xyz(MALE,yr));
        effnGTF_tot_y(yr) = num/den;
    }// i loop
    
    //survey size comps
    if (debug) cout<<" to eff N: survey "<<endl;
    effnSrv_y.initialize();
    for (int i=1;i<=nObsZCsSrv;i++) {
        int yr=yrsObsZCsSrv_n(i);
        num.initialize(); den.initialize();
        for (int m=1;m<=nMSs;m++){
            for (int x=1;x<=nSXs;x++){
                num += modSrvPrNatZ_mxyz(m,x,yr)*(1-modSrvPrNatZ_mxyz(m,x,yr));
                den += norm2(obsSrvPrNatZ_mxnz(m,x,i)-modSrvPrNatZ_mxyz(m,x,yr));
            } // x loop
        } // m loop
        effnSrv_y(yr) = num/den;
    } // i loop
    
    // spawning biomass and related outputs
    modSpNumMateTime_xsy.initialize();
    modSpNumMateTime_xy.initialize();
    modSpBioMateTime_xy.initialize();
    mspbio_old_matetime.initialize();
    fspbio_new_matetime.initialize();
 
    if (debug) cout<<" to sp bio matetime "<<endl;
    //spawning occurs AFTER fisheries, so cannot be evaluated in endyr
    dvar_vector tmpNS(1,nZBs);
    dvar_vector tmpOS(1,nZBs);
    for (int yr=styr;yr<endyr;yr++) {  //IMPORTANT CHANGE: was <=endyr
        if (debug) cout<<"yr = "<<yr<<endl;

        if (yr>=lyr_mort && yr<=uyr_mort && mort_switch==1){
            tmpNS = elem_prod(S_xsyz(MALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,MALE)*pMfac_Big(MALE))*modNum_yxmsz(yr,MALE,MATURE,NEW_SHELL));
            tmpOS = elem_prod(S_xsyz(MALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,MALE)*pMfac_Big(MALE))*modNum_yxmsz(yr,MALE,MATURE,OLD_SHELL));
            modSpNumMateTime_xsy(MALE,NEW_SHELL,yr) = sum(tmpNS);
            modSpNumMateTime_xsy(MALE,OLD_SHELL,yr) = sum(tmpOS);
            modSpNumMateTime_xy( MALE,yr) = modSpNumMateTime_xsy(MALE,NEW_SHELL,yr)+modSpNumMateTime_xsy(MALE,OLD_SHELL,yr);
            mspbio_old_matetime(yr)       = tmpOS*wt_xmz(MALE,  MATURE);         //dot product
            modSpBioMateTime_xy( MALE,yr) = (tmpNS+tmpOS)*wt_xmz(MALE,  MATURE); //dot product
            
            tmpNS = elem_prod(S_xsyz(FEMALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,FEMALE)*pMfac_Big(FEMALE))*modNum_yxmsz(yr,FEMALE,MATURE,NEW_SHELL));
            tmpOS = elem_prod(S_xsyz(FEMALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,FEMALE)*pMfac_Big(FEMALE))*modNum_yxmsz(yr,FEMALE,MATURE,OLD_SHELL));
            modSpNumMateTime_xsy(FEMALE,NEW_SHELL,yr) = sum(tmpNS);
            modSpNumMateTime_xsy(FEMALE,OLD_SHELL,yr) = sum(tmpOS);
            modSpNumMateTime_xy( FEMALE,yr) = modSpNumMateTime_xsy(FEMALE,NEW_SHELL,yr)+modSpNumMateTime_xsy(FEMALE,OLD_SHELL,yr);
            fspbio_new_matetime(yr)         = tmpNS*wt_xmz(FEMALE)(MATURE);     //dot product
            modSpBioMateTime_xy( FEMALE,yr) = (tmpNS+tmpOS)*wt_xmz(FEMALE)(MATURE); //dot product
        } else {
            tmpNS = elem_prod(S_xsyz(MALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,MALE))*modNum_yxmsz(yr,MALE,MATURE,NEW_SHELL));
            tmpOS = elem_prod(S_xsyz(MALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,MALE))*modNum_yxmsz(yr,MALE,MATURE,OLD_SHELL));
            modSpNumMateTime_xsy(MALE,NEW_SHELL,yr) = sum(tmpNS);
            modSpNumMateTime_xsy(MALE,OLD_SHELL,yr) = sum(tmpOS);
            modSpNumMateTime_xy(MALE,yr) = modSpNumMateTime_xsy(MALE,NEW_SHELL,yr)+modSpNumMateTime_xsy(MALE,OLD_SHELL,yr);
            mspbio_old_matetime(yr)      = tmpOS*wt_xmz(MALE,  MATURE);         //dot product
            modSpBioMateTime_xy(MALE,yr) = (tmpNS+tmpOS)*wt_xmz(MALE,  MATURE); //dot product
            
            tmpNS = elem_prod(S_xsyz(FEMALE,NEW_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,FEMALE))*modNum_yxmsz(yr,FEMALE,MATURE,NEW_SHELL));
            tmpOS = elem_prod(S_xsyz(FEMALE,OLD_SHELL,yr),mfexp(-(spmo+mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,FEMALE))*modNum_yxmsz(yr,FEMALE,MATURE,OLD_SHELL));
            modSpNumMateTime_xsy(FEMALE,NEW_SHELL,yr) = sum(tmpNS);
            modSpNumMateTime_xsy(FEMALE,OLD_SHELL,yr) = sum(tmpOS);
            modSpNumMateTime_xy(FEMALE,yr) = modSpNumMateTime_xsy(FEMALE,NEW_SHELL,yr)+modSpNumMateTime_xsy(FEMALE,OLD_SHELL,yr);
            fspbio_new_matetime(yr)        = tmpNS*wt_xmz(FEMALE)(MATURE);        //dot product
            modSpBioMateTime_xy(FEMALE,yr) = (tmpNS+tmpOS)*wt_xmz(FEMALE)(MATURE);//dot product
        }
        // effective sp numbers
        emspbio_matetime(yr) = mspbio_old_matetime(yr);
        
        // for male old shell mating only (AEP ERROR IN OLD CODE HAS >=)
        if (debug) cout<<" to efspbio_matetime "<<endl;
        efspbio_matetime(yr) = modSpBioMateTime_xy(FEMALE,yr);
        if (modSpNumMateTime_xsy(MALE,OLD_SHELL,yr) < modSpNumMateTime_xy(FEMALE,yr)/mate_ratio) {
            efspbio_matetime(yr) = modSpBioMateTime_xy(FEMALE,yr)*((modSpNumMateTime_xsy(MALE,OLD_SHELL,yr)*mate_ratio)/modSpNumMateTime_xy(FEMALE,yr));
        }
        
        // effective sp numbers for new shell females
        if (debug) cout<<" to efspbio_new_matetime "<<endl;
        efspbio_new_matetime(yr) = fspbio_new_matetime(yr);
        if (modSpNumMateTime_xsy(MALE,OLD_SHELL,yr) < modSpNumMateTime_xsy(FEMALE,NEW_SHELL,yr)/mate_ratio) {
            efspbio_new_matetime(yr) = fspbio_new_matetime(yr)*((modSpNumMateTime_xsy(MALE,OLD_SHELL,yr)*mate_ratio)/modSpNumMateTime_xsy(FEMALE,NEW_SHELL,yr));
        }
    }//year loop
    
    //spawning biomass prior to fisheries CAN be evaluated in endyr
    if (debug) cout<<" to sp bio fishtime "<<endl;
    modSpBioFishTime_xy.initialize();
    for (int yr=styr;yr<=endyr;yr++) {
        for (int x=1;x<=nSXs;x++){
            for (int s=1;s<=nSCs;s++) modSpBioFishTime_xy(x,yr) += modFT_PopNum_yxmsz(yr,x,MATURE,s)*wt_xmz(x,MATURE);//dot product
        }//x
    }//yr
    if (debug) cout<<" sex ratio "<<endl;
    // Sex ratio
    for (int yr=styr;yr<=endyr;yr++){
        if ((sum(modNum_xyz(FEMALE,yr))+sum(modNum_xyz(MALE,yr)))<0.01) { 
            modPopXR_y(yr)=0.0;
        } else {
            modPopXR_y(yr)=sum(modNum_xyz(FEMALE,yr))/(sum(modNum_xyz(FEMALE,yr))+sum(modNum_xyz(  MALE,yr)));
        }
    }//yr
    
    // Age-structure
    if (debug) cout<<" to age - struct "<<endl;
    modNumAtAZ_xyaz.initialize();
    
    // initialize
    dvariable tmpi = 1.0;
    for(int j=1;j<=(nages-3);j++) tmpi += mfexp(-j*M_msx(MATURE,OLD_SHELL,FEMALE));
    modNumAtAZ_xyaz(FEMALE,styr,1) = modNum_yxmsz(styr,FEMALE,  MATURE,OLD_SHELL)/(tmpi+(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,FEMALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,FEMALE)))));
    for(int j=1;j<=(nages-2);j++) modNumAtAZ_xyaz(FEMALE,styr,j+1) = modNumAtAZ_xyaz(FEMALE,styr,1)*mfexp(-j*M_msx(MATURE,OLD_SHELL,FEMALE));
    modNumAtAZ_xyaz(FEMALE,styr,nages) = modNumAtAZ_xyaz(FEMALE,styr,1)*(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,FEMALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,FEMALE))));
    if (debug) cout<<" 1 "<<endl;
    
    tmpi = 1.0;
    for(int j=1;j<=(nages-3);j++) tmpi += mfexp(-j*M_msx(MATURE,OLD_SHELL,  MALE));
    modNumAtAZ_xyaz(MALE,styr,1) = modNum_yxmsz(styr,MALE,MATURE,OLD_SHELL)/(tmpi+(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,  MALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,  MALE)))));
    for(int j=1;j<=(nages-2);j++) modNumAtAZ_xyaz(MALE,styr,j+1) = modNumAtAZ_xyaz(MALE,styr,1)*mfexp(-j*M_msx(MATURE,OLD_SHELL,  MALE));
    modNumAtAZ_xyaz(MALE,styr,nages) = modNumAtAZ_xyaz(MALE,styr,1)*(mfexp(-(nages-2)*M_msx(MATURE,OLD_SHELL,  MALE))/(1-mfexp(-M_msx(MATURE,OLD_SHELL,  MALE))));
    if (debug) cout<<" 2 "<<endl;
    
    //numbers at length from styr to endyr
    for (int x=1;x<=nSXs;x++) {
        for (int yr=styr;yr<endyr;yr++) {
            // for numbers by length and age assumes no molting after maturity
            modNumAtAZ_xyaz(x,yr+1,1) = mfexp(-(1-mdptFshs_y(yr))*M_msx(MATURE,NEW_SHELL,x)) * 
                                       elem_prod(S_xsyz(x,NEW_SHELL,yr),mfexp(-mdptFshs_y(yr)*M_msx(MATURE,NEW_SHELL,x))*modNum_yxmsz(yr,x,MATURE,NEW_SHELL));
            for(int j=1;j<=(nages-1);j++){
                modNumAtAZ_xyaz(x,yr+1,j+1) = (mfexp(-(1-mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,x)) * 
                                       elem_prod(S_xsyz(x,OLD_SHELL,yr),mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,x))*modNumAtAZ_xyaz(x,yr,j)));
            }
            modNumAtAZ_xyaz(x,yr+1,nages) += (mfexp(-(1-mdptFshs_y(yr))*M_msx(MATURE,OLD_SHELL,x)) * 
                                   elem_prod(S_xsyz(x,OLD_SHELL,yr),mfexp(-mdptFshs_y(yr)*M_msx(MATURE,OLD_SHELL,x))*modNumAtAZ_xyaz(x,yr,nages)));
            if (debug) cout<<" 3 "<<x<<tb<<yr<<endl;
        }//yr
    }//x
    
    // Legal males at time of survey
    if (debug) cout<<" to legal males "<<endl;
    modSrvNumLegal_sy.initialize();
    modSrvNumLegal_y.initialize();
    modSrvBioLegal_y.initialize();
    modPopNumLegal_y.initialize();
    modPopBioLegal_y.initialize();
    modSrvNum_xyz.initialize();
    modSrvBio_xy.initialize();
    for (int yr=styr;yr<=endyr;yr++) {
        // Selection pattern
        if (yr<1982)            useSelSrv = selSrv1_xz;
        if (yr>1981 && yr<1988) useSelSrv = selSrv2_xz;
        if (yr>1987)            useSelSrv = selSrv3_xz;
        
        // legal male size is based on zLegal
        modSrvNumLegal_sy(NEW_SHELL,yr) = natlength_new(MALE,yr)(iZLegal,nZBs)*useSelSrv(MALE)(iZLegal,nZBs);
        modSrvNumLegal_sy(OLD_SHELL,yr) = natlength_old(MALE,yr)(iZLegal,nZBs)*useSelSrv(MALE)(iZLegal,nZBs);
        modSrvNumLegal_y(yr) = modNum_xyz(MALE,yr)(iZLegal,nZBs)*useSelSrv(MALE)(iZLegal,nZBs);
        modSrvBioLegal_y(yr) = elem_prod(modNum_xyz(MALE,yr),useSelSrv(MALE))(iZLegal,nZBs)*wt_xmz(MALE,MATURE)(iZLegal,nZBs);
        modPopNumLegal_y(yr) = sum(modNum_xyz(MALE,yr)(iZLegal,nZBs));
        modPopBioLegal_y(yr) = modNum_xyz(MALE,yr)(iZLegal,nZBs)*wt_xmz(MALE,MATURE)(iZLegal,nZBs);
        
         for (int x=1;x<=nSXs;x++){
            // survey numbers
            modSrvImmNum_xsy(x,NEW_SHELL,yr) = multQ*modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL)*useSelSrv(x);
            modSrvImmNum_xsy(x,OLD_SHELL,yr) = multQ*modNum_yxmsz(yr,x,IMMATURE,OLD_SHELL)*useSelSrv(x);
            modSrvMatNum_xsy(x,NEW_SHELL,yr) = multQ*modNum_yxmsz(yr,x,  MATURE,NEW_SHELL)*useSelSrv(x);
            modSrvMatNum_xsy(x,OLD_SHELL,yr) = multQ*modNum_yxmsz(yr,x,  MATURE,OLD_SHELL)*useSelSrv(x);
        
            // total survey summaries
            if(x==FEMALE) {
                modSrvBio_xy(x,yr) = multQ*((modNum_yxmsz(yr,x,IMMATURE,NEW_SHELL)                                      )*elem_prod(useSelSrv(x),wt_xmz(FEMALE)(IMMATURE))+
                                            (modNum_yxmsz(yr,x,  MATURE,NEW_SHELL)+modNum_yxmsz(yr,x,  MATURE,OLD_SHELL))*elem_prod(useSelSrv(x),wt_xmz(FEMALE)(  MATURE)));
            } else {
                modSrvBio_xy(x,yr) = multQ*(modNum_xyz(x,yr)*elem_prod(useSelSrv(x),wt_xmz(MALE,  MATURE)));
            }
            modSrvNum_xyz(x,yr) = multQ*elem_prod(modNum_xyz(x,yr),useSelSrv(x));
            modPopNum_y(yr)    += sum(modNum_xyz(x,yr));
        }//x
    }//yr
    
    
    if (debug) cout<<"before sd"<<endl;
    if(sd_phase()){
        sdrDepletion = modPopBio_y(endyr) / modPopBio_y(styr);
        //cout<<"2"<<endl;
        sdrMMB      = modSpBioMateTime_xy( MALE)(styr,endyr-1);
        sdrMFB      = modSpBioMateTime_xy(FEMALE)(styr,endyr-1);
        sdrLnRec    = log(rec_y);
        sdrPrM2M_F = modPrM2M(FEMALE);
        sdrPrM2M_M = modPrM2M(  MALE);
        sdrMnGrw_F = meanPostMoltSize(FEMALE);
        sdrMnGrw_M = meanPostMoltSize(  MALE);
        sdrNatMort_INF(styr,endyr-1) = M_msx(IMMATURE,NEW_SHELL,FEMALE);
        sdrNatMort_INM(styr,endyr-1) = M_msx(IMMATURE,NEW_SHELL,  MALE);
        sdrNatMort_MNF(styr,endyr-1) = M_msx(  MATURE,NEW_SHELL,FEMALE);
        sdrNatMort_MNM(styr,endyr-1) = M_msx(  MATURE,NEW_SHELL,  MALE);
        sdrNatMort_MOF(styr,endyr-1) = M_msx(  MATURE,OLD_SHELL,FEMALE);
        sdrNatMort_MOM(styr,endyr-1) = M_msx(  MATURE,OLD_SHELL,  MALE);
        if (mort_switch==1){
            sdrNatMort_MNF(lyr_mort,uyr_mort) *= pMfac_Big(FEMALE);
            sdrNatMort_MNM(lyr_mort,uyr_mort) *= pMfac_Big(  MALE);
            sdrNatMort_MOF(lyr_mort,uyr_mort) *= pMfac_Big(FEMALE);
            sdrNatMort_MOM(lyr_mort,uyr_mort) *= pMfac_Big(  MALE);
        }
    }
    if (debug) cout<<"done MiscOutput"<<endl;
    //exit(-1);

  
// ==========================================================================
FUNCTION void writeMyProjectionFile(ofstream& os)
//    cout<<"starting writeMyProjectionFile"<<endl;
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
      os<<mean(rec_y(1982,endyr))*1000<<tb<<tb<<"#total average recruitment for BXX/Bmsy calculation (1000's of recruits)"<<endl;
      os<<modSpBioMateTime_xy(   MALE,endyr-1)<<tb<<tb<<"#'current' spawning biomass (MMB, 1000's t) "<<endl;
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
      os << 0.5*rec_y*1000 <<endl;//females; was endyr-1
      os << 0.5*rec_y*1000 <<endl;//  males; was endyr-1     
      os<<"#asmtModSpB(mMnYr,mMxYr-1): male spawning biomass at matetime (1000's t) for str year to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
      os<<modSpBioMateTime_xy(   MALE)(styr,endyr-1)<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#Pop info in last year of assessment model"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#numbers at size immature new shell female, male in final year (1000's)"<<endl;
      os<<modNum_yxmsz(endyr,FEMALE,IMMATURE,NEW_SHELL)*1000<<endl;
      os<<modNum_yxmsz(endyr,  MALE,IMMATURE,NEW_SHELL)*1000<<endl;
      os<<"#numbers at length immature old shell female male last year (1000's)"<<endl;
      os<<modNum_yxmsz(endyr,FEMALE,IMMATURE,OLD_SHELL)*1000<<endl;
      os<<modNum_yxmsz(endyr,  MALE,IMMATURE,OLD_SHELL)*1000<<endl;
      os<<"#numbers at length mature new shell female male last year (1000's)"<<endl;
      os<<modNum_yxmsz(endyr,FEMALE,  MATURE,NEW_SHELL)*1000<<endl;
      os<<modNum_yxmsz(endyr,  MALE,  MATURE,NEW_SHELL)*1000<<endl;
      os<<"#numbers at length mature old shell female male last year (1000's) "<<endl;
      os<<modNum_yxmsz(endyr,FEMALE,  MATURE,OLD_SHELL)*1000<<endl;
      os<<modNum_yxmsz(endyr,  MALE,  MATURE,OLD_SHELL)*1000<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#Fisheries info"<<endl;
      os<<"#---------------------------"<<endl;
      os<<"#time of catch in fraction of year from survey - 7 months"<<endl;
      os<<mdptFshs_y(endyr-1)<<endl;//IMPORTANT CHANGE; was endyr
      
      os<<0<<tb<<tb<<"#inpFmTCF: input F for directed Tanner crab fishing mortality"<<endl;
      os<<mean(fSCF_xy(MALE)(endyr-5,endyr-1))<<tb<<tb<<"#inpFmSCF: input male F for snow crab fishing mortality"<<endl;
      os<<mean(fRKF_xy(MALE)(endyr-5,endyr-1))<<tb<<tb<<"#inpFmRKF: input male F for BBRKC  fishing mortality"<<endl;
      os<<mean(fGTF_xy(MALE)(endyr-5,endyr-1))<<tb<<tb<<"#inpFmGTF: input male F for groundfish fishery fishing mortality"<<endl;
      
      os<<pAvgLnF_TCFF<<tb<<tb<<"#pAvgLnF_TCFF: ln-scale offset to F for female bycatch in the directed Tanner crab fishery"<<endl;
      os<<pAvgLnF_SCFF<<tb<<tb<<"#pAvgLnF_SCFF: ln-scale offset to F for female bycatch in the snow crab fishery"<<endl;
      os<<pAvgLnF_RKFF<<tb<<tb<<"#pAvgLnF_RKFF: ln-scale offset to F for female bycatch in the BBRKC fishery"<<endl;
      os<<pAvgLnF_GTFF<<tb<<tb<<"#pAvgLnF_GTFF: ln-scale offset to F for female bycatch in the groundfish fishery"<<endl;
      
      os<<"#selTCF_TotMale(nSCs,nSXs): average of last 4 years selTCFM_syz total male new old shell"<<endl;
      os<<(selTCFM_syz(NEW_SHELL,endyr-4)+selTCFM_syz(NEW_SHELL,endyr-3)+selTCFM_syz(NEW_SHELL,endyr-2)+selTCFM_syz(NEW_SHELL,endyr-1))/4.0<<endl;
      os<<(selTCFM_syz(OLD_SHELL,endyr-4)+selTCFM_syz(OLD_SHELL,endyr-3)+selTCFM_syz(OLD_SHELL,endyr-2)+selTCFM_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      if (optFM==0){
            os<<"#selTCF_RetMale(nSCs): average of last 4 years selTCFM_syz retained curve male new old shell"<<endl;
            os<<(selTCFR_syz(NEW_SHELL,endyr-4)+selTCFR_syz(NEW_SHELL,endyr-3)+selTCFR_syz(NEW_SHELL,endyr-2)+selTCFR_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
            os<<(selTCFR_syz(OLD_SHELL,endyr-4)+selTCFR_syz(OLD_SHELL,endyr-3)+selTCFR_syz(OLD_SHELL,endyr-2)+selTCFR_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      } else {
            os<<"#retFcn_syz(nSCs): average of last four years"<<endl;
            os<<(retFcn_syz(NEW_SHELL,endyr-4)+retFcn_syz(NEW_SHELL,endyr-3)+retFcn_syz(NEW_SHELL,endyr-2)+retFcn_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
            os<<(retFcn_syz(OLD_SHELL,endyr-4)+retFcn_syz(OLD_SHELL,endyr-3)+retFcn_syz(OLD_SHELL,endyr-2)+retFcn_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      }
      os<<"#selTCF_TotMaleEast(nSCs,nSXs): set same as average total"<<endl;
      os<<(selTCFM_syz(NEW_SHELL,endyr-4)+selTCFM_syz(NEW_SHELL,endyr-3)+selTCFM_syz(NEW_SHELL,endyr-2)+selTCFM_syz(NEW_SHELL,endyr-1))/4.0<<endl;
      os<<(selTCFM_syz(OLD_SHELL,endyr-4)+selTCFM_syz(OLD_SHELL,endyr-3)+selTCFM_syz(OLD_SHELL,endyr-2)+selTCFM_syz(OLD_SHELL,endyr-1))/4.0<<endl;
//      os<<"#selTCF_RetMaleEast(nSCs,nSXs): set same as avg retained"<<endl;
//      os<<(selTCFR_syz(NEW_SHELL,endyr-4)+selTCFR_syz(NEW_SHELL,endyr-3)+selTCFR_syz(NEW_SHELL,endyr-2)+selTCFR_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
//      os<<(selTCFR_syz(OLD_SHELL,endyr-4)+selTCFR_syz(OLD_SHELL,endyr-3)+selTCFR_syz(OLD_SHELL,endyr-2)+selTCFR_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      if (optFM==0){
            os<<"#selTCF_RetMaleEast(nSCs,nSXs): set same as avg retained"<<endl;
            os<<(selTCFR_syz(NEW_SHELL,endyr-4)+selTCFR_syz(NEW_SHELL,endyr-3)+selTCFR_syz(NEW_SHELL,endyr-2)+selTCFR_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
            os<<(selTCFR_syz(OLD_SHELL,endyr-4)+selTCFR_syz(OLD_SHELL,endyr-3)+selTCFR_syz(OLD_SHELL,endyr-2)+selTCFR_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      } else {
            os<<"#retFcn_syz(nSCs) for East: same as average retained"<<endl;
            os<<(retFcn_syz(NEW_SHELL,endyr-4)+retFcn_syz(NEW_SHELL,endyr-3)+retFcn_syz(NEW_SHELL,endyr-2)+retFcn_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
            os<<(retFcn_syz(OLD_SHELL,endyr-4)+retFcn_syz(OLD_SHELL,endyr-3)+retFcn_syz(OLD_SHELL,endyr-2)+retFcn_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      }
      os<<"#selTCF_TotMaleWest(nSCs,nSXs): set same as average total"<<endl;
      os<<(selTCFM_syz(NEW_SHELL,endyr-4)+selTCFM_syz(NEW_SHELL,endyr-3)+selTCFM_syz(NEW_SHELL,endyr-2)+selTCFM_syz(NEW_SHELL,endyr-1))/4.0<<endl;
      os<<(selTCFM_syz(OLD_SHELL,endyr-4)+selTCFM_syz(OLD_SHELL,endyr-3)+selTCFM_syz(OLD_SHELL,endyr-2)+selTCFM_syz(OLD_SHELL,endyr-1))/4.0<<endl;
//      os<<"#selTCF_RetMaleWest(nSCs,nSXs): SET SAME AS AVG RETAINED, BUT SHIFTED TO LOWER END BY 10 mm"<<endl;
//      os<<(selTCFR_syz(NEW_SHELL,endyr-4)+selTCFR_syz(NEW_SHELL,endyr-3)+selTCFR_syz(NEW_SHELL,endyr-2)+selTCFR_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
//      os<<(selTCFR_syz(OLD_SHELL,endyr-4)+selTCFR_syz(OLD_SHELL,endyr-3)+selTCFR_syz(OLD_SHELL,endyr-2)+selTCFR_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      if (optFM==0){
            os<<"#selTCF_RetMaleWest(nSCs,nSXs): SET SAME AS AVG RETAINED, BUT SHIFTED TO LOWER END BY 10 mm"<<endl;
            os<<(selTCFR_syz(NEW_SHELL,endyr-4)+selTCFR_syz(NEW_SHELL,endyr-3)+selTCFR_syz(NEW_SHELL,endyr-2)+selTCFR_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
            os<<(selTCFR_syz(OLD_SHELL,endyr-4)+selTCFR_syz(OLD_SHELL,endyr-3)+selTCFR_syz(OLD_SHELL,endyr-2)+selTCFR_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      } else {
            os<<"#retFcn_syz(nSCs) for West: same as average retained"<<endl;
            os<<(retFcn_syz(NEW_SHELL,endyr-4)+retFcn_syz(NEW_SHELL,endyr-3)+retFcn_syz(NEW_SHELL,endyr-2)+retFcn_syz(NEW_SHELL,endyr-1))/4.0<<endl;//IMPORTANT CHANGE: was only over last 3 years (but said 4)
            os<<(retFcn_syz(OLD_SHELL,endyr-4)+retFcn_syz(OLD_SHELL,endyr-3)+retFcn_syz(OLD_SHELL,endyr-2)+retFcn_syz(OLD_SHELL,endyr-1))/4.0<<endl;
      }
      os<<"#selTCF_Female(nZs): selectivity for females in directed fishery"<<endl;
      os<<selTCFF_z<<endl;
      os<<"#selSCF_cxz(nSXs,nZs): selectivity in snow crab fishery"<<endl;
      os<<selSCF_cxz(3,FEMALE)<<endl;
      os<<selSCF_cxz(3,  MALE)<<endl;
      os<<"#selRKF_cxz(nSXs,nZs): selectivity in BBRKC fishery"<<endl;
      os<<selRKF_cxz(3,FEMALE)<<endl;
      os<<selRKF_cxz(3,  MALE)<<endl;      
      os<<"#selGTF_cxz(nSXs,nZs): selectivity in groundfish fishery"<<endl;
      os<<selGTF_cxz(3)<<endl;
      
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
      os<<wt_xmz(FEMALE)(IMMATURE)*0.001<<endl;
      os<<"#weight at length female mature (t)"<<endl;
      os<<wt_xmz(FEMALE)(MATURE)*0.001<<endl;
      os<<"#weight at length male (t)"<<endl;
      os<<wt_xmz(MALE,  MATURE)*0.001<<endl;
      os<<"#tmZtoZ_xzz: size transition matrix"<<endl;
      os<<prGr_xzz<<endl;      
      os<<"#prMatNS(nSXs,nZs): maturity curve new shell female male"<<endl;
      os<<modPrM2M(FEMALE)<<endl;
      os<<modPrM2M(MALE)<<endl;
      os<<"#prMoltImm(nSXs,nZs): molting probability immature female male"<<endl;
      os<<prMoltImm_xz<<endl;
      os<<"#prMoltMat(nSXs,nZs): molting probability mature female male"<<endl;
      os<<prMoltMat_xz<<endl;
      os<<0.5     <<tb<<tb<<"#recPropAsMale: proportion recruiting as males"<<endl;
      os<<1.0     <<tb<<tb<<"#recPropAsNewShell: prop recruits to new shell"<<endl;
      os<<"#recPropAtZ(nZs): distribution of recruits to length bins"<<endl;
      os<<prRec_z<<endl;
      os<<"#propEast(nZs): proportion of population at size east of 166W"<<endl;
      os<<"???????????????????????????????"<<endl;
  
// ==========================================================================
//"Old style" R file
FUNCTION void writeToR_OLD(ofstream& R_out)
//    cout<<"starting writeToR_OLD"<<endl;

        int ii;
        dvar_vector preds_sexr(styr,endyr);
        dvar_matrix tmpo(1,2,styr,endyr);
        dvar_matrix tmpp(1,2,styr,endyr);
        dvar_vector obs_tmp(styr,endyr);
        dvariable ghl,ghl_number;
        dvariable hrate;
        dvar_vector totcatch(styr, endyr-1);  
        dvar_vector tmpp1(1,nZBs);
        dvar_vector tmpp2(1,nZBs);
        dvar_vector tmpp3(1,nZBs);
        dvar_vector tmpp4(1,nZBs);
        
        Misc_output();
        
        REP2R2(mod.styr,styr);
        REP2R2(mod.endyr,endyr);
        REP2R2(mod.obsyr,ptrMDS->pTSD->yrsAbund[1]);
        REP2R2(mod.pltyr,1969);
        REP2R2(mod.zBs,zBs);        
        REP2R2(mod.LZ,zLegal);
        
        REP2R2(mod.recLag,recLag);
        REP2R2(mod.mnYrRecDevsHist,mnYrRecDevsHist);
        REP2R2(mod.mnYrRecCurr,mnYrRecCurr);
        
        //model options
        REP2R2(mod.optFM,optFM);
        
        //model processes
        //--molting
        REP2R2(pop.prMolt.IF,prMoltImm_xz(FEMALE));
        REP2R2(pop.prMolt.IM,prMoltImm_xz(  MALE));
        REP2R2(pop.prMolt.MM,prMoltMat_xz(  MALE));
        //--terminal molt
        REP2R2(pop.prM2M.F,modPrM2M(FEMALE));
        REP2R2(pop.prM2M.M,modPrM2M(  MALE));
        //--growth
        REP2R2(pop.grw.AF1,pGrAF1);
        REP2R2(pop.grw.BF1,pGrBF1);
        REP2R2(pop.grw.mnPMZ.F,meanPostMoltSize(FEMALE));
        REP2R2(pop.grw.prGr_xzz.F,prGr_xzz(FEMALE));
        REP2R2(pop.grw.AM1,pGrAM1);
        REP2R2(pop.grw.BM1,pGrBM1);
        REP2R2(pop.grw.mnPMZ.M,meanPostMoltSize(MALE)); 
        REP2R2(pop.grw.prGr_xzz.M,prGr_xzz(  MALE));
        //--natural mortality
        dvar_vector mrt(styr,endyr-1);
        mrt = M_msx(IMMATURE,NEW_SHELL,FEMALE);
        REP2R2(pop.M.INF,mrt);
        mrt = M_msx(IMMATURE,NEW_SHELL,  MALE);
        REP2R2(pop.M.INM,mrt);
        mrt = M_msx(MATURE,NEW_SHELL,FEMALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= pMfac_Big(FEMALE);
        REP2R2(pop.M.MNF,mrt);
        mrt = M_msx(MATURE,OLD_SHELL,FEMALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= pMfac_Big(FEMALE);
        REP2R2(pop.M.MOF,mrt);
        mrt = M_msx(MATURE,NEW_SHELL,  MALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= pMfac_Big(  MALE);
        REP2R2(pop.M.MNM,mrt);
        mrt = M_msx(MATURE,OLD_SHELL,  MALE); if (mort_switch==1) mrt(lyr_mort,uyr_mort) *= pMfac_Big(  MALE);
        REP2R2(pop.M.MOM,mrt);
        //--recruitment
        REP2R2(pop.R,rec_y);    //total recruitment (millions))
        REP2R2(pop.prR_z,prRec_z);//recruitment size distribution
        
        
        //population quantities
        //--abundance (millions) and biomass (1000's t))
        REP2R2(pop.num,modPopNum_y);
        REP2R2(pop.bio,modPopBio_y);
        //--spawning abundance and biomass 
        REP2R2(pop.spnum.MF,modSpNumMateTime_xy(FEMALE));
        REP2R2(pop.spnum.MM,modSpNumMateTime_xy(  MALE));
        REP2R2(pop.spnum.MNF,modSpNumMateTime_xsy(FEMALE,NEW_SHELL));
        REP2R2(pop.spnum.MOM,modSpNumMateTime_xsy(  MALE,OLD_SHELL));
        REP2R2(pop.bio.MF,fspbio);
        REP2R2(pop.bio.MM,mspbio);
        REP2R2(pop.MFB,modSpBioMateTime_xy(FEMALE));
        REP2R2(pop.MMB,modSpBioMateTime_xy(  MALE));
        REP2R2(pop.MFB.NS,fspbio_new_matetime);
        REP2R2(pop.MMB.OS,mspbio_old_matetime);
        REP2R2(pop.effMFB,efspbio_matetime);
        REP2R2(pop.effMMB,emspbio_matetime);
        REP2R2(pop.effMFB.NS,efspbio_new_matetime);
        //--legal male numbers and biomass
        REP2R2(pop.num.LMs,modPopNumLegal_y);
        REP2R2(pop.bio.LMs,modPopBioLegal_y);
        REP2R2(pop.bio.LMs.FT,modFT_PopBioLegal_y);
        //--numbers-at-size
        ivector prm5(1,5); prm5[1]=4; prm5[2]=1; prm5[3]=2; prm5[4]=3; prm5[5]=5;
        myWriteN_yxmsz(R_out,modNum_yxmsz);
        
//        d5_array modNum_xmsyz  = wts::permuteDims(prm5, dmodNum_yxmsz);
//        REP2R2(pop.NatZ.INF,modNum_xmsyz(FEMALE,IMMATURE,NEW_SHELL));
//        REP2R2(pop.NatZ.IOF,modNum_xmsyz(FEMALE,IMMATURE,OLD_SHELL));
//        REP2R2(pop.NatZ.MNF,modNum_xmsyz(FEMALE,  MATURE,NEW_SHELL));
//        REP2R2(pop.NatZ.MOF,modNum_xmsyz(FEMALE,  MATURE,OLD_SHELL));
//        REP2R2(pop.NatZ.INM,modNum_xmsyz(  MALE,IMMATURE,NEW_SHELL));
//        REP2R2(pop.NatZ.IOM,modNum_xmsyz(  MALE,IMMATURE,OLD_SHELL));
//        REP2R2(pop.NatZ.MNM,modNum_xmsyz(  MALE,  MATURE,NEW_SHELL));
//        REP2R2(pop.NatZ.MOM,modNum_xmsyz(  MALE,  MATURE,OLD_SHELL));
//        REP2R2(pop.NatZ.F,modNum_xyz(FEMALE));
//        REP2R2(pop.NatZ.M,modNum_xyz(  MALE));
        
        //survey quantities
        //--numbers-at-size (millions), biomass-at-size (1000's t)
        {
            double totSrvNum;
            d5_array B_xmsnz(1,nSXs,1,nMSs,1,nSCs,1,ptrMDS->pTSD->nyNatZ,1,nZBs);
            d5_array n_xmsyz(1,nSXs,1,nMSs,1,nSCs,styr,endyr,1,nZBs);
            d5_array b_xmsyz(1,nSXs,1,nMSs,1,nSCs,styr,endyr,1,nZBs);
            B_xmsnz.initialize();
            n_xmsyz.initialize();
            b_xmsyz.initialize();
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                        for (int n=1;n<=ptrMDS->pTSD->nyNatZ;n++){
                            B_xmsnz(x,m,s,n) = elem_prod(wt_xmz(x,m),obsSrvNatZs_msxnz(m,s,x,n));
                        }
                        for (int y=styr;y<=endyr;y++){
                            totSrvNum = value(modSrvNum_xy(FEMALE,y) + modSrvNum_xy(MALE,y));
                            n_xmsyz(x,m,s,y) = totSrvNum*value(modSrvPrNatZ_msxyz(m,s,x,y));
                            b_xmsyz(x,m,s,y) = elem_prod(wt_xmz(x,m),n_xmsyz(x,m,s,y));
                        }
                    }
                }
            }
            //--observed abundance (millions)
            REP2R2(srv.obs.num.yrs,ptrMDS->pTSD->yrsAbund);
            REP2R2(srv.obs.NatZ.INF,obsSrvNatZs_msxnz(IMMATURE,NEW_SHELL,FEMALE));
            REP2R2(srv.obs.NatZ.IOF,obsSrvNatZs_msxnz(IMMATURE,OLD_SHELL,FEMALE));
            REP2R2(srv.obs.NatZ.MNF,obsSrvNatZs_msxnz(  MATURE,NEW_SHELL,FEMALE));
            REP2R2(srv.obs.NatZ.MOF,obsSrvNatZs_msxnz(  MATURE,OLD_SHELL,FEMALE));
            REP2R2(srv.obs.NatZ.INM,obsSrvNatZs_msxnz(IMMATURE,NEW_SHELL,  MALE));
            REP2R2(srv.obs.NatZ.IOM,obsSrvNatZs_msxnz(IMMATURE,OLD_SHELL,  MALE));
            REP2R2(srv.obs.NatZ.MNM,obsSrvNatZs_msxnz(  MATURE,NEW_SHELL,  MALE));
            REP2R2(srv.obs.NatZ.MOM,obsSrvNatZs_msxnz(  MATURE,OLD_SHELL,  MALE));        
            //predicted abundance
            REP2R2(srv.mod.NatZ.INF,n_xmsyz(FEMALE,IMMATURE,NEW_SHELL));
            REP2R2(srv.mod.NatZ.MNF,n_xmsyz(FEMALE,  MATURE,NEW_SHELL));
            REP2R2(srv.mod.NatZ.MOF,n_xmsyz(FEMALE,  MATURE,OLD_SHELL));
            REP2R2(srv.mod.NatZ.INM,n_xmsyz(  MALE,IMMATURE,NEW_SHELL));
            REP2R2(srv.mod.NatZ.MNM,n_xmsyz(  MALE,  MATURE,NEW_SHELL));
            REP2R2(srv.mod.NatZ.MOM,n_xmsyz(  MALE,  MATURE,OLD_SHELL));
            //--observed biomass-at-size (1000's t)
            REP2R2(srv.obs.BatZ.INF,B_xmsnz(FEMALE,IMMATURE,NEW_SHELL));
            REP2R2(srv.obs.BatZ.MNF,B_xmsnz(FEMALE,  MATURE,NEW_SHELL));
            REP2R2(srv.obs.BatZ.MOF,B_xmsnz(FEMALE,  MATURE,OLD_SHELL));
            REP2R2(srv.obs.BatZ.INM,B_xmsnz(  MALE,IMMATURE,NEW_SHELL));
            REP2R2(srv.obs.BatZ.MNM,B_xmsnz(  MALE,  MATURE,NEW_SHELL));
            REP2R2(srv.obs.BatZ.MOM,B_xmsnz(  MALE,  MATURE,OLD_SHELL));
            //--predicted biomass-at-size (1000's t)
            REP2R2(srv.mod.BatZ.INF,b_xmsyz(FEMALE,IMMATURE,NEW_SHELL));
            REP2R2(srv.mod.BatZ.MNF,b_xmsyz(FEMALE,  MATURE,NEW_SHELL));
            REP2R2(srv.mod.BatZ.MOF,b_xmsyz(FEMALE,  MATURE,OLD_SHELL));
            REP2R2(srv.mod.BatZ.INM,b_xmsyz(  MALE,IMMATURE,NEW_SHELL));
            REP2R2(srv.mod.BatZ.MNM,b_xmsyz(  MALE,  MATURE,NEW_SHELL));
            REP2R2(srv.mod.BatZ.MOM,b_xmsyz(  MALE,  MATURE,OLD_SHELL));
        }
        //--observed biomass (1000's t)
        REP2R2(srv.obs.bio.yrs,ptrMDS->pTSD->yrsAbund);
        REP2R2(srv.obs.bio.F,obsSrvBio_xy(FEMALE));
        REP2R2(srv.obs.bio.M,obsSrvBio_xy(  MALE));
        REP2R2(srv.obs.bio.MF,obsSrvMatBio_xy(FEMALE));
        REP2R2(srv.obs.bio.MM,obsSrvMatBio_xy(  MALE));
        REP2R2(srv.obs.bio.cv.MF,obsSrvCV_xn(FEMALE));
        REP2R2(srv.obs.bio.cv.MM,obsSrvCV_xn(  MALE));
        //--predicted biomass (1000's t)
        REP2R2(srv.mod.bio.F,modSrvBio_xy(FEMALE));
        REP2R2(srv.mod.bio.M,modSrvBio_xy(  MALE));
        REP2R2(srv.mod.bio.MF,modSrvMatBio_xy(FEMALE));
        REP2R2(srv.mod.bio.MM,modSrvMatBio_xy(  MALE));
        
        //--proportions-at-size
        REP2R2(srv.obs.PrNatZ.yrs,yrsObsZCsSrv_n);
        REP2R2(srv.obs.PrNatZ.INF,obsSrvPrNatZ_msxnz(IMMATURE,NEW_SHELL,FEMALE));
        REP2R2(srv.obs.PrNatZ.IOF,obsSrvPrNatZ_msxnz(IMMATURE,OLD_SHELL,FEMALE));
        REP2R2(srv.obs.PrNatZ.MNF,obsSrvPrNatZ_msxnz(  MATURE,NEW_SHELL,FEMALE));
        REP2R2(srv.obs.PrNatZ.MOF,obsSrvPrNatZ_msxnz(  MATURE,OLD_SHELL,FEMALE));
        REP2RTS(srv.mod.PrNatZ.INF,modSrvPrNatZ_msxyz(IMMATURE,NEW_SHELL,FEMALE),yrsObsZCsSrv_n);
        REP2RTS(srv.mod.PrNatZ.IOF,modSrvPrNatZ_msxyz(IMMATURE,OLD_SHELL,FEMALE),yrsObsZCsSrv_n);
        REP2RTS(srv.mod.PrNatZ.MNF,modSrvPrNatZ_msxyz(  MATURE,NEW_SHELL,FEMALE),yrsObsZCsSrv_n);
        REP2RTS(srv.mod.PrNatZ.MOF,modSrvPrNatZ_msxyz(  MATURE,OLD_SHELL,FEMALE),yrsObsZCsSrv_n);
        
        REP2R2(srv.obs.PrNatZ.INM,obsSrvPrNatZ_msxnz(IMMATURE,NEW_SHELL,  MALE));
        REP2R2(srv.obs.PrNatZ.IOM,obsSrvPrNatZ_msxnz(IMMATURE,OLD_SHELL,  MALE));
        REP2R2(srv.obs.PrNatZ.MNM,obsSrvPrNatZ_msxnz(  MATURE,NEW_SHELL,  MALE));
        REP2R2(srv.obs.PrNatZ.MOM,obsSrvPrNatZ_msxnz(  MATURE,OLD_SHELL,  MALE));
        REP2RTS(srv.mod.PrNatZ.INM,modSrvPrNatZ_msxyz(IMMATURE,NEW_SHELL,  MALE),yrsObsZCsSrv_n);
        REP2RTS(srv.mod.PrNatZ.IOM,modSrvPrNatZ_msxyz(IMMATURE,OLD_SHELL,  MALE),yrsObsZCsSrv_n);
        REP2RTS(srv.mod.PrNatZ.MNM,modSrvPrNatZ_msxyz(  MATURE,NEW_SHELL,  MALE),yrsObsZCsSrv_n);
        REP2RTS(srv.mod.PrNatZ.MOM,modSrvPrNatZ_msxyz(  MATURE,OLD_SHELL,  MALE),yrsObsZCsSrv_n);
        
//        //--pearsons residuals for survey size comps
//        dmatrix mod_nz(1,nObsZCsSrv,1,nZBs);
//        dmatrix prs_nz(1,nObsZCsSrv,1,nZBs);
//        mod_nz = value(modSrvPrNatZ_mxyz(IMMATURE,FEMALE)(yrsObsZCsSrv_n));
//        prs_nz = elem_div(obsSrvPrNatZ_mxnz(IMMATURE,FEMALE)-mod_nz,sqrt(elem_prod(mod_nz,1.0-mod_nz)));
//        REP2R2(srv.prs.PrNatZ.IF,prs_nz);
        
        //--legal males in survey (in millions and 1000's t))
        REP2R2(srv.obs.num.LMs,obsSrvNumLegal_n);
        REP2R2(srv.mod.num.LMs,modSrvNumLegal_y);
        REP2R2(srv.obs.bio.LMs,obsSrvBioLegal_n);
        REP2R2(srv.mod.bio.LMs,modSrvBioLegal_y);        
        
        //Fisheries 
        REP2R2(qSCF,qSCF);
        REP2R2(qRKF,qRKF);
        //--TCF
        ivector yrs_ret(1965,endyr-1);
        yrs_ret.fill_seqadd(1965,1);
        REP2R2(fsh.obs.ret.bio.yrs.TCF,yrs_ret);
        REP2R2(fsh.obs.ret.bio.TCF.M,obsRetCatchBio_y);
        REP2R2(fsh.mod.ret.bio.TCF.M,predRetBioMortTCFM_y);
        REP2R2(fsh.mod.ret.bio.TCF.NM,predRetNumMortTCFM_syz(NEW_SHELL)*wt_xmz(MALE,  MATURE));
        REP2R2(fsh.mod.ret.bio.TCF.OM,predRetNumMortTCFM_syz(OLD_SHELL)*wt_xmz(MALE,  MATURE));
        
        REP2R2(fsh.obs.totm.bio.yrs.TCF,yrsObsDscTCF_n);
        REP2R2(fsh.obs.totm.bio.TCF.F,obsDscBioMortTCF_xn(FEMALE));
        REP2R2(fsh.obs.totm.bio.TCF.M,obsTotBioMortTCFM_n);
        REP2R2(fsh.mod.totm.bio.TCF.F,predDscBioMortTCF_xy(FEMALE));
        REP2R2(fsh.mod.totm.bio.TCF.M,predTotBioMortTCFM_y);
        REP2R2(fsh.mod.totm.bio.TCF.NM,(predTotNumMortTCFM_syz(NEW_SHELL)*wt_xmz(MALE,  MATURE)));
        REP2R2(fsh.mod.totm.bio.TCF.OM,(predTotNumMortTCFM_syz(OLD_SHELL)*wt_xmz(MALE,  MATURE)));
        
        REP2R2(fsh.obs.dscm.bio.TCF.M.chk,(obsTotBioMortTCFM_n-obsRetCatchBio_y(yrsObsDscTCF_n)));
        REP2R2(fsh.obs.dscm.bio.TCF.M,obsDscBioMortTCF_xn(  MALE));
        REP2R2(fsh.mod.dscm.bio.TCF.M,predTotBioMortTCFM_y-predRetBioMortTCFM_y);
        
        //----retained catch size comps
        REP2R2(fsh.obs.ret.PrNatZ.yrs.TCF,yrsObsRetZCsTCF_n);
        REP2R2(fsh.obs.ret.PrNatZ.TCF.M, obsPrNatZ_TCFR_sn(NEW_SHELL)+obsPrNatZ_TCFR_sn(OLD_SHELL));
        REP2R2(fsh.obs.ret.PrNatZ.TCF.NM,obsPrNatZ_TCFR_sn(NEW_SHELL));
        REP2R2(fsh.obs.ret.PrNatZ.TCF.OM,obsPrNatZ_TCFR_sn(OLD_SHELL));
        REP2R2(fsh.mod.ret.PrNatZ.TCF.M, modPrNatZ_TCFR_syz(NEW_SHELL)+modPrNatZ_TCFR_syz(OLD_SHELL));
        REP2R2(fsh.mod.ret.PrNatZ.TCF.NM,modPrNatZ_TCFR_syz(NEW_SHELL));
        REP2R2(fsh.mod.ret.PrNatZ.TCF.OM,modPrNatZ_TCFR_syz(OLD_SHELL));
        //----total catch size comps
        //------males
        REP2R2(fsh.obs.tot.PrNatZ.yrs.TCF,yrsObsZCsTCFM_n);
        REP2R2(fsh.obs.tot.PrNatZ.TCF.M, obsPrNatZ_TCFM_snz(NEW_SHELL)+obsPrNatZ_TCFM_snz(OLD_SHELL));
        REP2R2(fsh.obs.tot.PrNatZ.TCF.NM,obsPrNatZ_TCFM_snz(NEW_SHELL));
        REP2R2(fsh.obs.tot.PrNatZ.TCF.OM,obsPrNatZ_TCFM_snz(OLD_SHELL));
        REP2R2(fsh.mod.tot.PrNatZ.TCF.M, modPrNatZ_TCFM_syz(NEW_SHELL)+modPrNatZ_TCFM_syz(OLD_SHELL));
        REP2R2(fsh.mod.tot.PrNatZ.TCF.NM,modPrNatZ_TCFM_syz(NEW_SHELL));
        REP2R2(fsh.mod.tot.PrNatZ.TCF.OM,modPrNatZ_TCFM_syz(OLD_SHELL));
        //------females
        REP2R2(fsh.obs.tot.PrNatZ.TCF.F,obsPrNatZ_TCFF_nz);
        REP2R2(fsh.mod.tot.PrNatZ.TCF.F,modPrNatZ_TCFF_yz);
        //--SCF
        REP2R2(fsh.obs.totm.bio.yrs.SCF,yrsObsDscSCF);
        REP2R2(fsh.obs.totm.bio.SCF.M,obsDscBioMortSCF_xn(  MALE));
        REP2R2(fsh.obs.totm.bio.SCF.F,obsDscBioMortSCF_xn(FEMALE));
        REP2R2(fsh.mod.totm.bio.SCF.M,predDscBioMortSCF_xy(  MALE));
        REP2R2(fsh.mod.totm.bio.SCF.F,predDscBioMortSCF_xy(FEMALE));
        //----total catch size comps 
        REP2R2(fsh.obs.tot.PrNatZ.yrs.SCF,yrsObsZCsSCF_n);
        REP2R2(fsh.obs.tot.PrNatZ.SCF.M,obsPrNatZ_SCF_xnz(  MALE));
        REP2R2(fsh.obs.tot.PrNatZ.SCF.F,obsPrNatZ_SCF_xnz(FEMALE));
        REP2R2(fsh.mod.tot.PrNatZ.SCF.M,modPrNatZ_SCF_xyz(  MALE));
        REP2R2(fsh.mod.tot.PrNatZ.SCF.F,modPrNatZ_SCF_xyz(FEMALE));
        //--RKF
        REP2R2(fsh.obs.totm.bio.yrs.RKF,yrsObsDscRKF);
        REP2R2(fsh.obs.totm.bio.RKF.M,obsDscBioMortRKF_xn(  MALE));
        REP2R2(fsh.obs.totm.bio.RKF.F,obsDscBioMortRKF_xn(FEMALE));
        REP2R2(fsh.mod.totm.bio.RKF.M,predDscBioMortRKF_xy(  MALE));
        REP2R2(fsh.mod.totm.bio.RKF.F,predDscBioMortRKF_xy(FEMALE));
        //----total catch size comps
        REP2R2(fsh.obs.tot.PrNatZ.yrs.RKF,yrsObsZCsRKF_n);
        REP2R2(fsh.obs.tot.PrNatZ.RKF.M,obsPrNatZ_RKF_xnz(  MALE));
        REP2R2(fsh.obs.tot.PrNatZ.RKF.F,obsPrNatZ_RKF_xnz(FEMALE));
        REP2R2(fsh.mod.tot.PrNatZ.RKF.M,modPrNatZ_RKF_xyz(  MALE));
        REP2R2(fsh.mod.tot.PrNatZ.RKF.F,modPrNatZ_RKF_xyz(FEMALE));
        //--GTF
        REP2R2(fsh.obs.totm.bio.yrs.GTF,yrsObsDscGTF);
        REP2R2(fsh.obs.totm.bio.GTF,obsDscBioMortGTF_n);
        REP2R2(fsh.mod.totm.bio.GTF,predDscBioMortGTF_xy(  MALE)+predDscBioMortGTF_xy(FEMALE));
        REP2R2(fsh.mod.totm.bio.GTF.M,predDscBioMortGTF_xy(  MALE));
        REP2R2(fsh.mod.totm.bio.GTF.F,predDscBioMortGTF_xy(FEMALE));
        //----total catch size comps
        REP2R2(fsh.obs.tot.PrNatZ.yrs.GTF,yrsObsZCsGTF_n);
        REP2R2(fsh.obs.tot.PrNatZ.GTF.M,obsPrNatZ_GTF_xnz(  MALE));
        REP2R2(fsh.obs.tot.PrNatZ.GTF.F,obsPrNatZ_GTF_xnz(FEMALE));
        REP2R2(fsh.mod.tot.PrNatZ.GTF.M,modPrNatZ_GTF_xyz(  MALE));
        REP2R2(fsh.mod.tot.PrNatZ.GTF.F,modPrNatZ_GTF_xyz(FEMALE));
        
        //total fishing mortality
        REP2R2(fsh.mod.totm.bio.All.M,predTotBioMortTCFM_y        +predDscBioMortRKF_xy(  MALE)+predDscBioMortSCF_xy(  MALE)+predDscBioMortGTF_xy(  MALE));
        REP2R2(fsh.mod.totm.bio.All.F,predDscBioMortTCF_xy(FEMALE)+predDscBioMortRKF_xy(FEMALE)+predDscBioMortSCF_xy(FEMALE)+predDscBioMortGTF_xy(FEMALE));
        
        adstring ftype;
        if (optFM==1) ftype="fcr"; else ftype="fmr";
        R_out << "$fsh.mod."<<ftype<<".fully.selected.TCF"<<endl;
        R_out << fTCF_xy(MALE)(styr,endyr-1) << endl;
        R_out << "$fsh.mod."<<ftype<<".fully.selected.SCF"<<endl;
        R_out << fSCF_xy(MALE)(styr,endyr-1) << endl;
        R_out << "$fsh.mod."<<ftype<<".fully.selected.RKF"<<endl;
        R_out << fRKF_xy(MALE)(styr,endyr-1) << endl;
        R_out << "$fsh.mod."<<ftype<<".fully.selected.GTF"<<endl;
        R_out << fGTF_xy(MALE)(styr,endyr-1) <<endl;

        //wts: 20150601: fc's are CAPTURE rates (ONLY calculated if using gmacs calculations)
        if (optFM==1){
            //max rates
            R_out<<"$fsh.fcr.max.TCF.M"<<endl;//new shell, old shell same
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.TCF.NM"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.TCF.OM"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.TCF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcTCFF_yz(i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.SCF.M"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcSCF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.SCF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.RKF.M"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcRKF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.RKF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.GTF.M"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcGTF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.max.GTF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << max(fcGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            //mean rates
            R_out<<"$fsh.fcr.mean.TCF.M"<<endl;//new shell, old shell rates same
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.TCF.NM"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.TCF.OM"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.TCF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcTCFF_yz(i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.SCF.M"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcSCF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.SCF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.RKF.M"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcRKF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.RKF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.GTF.M"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcGTF_xyz(MALE,i)) <<" "; R_out<< endl;
            R_out<<"$fsh.fcr.mean.GTF.F"<<endl;
            for (int i=styr;i<=(endyr-1);i++) R_out << mean(fcGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        }

        //max fishing MORTALITY RATES (changed f... to fm... to clarify: 20150601)
        R_out<<"$fsh.fmr.max.TCF.M"<<endl; //new shell, old shell are identical
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.TCF.NM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.TCF.OM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.TCF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTCFF_yz(i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.SCF.M"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmSCF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.SCF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.RKF.M"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmRKF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.RKF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.GTF.M"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmGTF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.GTF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.TOT.NM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(MALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.TOT.OM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(MALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.TOT.NF"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(FEMALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.max.TOT.OF"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << max(fmTOT_xsyz(FEMALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        //mean fishing MORTALITY RATES (changed f... to fm... to clarify: 20150601)
        R_out<<"$fsh.fmr.mean.TCF.M"<<endl;//new shell, old shell rates same
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.TCF.NM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTCFM_syz(NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.TCF.OM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTCFM_syz(OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.TCF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTCFF_yz(i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.SCF.M"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmSCF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.SCF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmSCF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.RKF.M"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmRKF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.RKF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmRKF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.GTF.M"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmGTF_xyz(MALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.GTF.F"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmGTF_xyz(FEMALE,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.TOT.NM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(MALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.TOT.OM"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(MALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.TOT.NF"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(FEMALE,NEW_SHELL,i)) <<" "; R_out<< endl;
        R_out<<"$fsh.fmr.mean.TOT.OF"<<endl;
        for (int i=styr;i<=(endyr-1);i++) R_out << mean(fmTOT_xsyz(FEMALE,OLD_SHELL,i)) <<" "; R_out<< endl;
        
        R_out <<"$fsh.rmr.max" << endl;
        for (int i=styr;i<=(endyr-1);i++) R_out <<max(fmTCFR_syz(NEW_SHELL,i))<<" "; R_out<<endl; //same as old shell
        R_out <<"$fsh.rmr.mean" << endl;
        for (int i=styr;i<=(endyr-1);i++) R_out <<mean(fmTCFR_syz(NEW_SHELL,i))<<" "; R_out<<endl; //same as old shell
        
        if (optFM==1){
            dvar_vector ratio1(1,nZBs);
            
            //total numbers-at-size captured
            cpN_fyxmsz.initialize();
            for (int yr=styr;yr<=(endyr-1);yr++) {
                //cout<<"yr = "<<yr<<endl;
                //numbers of males captured in TCF
                ratio1 = elem_prod(elem_div(fcTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
                //cout<<"--TCF: males"<<endl;
                //numbers of females captured in TCF
                ratio1 = elem_prod(elem_div(fmTCFF_yz(yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcTCFF_yz(yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcTCFF_yz(yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcTCFF_yz(yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iTCF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
                //cout<<"--TCF: females"<<endl;

                //numbers of males captured in SCF
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
                //cout<<"--SCF: males"<<endl;
                //numbers of females captured in SCF
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcSCF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iSCF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
                //cout<<"--SCF: females"<<endl;

                //numbers of males captured in RKF
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
                //cout<<"--RKF: males"<<endl;
                //numbers of females captured in RKF
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcRKF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iRKF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
                //cout<<"--RKF: females"<<endl;

                //numbers of males captured in GTF
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(MALE,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
                //cout<<"--GTF: males"<<endl;
                //numbers of females captured in GTF
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,NEW_SHELL,yr)),1.0-S_xsyz(FEMALE,NEW_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fcGTF_xyz(FEMALE,yr),fmTOT_xsyz(FEMALE,OLD_SHELL,yr)),1.0-S_xsyz(FEMALE,OLD_SHELL,yr));
                cpN_fyxmsz(iGTF,yr,FEMALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL));
                //cout<<"--GTF: females"<<endl;
            }//yr

            //write out fishery captures for TCF
            //--males
            R_out<<"$fsh.mod.cap.num.TCF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.TCF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.TCF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iTCF,yr,MALE,m,NEW_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.TCF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iTCF,yr,MALE,m,OLD_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            //--females
            R_out<<"$fsh.mod.cap.num.TCF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iTCF,yr,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.TCF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iTCF,yr,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iTCF,yr,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.TCF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iTCF,yr,FEMALE,m,NEW_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.TCF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iTCF,yr,FEMALE,m,OLD_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;

            //write out fishery captures for SCF
            //--males
            R_out<<"$fsh.mod.cap.num.SCF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iSCF,yr,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iSCF,yr,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.SCF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iSCF,yr,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iSCF,yr,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.SCF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iSCF,yr,MALE,m,NEW_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.SCF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iSCF,yr,MALE,m,OLD_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            //--females
            R_out<<"$fsh.mod.cap.num.SCF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iSCF,yr,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.SCF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iSCF,yr,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iSCF,yr,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.SCF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iSCF,yr,FEMALE,m,NEW_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.SCF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iSCF,yr,FEMALE,m,OLD_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            
            //write out fishery captures for RKF
            //--males
            R_out<<"$fsh.mod.cap.num.RKF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iRKF,yr,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iRKF,yr,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.RKF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iRKF,yr,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iRKF,yr,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.RKF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iRKF,yr,MALE,m,NEW_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.RKF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iRKF,yr,MALE,m,OLD_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            //--females
            R_out<<"$fsh.mod.cap.num.RKF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iRKF,yr,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.RKF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iRKF,yr,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iRKF,yr,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.RKF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iRKF,yr,FEMALE,m,NEW_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.RKF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += cpN_fyxmsz(iRKF,yr,FEMALE,m,OLD_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;

            //write out fishery captures for GTF
            //--males
            R_out<<"$fsh.mod.cap.num.GTF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iGTF,yr,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iGTF,yr,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.GTF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iGTF,yr,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iGTF,yr,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.GTF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++)wt += cpN_fyxmsz(iGTF,yr,MALE,m,NEW_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.GTF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++)wt += cpN_fyxmsz(iGTF,yr,MALE,m,OLD_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            //--females
            R_out<<"$fsh.mod.cap.num.GTF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iGTF,yr,FEMALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.num.GTF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iGTF,yr,FEMALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iGTF,yr,FEMALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.GTF.NF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++)wt += cpN_fyxmsz(iGTF,yr,FEMALE,m,NEW_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.cap.bio.GTF.OF"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++)wt += cpN_fyxmsz(iGTF,yr,FEMALE,m,OLD_SHELL)*wt_xmz(FEMALE)(m);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            
        
            //discard rates
            for (int s=NEW_SHELL;s<=OLD_SHELL;s++) {
                for (int yr=styr;yr<=(endyr-1);yr++) fdTCFM_syz(s,yr) = elem_prod((1.0-retFcn_syz(s,yr)),fcTCFM_syz(s,yr));  //discard mortality rate
            }//shell category
            
            //total numbers-at-size discarded
            dsN_fyxmsz.initialize();
            for (int yr=styr;yr<=(endyr-1);yr++) {
                //numbers of males discarded in TCF
                ratio1 = elem_prod(elem_div(fdTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fdTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,IMMATURE,OLD_SHELL));
                ratio1 = elem_prod(elem_div(fdTCFM_syz(NEW_SHELL,yr),fmTOT_xsyz(MALE,NEW_SHELL,yr)),1.0-S_xsyz(MALE,NEW_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,NEW_SHELL));
                ratio1 = elem_prod(elem_div(fdTCFM_syz(OLD_SHELL,yr),fmTOT_xsyz(MALE,OLD_SHELL,yr)),1.0-S_xsyz(MALE,OLD_SHELL,yr));
                dsN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL) = elem_prod(ratio1,modFT_PopNum_yxmsz(yr,MALE,  MATURE,OLD_SHELL));
                dsN_fyxmsz(iTCF,yr,FEMALE) = cpN_fyxmsz(iTCF,yr,FEMALE);
                dsN_fyxmsz(iSCF,yr,  MALE) = cpN_fyxmsz(iSCF,yr,  MALE);
                dsN_fyxmsz(iSCF,yr,FEMALE) = cpN_fyxmsz(iSCF,yr,FEMALE);
                dsN_fyxmsz(iRKF,yr,  MALE) = cpN_fyxmsz(iRKF,yr,  MALE);
                dsN_fyxmsz(iRKF,yr,FEMALE) = cpN_fyxmsz(iRKF,yr,FEMALE);
                dsN_fyxmsz(iGTF,yr,  MALE) = cpN_fyxmsz(iGTF,yr,  MALE);
                dsN_fyxmsz(iGTF,yr,FEMALE) = cpN_fyxmsz(iGTF,yr,FEMALE);
            }
            
            //write out fishery discards for TCF
            //--males
            R_out<<"$fsh.mod.dsc.num.TCF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,NEW_SHELL))+sum(cpN_fyxmsz(iTCF,yr,MALE,MATURE,NEW_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.dsc.num.TCF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) R_out << sum(cpN_fyxmsz(iTCF,yr,MALE,IMMATURE,OLD_SHELL))+sum(cpN_fyxmsz(iTCF,yr,MALE,MATURE,OLD_SHELL)) <<" "; R_out<< endl;
            R_out<<"$fsh.mod.dsc.bio.TCF.NM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += dsN_fyxmsz(iTCF,yr,MALE,m,NEW_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
            R_out<<"$fsh.mod.dsc.bio.TCF.OM"<<endl;
            for (int yr=styr;yr<=(endyr-1);yr++) {
                dvariable wt; wt.initialize();
                for (int m=1;m<=nMSs;m++) wt += dsN_fyxmsz(iTCF,yr,MALE,m,OLD_SHELL)*wt_xmz(MALE,  MATURE);///dot product on z
                R_out << wt <<" "; 
            }  R_out<< endl;
       }//optFM==1
        
        //selectivity and retention curves
        //--surveys
        R_out << "$srv.mod.sel.F"<< endl;
        R_out << selSrv1_xz(FEMALE) << endl; //1974.to.1981
        R_out << selSrv2_xz(FEMALE)<< endl;  //1982.to.1987
        R_out << selSrv3_xz(FEMALE) << endl; //1988.to.endyr
        R_out << "$srv.mod.sel.M"<< endl;
        R_out << selSrv1_xz(MALE) << endl;
        R_out << selSrv2_xz(MALE)<< endl;
        R_out << selSrv3_xz(MALE) << endl;
        //--fisheries
        REP2R2(fsh.mod.sel.TCF.M,selTCFM_syz(NEW_SHELL));
        REP2R2(fsh.mod.selr.TCF.M,selTCFR_syz(NEW_SHELL));
        REP2R2(fsh.mod.ret.TCF.M,retFcn_syz(NEW_SHELL));
        REP2R2(fsh.mod.sel.TCF.F,selTCFF_z);
        R_out << "$fsh.mod.sel.SCF.F"<< endl;
        R_out <<selSCF_cxz(1,FEMALE)<<endl;
        R_out <<selSCF_cxz(2,FEMALE)<<endl;
        R_out <<selSCF_cxz(3,FEMALE)<<endl;  
        R_out << "$fsh.mod.sel.SCF.M"<< endl;
        R_out <<selSCF_cxz(1,MALE)<<endl;
        R_out <<selSCF_cxz(2,MALE)<<endl;
        R_out <<selSCF_cxz(3,MALE)<<endl;
        R_out << "$fsh.mod.sel.RKF.F"<< endl;
        R_out <<selRKF_cxz(1,FEMALE)<<endl;
        R_out <<selRKF_cxz(2,FEMALE)<<endl;
        R_out <<selRKF_cxz(3,FEMALE)<<endl;  
        R_out << "$fsh.mod.sel.RKF.M"<< endl;
        R_out <<selRKF_cxz(1,MALE)<<endl;
        R_out <<selRKF_cxz(2,MALE)<<endl;
        R_out <<selRKF_cxz(3,MALE)<<endl;  
        R_out << "$fsh.mod.sel.GTF.F"<< endl;
        R_out <<selGTF_cxz(1,FEMALE)<<endl;
        R_out <<selGTF_cxz(2,FEMALE)<<endl;
        R_out <<selGTF_cxz(3,FEMALE)<<endl;  
        R_out << "$fsh.mod.sel.GTF.M"<< endl;
        R_out <<selGTF_cxz(1,MALE)<<endl;
        R_out <<selGTF_cxz(2,MALE)<<endl;
        R_out <<selGTF_cxz(3,MALE)<<endl;
        
        //input & effective sample sizes
        REP2R2(srv.inpSS,ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,FEMALE));//same for all m,s,x for survey
//        REP2R2(srv.inpSS.IM,ssObsZCsSrv_msxn(IMMATURE,NEW_SHELL,  MALE));
//        REP2R2(srv.inpSS.MF,ssObsZCsSrv_msxn(  MATURE,NEW_SHELL,FEMALE));
//        REP2R2(srv.inpSS.MM,ssObsZCsSrv_msxn(  MATURE,NEW_SHELL,  MALE));
        REP2R2(srv.effSS.McI,effnSrv_y);
        
        REP2R2(fsh.inpSS.ret.TCF,ssRetZCsTCF_sn(NEW_SHELL));
        REP2R2(fsh.inpSS.tot.TCF.M,ssTotZCsTCFM_sn(NEW_SHELL));  
        REP2R2(fsh.inpSS.tot.TCF.F,ssZCsTCFF_n);  
        REP2R2(fsh.inpSS.tot.SCF.M,ssZCsSCFM_sn(NEW_SHELL));  
        REP2R2(fsh.inpSS.tot.SCF.F,ssZCsSCFF_n);  
        REP2R2(fsh.inpSS.tot.RKF.M,ssZCsRKFM_sn(NEW_SHELL));  
        REP2R2(fsh.inpSS.tot.RKF.F,ssZCsRKFF_n);
        REP2R2(fsh.inpSS.tot.GTF,ssObsZCsGTF_n);  
        
        REP2R2(fsh.effSS.McI.ret.TCF,effnTCF_ret_y);
        REP2R2(fsh.effSS.McI.tot.TCF.F,effnTCF_tot_xy(FEMALE));
        REP2R2(fsh.effSS.McI.tot.TCF.M,effnTCF_tot_xy(  MALE));
        REP2R2(fsh.effSS.McI.tot.SCF.F,effnSCF_tot_xy(FEMALE));
        REP2R2(fsh.effSS.McI.tot.SCF.M,effnSCF_tot_xy(  MALE));
        REP2R2(fsh.effSS.McI.tot.RKF.F,effnRKF_tot_xy(FEMALE));
        REP2R2(fsh.effSS.McI.tot.RKF.M,effnRKF_tot_xy(  MALE));
        REP2R2(fsh.effSS.McI.tot.GTF,effnGTF_tot_y);
        
        //Z-scores
        //--surveys
        dmatrix zsc(1,nSXs,styr,endyr);
        zsc.initialize();
        for (int x=1;x<=nSXs;x++){
            for (int n=1;n<=nObsSrvBio;n++){
                int y = yrsObsSrvBio_n(n);
                zsc(x,y) = value(zsSrvMatBio_xn(x,n));
            }
        }
        REP2R2(srv.bio.zscr.F,zsc(FEMALE));
        REP2R2(srv.bio.zscr.M,zsc(  MALE));
        
        //--fisheries
        //cout<<"TCF"<<endl;
        REP2R2(fsh.ret.zscr.TCF,zsRetMortBio_TCFR_y);
        zsc.allocate(1,nSXs,styr,endyr-1);
        zsc.initialize();
        for (int n=1;n<=nObsDscTCF;n++){
            int y = yrsObsDscTCF_n(n);
            //cout<<n<<tb<<y<<endl;
            if (optTCFMfit==0){
                zsc(MALE,y) = value(zsTotMortBio_TCFM_n(n));
            } else {
                zsc(MALE,y) = value(zsDscMortBio_TCFM_n(n));
            }
            zsc(FEMALE,y) = value(zsDscMortBio_TCFF_n(n));
        }
        REP2R2(fsh.bio.zscr.TCF.F,zsc(FEMALE));
        REP2R2(fsh.bio.zscr.TCF.M,zsc(  MALE));
        //cout<<"SCF"<<endl;
        zsc.initialize();
        for (int n=1;n<=nObsDscSCF;n++){
            int y = yrsObsDscSCF(n);
            //cout<<n<<tb<<y<<endl;
            for (int x=1;x<=nSXs;x++) zsc(x,y) = value(zsDscMortBio_SCF_xn(x,n));
        }
        REP2R2(fsh.bio.zscr.SCF.F,zsc(FEMALE));
        REP2R2(fsh.bio.zscr.SCF.M,zsc(  MALE));
        //cout<<"RKF"<<endl;
        zsc.initialize();
        for (int n=1;n<=nObsDscRKF;n++){
            int y = yrsObsDscRKF(n);
            //cout<<n<<tb<<y<<endl;
            for (int x=1;x<=nSXs;x++) zsc(x,y) = value(zsDscMortBio_RKF_xn(x,n));
        }
        REP2R2(fsh.bio.zscr.RKF.F,zsc(FEMALE));
        REP2R2(fsh.bio.zscr.RKF.M,zsc(  MALE));
        //cout<<"GTF"<<endl;
        zsc.initialize();
        for (int n=1;n<=nObsDscGTF;n++){
            int y = yrsObsDscGTF(n);
            //cout<<n<<tb<<y<<endl;
            zsc(1,y) = value(zsDscMortBio_GTF_n(n));
        }
        REP2R2(fsh.bio.zscr.GTF,zsc(1));

        //write out fishery numbers/biomass-at-size
        if (optFM==1){
            myWriteN_fyxmsz(R_out,iTCF,"TCF","cap",cpN_fyxmsz);
            myWriteN_fyxmsz(R_out,iSCF,"SCF","cap",cpN_fyxmsz);
            myWriteN_fyxmsz(R_out,iRKF,"RKF","cap",cpN_fyxmsz);
            myWriteN_fyxmsz(R_out,iGTF,"GTF","cap",cpN_fyxmsz);
            myWriteN_fyxmsz(R_out,iTCF,"TCF","dsc",dsN_fyxmsz);
            myWriteN_fyxmsz(R_out,iSCF,"SCF","dsc",dsN_fyxmsz);
            myWriteN_fyxmsz(R_out,iRKF,"RKF","dsc",dsN_fyxmsz);
            myWriteN_fyxmsz(R_out,iGTF,"GTF","dsc",dsN_fyxmsz);
        }
        dvar6_array rmN_fyxmsz(1,nFsh,styr,endyr-1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
        rmN_fyxmsz.initialize();
        for (int y=styr;y<=(endyr-1);y++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++) rmN_fyxmsz(iTCF,y,MALE,m,s) = rmN_ymsz(y,m,s);
            }
        }
        myWriteN_fyxmsz(R_out,iTCF,"TCF","rm",rmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iSCF,"SCF","rm",rmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iRKF,"RKF","rm",rmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iGTF,"GTF","rm",rmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iTCF,"TCF","dm",dmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iSCF,"SCF","dm",dmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iRKF,"RKF","dm",dmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iGTF,"GTF","dm",dmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iTCF,"TCF","tm",tmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iSCF,"SCF","tm",tmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iRKF,"RKF","tm",tmN_fyxmsz);
        myWriteN_fyxmsz(R_out,iGTF,"GTF","tm",tmN_fyxmsz);
//    cout<<"done writeToR_OLD"<<endl;

// ==========================================================================
FUNCTION void myWriteN_yxmsz(ostream& os, dvar5_array& a_yxmsz)
        //write out fishery array
        //population abundance (millions)
        os<<"$pop.mod.NatZ.INM"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << a_yxmsz(yr,MALE,IMMATURE,NEW_SHELL)<<endl;
        os<<"$pop.mod.NatZ.MNM"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << a_yxmsz(yr,MALE,  MATURE,NEW_SHELL)<<endl;
        os<<"$pop.mod.NatZ.MOM"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << a_yxmsz(yr,MALE,  MATURE,OLD_SHELL)<<endl;
        os<<"$pop.mod.NatZ.INF"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << a_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL)<<endl;
        os<<"$pop.mod.NatZ.MNF"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << a_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL)<<endl;
        os<<"$pop.mod.NatZ.MOF"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << a_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL)<<endl;
        //population biomass (1000's t)
        os<<"$pop.mod.BatZ.INM"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << elem_prod(a_yxmsz(yr,MALE,IMMATURE,NEW_SHELL),wt_xmz(MALE,IMMATURE))<<endl;
        os<<"$pop.mod.BatZ.MNM"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << elem_prod(a_yxmsz(yr,MALE,  MATURE,NEW_SHELL),wt_xmz(MALE,  MATURE))<<endl;
        os<<"$pop.mod.BatZ.MOM"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << elem_prod(a_yxmsz(yr,MALE,  MATURE,OLD_SHELL),wt_xmz(MALE,  MATURE))<<endl;
        os<<"$pop.mod.BatZ.INF"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << elem_prod(a_yxmsz(yr,FEMALE,IMMATURE,NEW_SHELL),wt_xmz(FEMALE,IMMATURE))<<endl;
        os<<"$pop.mod.BatZ.MNF"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << elem_prod(a_yxmsz(yr,FEMALE,  MATURE,NEW_SHELL),wt_xmz(FEMALE,  MATURE))<<endl;
        os<<"$pop.mod.BatZ.MOF"<<endl;
        for (int yr=styr;yr<=endyr;yr++) os << elem_prod(a_yxmsz(yr,FEMALE,  MATURE,OLD_SHELL),wt_xmz(FEMALE,  MATURE))<<endl;
        
// ==========================================================================
FUNCTION void myWriteN_fyxmsz(ostream& os, int iF, adstring aF, adstring type, dvar6_array& a_fyxmsz)
        //write out fishery array
        //capture abundance (millions)
        os<<"$fsh.mod."<<type<<".NatZ."<<aF<<".INM"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << a_fyxmsz(iF,yr,MALE,IMMATURE,NEW_SHELL)<<endl;
        os<<"$fsh.mod."<<type<<".NatZ."<<aF<<".MNM"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << a_fyxmsz(iF,yr,MALE,  MATURE,NEW_SHELL)<<endl;
        os<<"$fsh.mod."<<type<<".NatZ."<<aF<<".MOM"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << a_fyxmsz(iF,yr,MALE,  MATURE,OLD_SHELL)<<endl;
        os<<"$fsh.mod."<<type<<".NatZ."<<aF<<".INF"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << a_fyxmsz(iF,yr,FEMALE,IMMATURE,NEW_SHELL)<<endl;
        os<<"$fsh.mod."<<type<<".NatZ."<<aF<<".MNF"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << a_fyxmsz(iF,yr,FEMALE,  MATURE,NEW_SHELL)<<endl;
        os<<"$fsh.mod."<<type<<".NatZ."<<aF<<".MOF"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << a_fyxmsz(iF,yr,FEMALE,  MATURE,OLD_SHELL)<<endl;
        //capture BatZmass (1000's t)
        os<<"$fsh.mod."<<type<<".BatZ."<<aF<<".INM"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << elem_prod(a_fyxmsz(iF,yr,MALE,IMMATURE,NEW_SHELL),wt_xmz(MALE,IMMATURE))<<endl;
        os<<"$fsh.mod."<<type<<".BatZ."<<aF<<".MNM"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << elem_prod(a_fyxmsz(iF,yr,MALE,  MATURE,NEW_SHELL),wt_xmz(MALE,  MATURE))<<endl;
        os<<"$fsh.mod."<<type<<".BatZ."<<aF<<".MOM"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << elem_prod(a_fyxmsz(iF,yr,MALE,  MATURE,OLD_SHELL),wt_xmz(MALE,  MATURE))<<endl;
        os<<"$fsh.mod."<<type<<".BatZ."<<aF<<".INF"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << elem_prod(a_fyxmsz(iF,yr,FEMALE,IMMATURE,NEW_SHELL),wt_xmz(FEMALE,IMMATURE))<<endl;
        os<<"$fsh.mod."<<type<<".BatZ."<<aF<<".MNF"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << elem_prod(a_fyxmsz(iF,yr,FEMALE,  MATURE,NEW_SHELL),wt_xmz(FEMALE,  MATURE))<<endl;
        os<<"$fsh.mod."<<type<<".BatZ."<<aF<<".MOF"<<endl;
        for (int yr=styr;yr<=(endyr-1);yr++) os << elem_prod(a_fyxmsz(iF,yr,FEMALE,  MATURE,OLD_SHELL),wt_xmz(FEMALE,  MATURE))<<endl;
        
// ==========================================================================
FUNCTION void myWriteParamsToR(ostream& os)
    cout<<"starting myWriteParamsToR"<<endl;
    adstring strp;
    os<<"params=list(";
        os<<"growth=list(";
            os<<"pGrAF1="<<pGrAF1<<cc;
            os<<"pGrBF1="<<pGrBF1<<cc;
            os<<"pGrAM1="<<pGrAM1<<cc;
            os<<"pGrBM1="<<pGrBM1;
        os<<")"<<cc;
        os<<"natMort=list(";
            os<<"Mmult.imm="<<pMfac_Imm<<cc;
            os<<"Mmult.m="<<pMfac_MatM<<cc;
            os<<"Mmult.f="<<pMfac_MatF<<cc;
            os<<"big.mort="; wts::writeToR(os,value(pMfac_Big),dmX); 
        os<<")"<<cc;
        os<<"recruitment=list(";
            strp = "y="+str(mnYrRecCurr)+":"+str(endyr);
            os<<"pMnLnRec="<<pMnLnRec<<cc;
            os<<"pRecDevs="; wts::writeToR(os,value(pRecDevs),strp);  os<<cc;
            strp = "y="+str(mnYrRecDevsHist)+":"+str(mnYrRecCurr-1);
            os<<"pMnLnRecInit="<<pMnLnRecInit<<cc;
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
                sel(1) = pRetTCFM_z50A1; sel(2) = pRetTCFM_z50A2;
                slp(1) = pRetTCFM_slpA1; slp(2) = pRetTCFM_slpA2;
                os<<"retention=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                os<<"male=list(";
                        sel = mfexp(pSelTCFM_mnLnZ50A2);
                        slp(1) = pSelTCFM_slpA1; slp(2) = pSelTCFM_slpA2;
                        strp = qt+str(1)+":"+str(nSelTCFM_devsZ50)+qt;
                        os<<"z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<cc<<"devs.lnSel50="; wts::writeToR(os,value(pSelTCFM_devsZ50),strp);
                os<<"),";
                os<<"female=list(z50="<<pSelTCFF_z50<<cc<<"slope="<<pSelTCFF_slp<<")";
            os<<"),";
            }
            {
            os<<"scf=list(";
                dvar_vector sel(1,3); dvar_vector slp(1,3);
                sel(1) = pSelSCFF_z50A1; sel(2) = pSelSCFF_z50A2; sel(3) = pSelSCFF_z50A3;
                slp(1) = pSelSCFF_slpA1; slp(2) = pSelSCFF_slpA2; slp(3) = pSelSCFF_slpA3;
                os<<"female=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                os<<"male=list(";
                    sel(1) = pSelSCFM_z50A1; sel(2) = pSelSCFM_z50A2; sel(3) = pSelSCFM_z50A3;
                    slp(1) = pSelSCFM_slpA1; slp(2) = pSelSCFM_slpA2; slp(3) = pSelSCFM_slpA3;
                    os<<"ascending.limb=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                    sel(1) = pSelSCFM_lnZ50D1; sel(2) = pSelSCFM_lnZ50D2; sel(3) = pSelSCFM_lnZ50D3;
                    slp(1) = pSelSCFM_slpD1; slp(2) = pSelSCFM_slpD2; slp(3) = pSelSCFM_slpD3;
                    os<<"descending.limb=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<")";
                os<<")";
            os<<"),"<<endl;
            }
            {
            os<<"rkf=list(";
                dvar_vector sel(1,3); dvar_vector slp(1,3);
                sel(1) = pSelRKFF_z50A1; sel(2) = pSelRKFF_z50A2; sel(3) = pSelRKFF_z50A3;
                slp(1) = pSelRKFF_slpA1; slp(2) = pSelRKFF_slpA2; slp(3) = pSelRKFF_slpA3;
                os<<"female=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                sel(1) = pSelRKFM_z50A1; sel(2) = pSelRKFM_z50A2; sel(3) = pSelRKFM_z50A3;
                slp(1) = pSelRKFM_slpA1; slp(2) = pSelRKFM_slpA2; slp(3) = pSelRKFM_slpA3;
                os<<"male=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<")";
            os<<"),";
            }
            {
            os<<"gtf=list(";
                dvar_vector sel(1,3); dvar_vector slp(1,3);
                sel(1) = pSelGTFF_z50A1; sel(2) = pSelGTFF_z50A2; sel(3) = pSelGTFF_z50A3;
                slp(1) = pSelGTFF_slpA1; slp(2) = pSelGTFF_slpA2; slp(3) = pSelGTFF_slpA3;
                os<<"female=list(z50="; wts::writeToR(os,value(sel)); os<<", slope="; wts::writeToR(os,value(slp)); os<<"),";
                sel(1) = pSelGTFM_z50A1; sel(2) = pSelGTFM_z50A2; sel(3) = pSelGTFM_z50A3;
                slp(1) = pSelGTFM_slpA1; slp(2) = pSelGTFM_slpA2; slp(3) = pSelGTFM_slpA3;
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
            nAtZ(x,NEW_SHELL,IMMATURE,y) = value(modNum_yxmsz(y,x,IMMATURE,NEW_SHELL));
            nAtZ(x,NEW_SHELL,  MATURE,y) = value(modNum_yxmsz(y,x,  MATURE,NEW_SHELL));
            nAtZ(x,OLD_SHELL,IMMATURE,y) = value(modNum_yxmsz(y,x,IMMATURE,OLD_SHELL));
            nAtZ(x,OLD_SHELL,  MATURE,y) = value(modNum_yxmsz(y,x,  MATURE,OLD_SHELL));
        }//yr
    }//x
    
    os<<"mod.pop=list("<<endl;
        os<<"mnYr="<<mnYr<<cc<<"mxYr="<<mxYr<<cc<<"recLag="<<recLag<<cc;
        os<<"mnYrRecDevsHist="<<mnYrRecDevsHist<<cc<<"mnYrRecCurr="<<mnYrRecCurr;   os<<cc<<endl;
        os<<"rec="; wts::writeToR(os,value(rec_y),dmY);                             os<<cc<<endl;
        os<<"MMB="; wts::writeToR(os,value(modSpBioMateTime_xy(   MALE)),dmYm1);    os<<cc<<endl;
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
    //cout<<"TCF"<<endl;
    for (int n=1;n<=nObsDscTCF;n++){
        int y = yrsObsDscTCF_n(n);
        //cout<<n<<tb<<y<<endl;
        if (optTCFMfit==0){
            zsTMB_TCFM_y(y) = value(zsTotMortBio_TCFM_n(n));
        } else {
            zsTMB_TCFM_y(y) = value(zsDscMortBio_TCFM_n(n));
        }
        zsDMB_TCFF_y(y) = value(zsDscMortBio_TCFF_n(n));
    }
    //cout<<"SCF"<<endl;
    for (int n=1;n<=nObsDscSCF;n++){
        int y = yrsObsDscSCF(n);
        //cout<<n<<tb<<y<<endl;
        for (int x=1;x<=nSXs;x++) zsDMB_SCF_xy(x,y) = value(zsDscMortBio_SCF_xn(x,n));
    }
    //cout<<"RKF"<<endl;
    for (int n=1;n<=nObsDscRKF;n++){
        int y = yrsObsDscRKF(n);
        //cout<<n<<tb<<y<<endl;
        for (int x=1;x<=nSXs;x++) zsDMB_RKF_xy(x,y) = value(zsDscMortBio_RKF_xn(x,n));
    }
    //cout<<"GTF"<<endl;
    for (int n=1;n<=nObsDscGTF;n++){
        int y = yrsObsDscGTF(n);
        //cout<<n<<tb<<y<<endl;
        zsDMB_GTF_y(y) = value(zsDscMortBio_GTF_n(n));
    }
    
    os<<"fsh=list("<<endl;
        os<<"legalSize="<<zLegal<<cc<<endl;
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
            os<<"eff=list("<<endl;
                os<<"q="<<qSCF<<cc<<endl;
                os<<"F="; wts::writeToR(os,value(mfexp(pAvgLnF_SCF+pF_DevsSCF))); os<<cc<<endl;
                os<<"E="; wts::writeToR(os,effSCF_y(yrsObsDscSCF)); os<<endl;
            os<<")"; os<<cc<<endl;
            os<<"fits=list("<<endl;
                os<<"zscr="; wts::writeToR(os,zsDMB_SCF_xy,dmX,dmYm1);        os<<cc<<endl;
                os<<"like="; wts::writeToR(os,value(lkDscMortBio_SCF_x),dmX); os<<cc<<"type='bycatch'"<<endl;
            os<<")"; os<<endl;
        os<<")"; os<<cc<<endl;
        os<<"RKF=list("<<endl;
            os<<"eff=list("<<endl;
                os<<"q="<<qRKF<<cc<<endl;
                os<<"F="; wts::writeToR(os,value(mfexp(pAvgLnF_RKF+pF_DevsRKF))); os<<cc<<endl;
                os<<"E="; wts::writeToR(os,effRKF_y(yrsObsDscRKF)); os<<endl;
            os<<")"; os<<cc<<endl;
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
FUNCTION void writeToR_NEW(ostream& os)
    cout<<"starting writeToR_NEW"<<endl;
    os<<"res<-list("<<endl;
        ptrMDS->writeToR(os,adstring("model.data")); os<<cc<<endl;
        myWriteParamsToR(os);                        os<<cc<<endl;
        myWriteModPopInfoToR(os);                    os<<cc<<endl;
        myWriteModFshInfoToR(os);                    os<<cc<<endl;
        myWriteModSrvInfoToR(os);                    os<<endl;
    os<<")"<<endl;
    
    cout<<"finished writeToR_NEW"<<endl;
    
// ==========================================================================
FUNCTION void myWriteMCMCToR(ostream& os)
    cout<<"starting myWriteMCMCToR"<<endl;
    
    d4_array nAtZlast(1,nSXs,1,nSCs,1,nMSs,1,nZBs);
    for (int x=1;x<=nSXs;x++){
        nAtZlast(x,NEW_SHELL,IMMATURE) = value(modNum_yxmsz(endyr,x,IMMATURE,NEW_SHELL));
        nAtZlast(x,NEW_SHELL,  MATURE) = value(modNum_yxmsz(endyr,x,  MATURE,NEW_SHELL));
        nAtZlast(x,OLD_SHELL,IMMATURE) = value(modNum_yxmsz(endyr,x,IMMATURE,OLD_SHELL));
        nAtZlast(x,OLD_SHELL,  MATURE) = value(modNum_yxmsz(endyr,x,  MATURE,OLD_SHELL));
    }
    
    os<<"list("<<endl;
        myWriteParamsToR(os); os<<cc<<endl;
        os<<"rec="; wts::writeToR(os,value(rec_y),dmY); os<<cc<<endl;
        os<<"MMB="; wts::writeToR(os,value(modSpBioMateTime_xy(   MALE)),dmYm1); os<<cc<<endl;
        os<<"nAtZ.last="; wts::writeToR(os,nAtZlast,dmX,dmS,dmM,dmZ); os<<cc<<endl;
        
    os<<"dummy=0),"<<endl;
    
    cout<<"finished myWriteMCMCToR"<<endl;
    
// ==========================================================================
FUNCTION void writeLikelihoodComponents(ostream& os, int toR)
    os<<"----------------------------------------------------------------------------"<<endl;
    os<<"Likelihood components"<<endl;
    os<<"----------------------------------------------------------------------------"<<endl;
    os<<f<<cc<<tb<<"objective function value"<<endl;
    os<<"idx,   weight,      likelihood,      objFun,    category, description"<<endl;
    for (int i=1;i<=NUM_FOUT;i++){
        os<<i<<cc<<wgtsOut(i)<<cc<<likeOut(i)<<cc<<objfOut(i)<<cc<<tb<<strFOUT(i)<<endl;
    }
  
// ==========================================================================
// ==========================================================================
REPORT_SECTION
    cout<<"starting REPORT_SECTION for phase "<<current_phase()<<endl;
    
    cout<<"qSCF = "<<qSCF<<tb<<"ln(qSCF) = "<<log(qSCF)<<endl;
    cout<<"qRKF = "<<qRKF<<tb<<"ln(qSCF) = "<<log(qRKF)<<endl;
    
    if (active(pSelTCFM_devsZ50)) { 
        double llw = 0.0;
        double red = 0.01;
        int phs = pSelTCFM_devsZ50.get_phase_start();
        llw = 1.0; llw = pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))*llw;
        cout<<"llw for pSelTCFM_devsZ50:"<<endl;
        cout<<current_phase()<<tb<<phs<<tb<<maxPhase<<endl;
        cout<<(current_phase()-phs)<<tb<<max(1.0,1.0*(maxPhase-phs))<<endl;
        cout<<pow(red,(current_phase()-phs)/max(1.0,1.0*(maxPhase-phs)))<<endl;
     }
    if (last_phase()) {
        ofstream os("TCSAM2013.NEWSTYLE.final.R", ios::trunc);    
        writeToR_NEW(os);
        os.close();

        os.open("TCSAM2013.OLDSTYLE.final.R");
        writeToR_OLD(os);
        os.close();
        
        os.open("TCSAM2013ProjMod.dat");
        writeMyProjectionFile(os);
        os.close();
    
        os.open("TCSAM2013.params.all.final.csv");
        writeParameters(os,0,0);
        os.close();
        os.open("TCSAM2013.params.active.final.csv");
        writeParameters(os,0,1);
        os.close();
        
        os.open("TCSAM2013.final_likelihood_components.csv");
        writeLikelihoodComponents(os,0);
        os.close();
        
        if (option_match(ad_comm::argc,ad_comm::argv,"-jitter")>-1) {
            ofstream fs("jitterInfo.csv");
            fs<<"seed"<<cc<<"objfun"<<endl;
            fs<<iSeed<<cc<<f<<endl;
        }
        
        //do OFL calculations
        if (doOFL){
            cout<<"ReportToR: starting OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
            calcOFL(asmtYr,debugOFL,echoOFL);//updates oflResults
            oflResults.writeCSVHeader(echoOFL); echoOFL<<endl;
            oflResults.writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL.close();
//            os<<","<<endl;
//            oflResults.writeToR(os,"oflResults",0);
            cout<<"ReportToR: finished OFL calculations"<<endl;
        }
    }
//  report << offset << endl;
    cout<<"finished REPORT_SECTION"<<endl;

// ===============================================================================
// ===============================================================================
BETWEEN_PHASES_SECTION
// ===============================================================================
// ===============================================================================
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 500,1000,3000,3000,5000,5000,10000
  convergence_criteria 1,1,.01,.001,.001,.001,1e-3,1e-4

// ===============================================================================
// ===============================================================================
TOP_OF_MAIN_SECTION
  arrmblsize = 10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(5000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(150000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(7000);
  time(&start);
  CheckFile.open("CheckFile.dat");
  

// ===============================================================================
// ===============================================================================
FINAL_SECTION
  
    if (sd_phase()>0){
        cout<<"nvarcalc                  = "<<initial_params::nvarcalc()<<endl;
        cout<<"nvar_calc_all             = "<<initial_params::nvarcalc_all()<<endl;
        cout<<"num_initial_params        = "<<initial_params::num_initial_params<<endl;
        cout<<"num_active_calc           = "<<initial_params::num_active_calc()<<endl;
        adlist_ptr varsptr = initial_params::varsptr;
        int ctr = 0;
        for (int i=0;i<initial_params::num_initial_params;i++){
            pinitial_params p = varsptr[i];
            int cnt = p->size_count();
            int phs = p->get_phase_start();
            cout<<i<<tb<<p->label()<<tb<<"count = "<<cnt<<tb<<"phase = "<<phs<<endl;
            dvector x(1,cnt); int ii=1;
            p->copy_value_to_vector(x,ii);
            ctr += ii;
            cout<<"value = "<<x<<endl;
        }
        cout<<"ctr = "<<ctr<<endl;
        cout<<"nvarcalc                  = "<<initial_params::nvarcalc()<<endl;
        cout<<"nvar_calc_all             = "<<initial_params::nvarcalc_all()<<endl;
        cout<<"num_initial_params        = "<<initial_params::num_initial_params<<endl;
        cout<<"num_active_calc           = "<<initial_params::num_active_calc()<<endl;
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

//-------------OFL Calculations--------------
FUNCTION void calcOFL(int yr, int debug, ostream& cout)
    if (debug) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcOFL(yr,debug,cout)"<<endl;
        cout<<"year for projection = "<<yr<<endl;
    }

    //1. get initial population for "upcoming" year, yr
    d4_array n_xmsz = wts::value(modNum_yxmsz(yr));
    if (debug) {cout<<"  males_msz:"<<endl; wts::print(n_xmsz(  MALE),cout,1);}
    if (debug) {cout<<"females_msz:"<<endl; wts::print(n_xmsz(FEMALE),cout,1);}
    
    //2. determine average recruitment
    dvector avgRec_x(1,nSXs);
    avgRec_x = 0.5*value(mean(rec_y(1982,yr)));
    if (debug) {
        cout<<"R_y( 1982:"<<yr<<")      = "<<rec_y(1982,yr)<<endl;
        cout<<"R_yx(1982:"<<yr<<",MALE) = "<<0.5<<endl;
        cout<<"Average recruitment = "<<avgRec_x<<endl;
    }
    
    //3. set yr back one year to get population rates, etc.
    yr = yr-1;//don't have pop rates, etc. for endyr
    if (debug) cout<<"year for pop rates = "<<yr<<endl;
    
    //4. Determine population rates for next year, using yr
    double dtF = mdptFshs_y(yr);//time at which fisheries occur
    double dtM = dtF+spmo;//time at which mating occurs (spmo is offset from dtF)
    
    PopDyInfo* pPIM = new PopDyInfo(nZBs,dtM);//  males info
    pPIM->R_z   = value(prRec_z);
    pPIM->w_mz  = ptrMDS->pBio->wAtZ_xmz(MALE);
    for (int m=1;m<=nMSs;m++) {
        for (int s=1;s<=nSCs;s++) pPIM->M_msz(m,s) = value(M_msx(m,s,MALE));
    }
    for (int s=1;s<=nSCs;s++) pPIM->T_szz(s) = trans(value(prGr_xzz(MALE)));//need to transpose growth matrix
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = value(modPrM2M(MALE));
    
    PopDyInfo* pPIF = new PopDyInfo(nZBs,dtM);//females info
    pPIF->R_z   = value(prRec_z);
    pPIF->w_mz  = ptrMDS->pBio->wAtZ_xmz(FEMALE);
    for (int m=1;m<=nMSs;m++) {
        for (int s=1;s<=nSCs;s++) pPIF->M_msz(m,s) = value(M_msx(m,s,FEMALE));
    }
    for (int s=1;s<=nSCs;s++) pPIF->T_szz(s) = trans(value(prGr_xzz(FEMALE)));//need to transpose growth matrix
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = value(modPrM2M(FEMALE));
    if (debug) cout<<"calculated pPIM, pPIF."<<endl;
    
    //5. Determine fishery conditions for next year based on averages for recent years
        int oflAvgPeriodYrs = 1;
        //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
        //a. Calculate average handling mortality, retention curves and capture rates
        int ny;   //number of years fishery is active
        dvector avgHM_f(1,nFsh);
        avgHM_f.initialize();
        avgHM_f(iTCF) = hm_pot;
        avgHM_f(iSCF) = hm_pot;
        avgHM_f(iRKF) = hm_pot;
        avgHM_f(iGTF) = hm_trawl;
        if (debug) cout<<"avgHm_f = "<<avgHM_f<<endl;
    if (debug) cout<<"5.1"<<endl;

        d5_array avgCapF_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged capture mortality
        d5_array avgRFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged retention function
        //no need to average retention function
        avgRFcn_xfmsz.initialize();
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++) avgRFcn_xfmsz(MALE,iTCF,m,s) = value(retFcn_syz(s,yr));
        }
        avgCapF_xfmsz.initialize();
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++) {
                for (int z=1;z<=nZBs;z++){
                    ny = 0;
                    for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                        int nchk = 0;
                        for (int ychk=yrsObsDscTCF_n.indexmin();ychk<=yrsObsDscTCF_n.indexmax();ychk++) {
                            if (y==yrsObsDscTCF_n(ychk)){
                                ny++;
                                avgCapF_xfmsz(MALE,iTCF,m,s,z) += value(fcTCFM_syz(s,y,z));
                            }
                        }
                    }//y
                    avgCapF_xfmsz(MALE,iTCF,m,s,z) /= 1.0*ny;
                }//z
            }//s
        }//m
    if (debug) cout<<"5.2"<<endl;
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++) {
                for (int z=1;z<=nZBs;z++){
                    ny = 0;
                    for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                        int nchk = 0;
                        for (int ychk=yrsObsDscTCF_n.indexmin();ychk<=yrsObsDscTCF_n.indexmax();ychk++) {
                            if (y==yrsObsDscTCF_n(ychk)){
                                ny++;
                                avgCapF_xfmsz(FEMALE,iTCF,m,s,z) += value(fcTCFF_yz(y,z));
                            }
                        }
                    }//y
                    avgCapF_xfmsz(FEMALE,iTCF,m,s,z) /= 1.0*ny;
                }//z
            }//s
        }//m
    if (debug) cout<<"5.3"<<endl;
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++) {
                    for (int z=1;z<=nZBs;z++){
                        ny = 0;
                        for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                            int nchk = 0;
                            for (int ychk=yrsObsDscSCF.indexmin();ychk<=yrsObsDscSCF.indexmax();ychk++) {
                                if (y==yrsObsDscSCF(ychk)){
                                    ny++;
                                    avgCapF_xfmsz(x,iSCF,m,s,z) += value(fcSCF_xyz(x,y,z));
                                }
                            }
                        }//y
                        avgCapF_xfmsz(x,iSCF,m,s,z) /= 1.0*ny;
                        ny = 0;
                        for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                            int nchk = 0;
                            for (int ychk=yrsObsDscRKF.indexmin();ychk<=yrsObsDscRKF.indexmax();ychk++) {
                                if (y==yrsObsDscRKF(ychk)){
                                    ny++;
                                    avgCapF_xfmsz(x,iRKF,m,s,z) += value(fcRKF_xyz(x,y,z));
                                }
                            }
                        }//y
                        avgCapF_xfmsz(x,iRKF,m,s,z) /= 1.0*ny;
                        ny = 0;
                        for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                            int nchk = 0;
                            for (int ychk=yrsObsDscGTF.indexmin();ychk<=yrsObsDscGTF.indexmax();ychk++) {
                                if (y==yrsObsDscGTF(ychk)){
                                    ny++;
                                    avgCapF_xfmsz(x,iGTF,m,s,z) += value(fcGTF_xyz(x,y,z));
                                }
                            }
                        }//y
                        avgCapF_xfmsz(x,iGTF,m,s,z) /= 1.0*ny;
                    }//z
                }//s
            }//m
        }//x
        if (debug){
            for (int f=1;f<=nFsh;f++){
                cout<<"avgCapF_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgCapF_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
            }
        }
        
        CatchInfo* pCIM = new CatchInfo(nZBs,nFsh,dtF);//male catch info
        pCIM->setCaptureRates(avgCapF_xfmsz(MALE));
        pCIM->setRetentionFcns(avgRFcn_xfmsz(MALE));
        pCIM->setHandlingMortality(avgHM_f);
        double maxCapF = pCIM->findMaxTargetCaptureRate(cout);
        if (debug) cout<<"maxCapF = "<<maxCapF<<endl;
        
        CatchInfo* pCIF = new CatchInfo(nZBs,nFsh,dtF);//female catch info
        pCIF->setCaptureRates(avgCapF_xfmsz(FEMALE));
        pCIF->setRetentionFcns(avgRFcn_xfmsz(FEMALE));
        pCIF->setHandlingMortality(avgHM_f);
        pCIF->maxF = maxCapF;//need to set this for females
    if (debug) cout<<"5.4"<<endl;
        
    //6. Create PopProjectors
        PopProjector* pPPM = new PopProjector(pPIM,pCIM);
        PopProjector* pPPF = new PopProjector(pPIF,pCIF);
        if (debug) cout<<"created pPPs."<<endl;
        
    //7. Create Equilibrium_Calculators
        Equilibrium_Calculator* pECM = new Equilibrium_Calculator(pPPM);
        Equilibrium_Calculator* pECF = new Equilibrium_Calculator(pPPF);
        if (debug) cout<<"created pECs."<<endl;
        
    //8. define OFL_Calculator
        OFL_Calculator*  pOC;
        if (debug) cout<<"declared pOC."<<endl;
        
    //9. Determine TIER LEVEL, create Tier_Calculators, calculate OFL
        int tier = 3;
        if (tier==3){
            Tier3_Calculator* pT3CM = new Tier3_Calculator(0.35,pECM);
            Tier3_Calculator* pT3CF = new Tier3_Calculator(0.35,pECF);
            if (debug) cout<<"created pT3Cs."<<endl;
            pOC = new OFL_Calculator(pT3CM,pT3CF);
            if (debug) {
                cout<<"created pOC."<<endl;
//                OFL_Calculator::debug=1;
//                Tier3_Calculator::debug=1;
//                Equilibrium_Calculator::debug=1;
            }
            oflResults = pOC->calcOFLResults(avgRec_x,n_xmsz,cout);
            if (debug) {
                cout<<"calculated oflResults."<<endl;
                oflResults.writeCSVHeader(cout); cout<<endl;
                oflResults.writeToCSV(cout); cout<<endl;
//                OFL_Calculator::debug=0;
//                Tier3_Calculator::debug=0;
//                Equilibrium_Calculator::debug=0;
            }
        }//Tier 3 calculation
    
    if (debug) {
        int n = 100;
        MultiYearPopProjector* pMYPPM = new MultiYearPopProjector(pPPM);
        MultiYearPopProjector* pMYPPF = new MultiYearPopProjector(pPPF);
        pMYPPM->projectUnFished(n,avgRec_x(  MALE),n_xmsz(  MALE),cout);
        double myPP_B0 = pMYPPM->matBio_y(n);
        pMYPPM->project(n,avgRec_x(  MALE),oflResults.Fmsy,n_xmsz(  MALE),cout);
        pMYPPF->project(n,avgRec_x(FEMALE),oflResults.Fmsy,n_xmsz(FEMALE),cout);
        double myPP_Bmsy = pMYPPM->matBio_y(n);
        double myPP_prjB = pMYPPM->matBio_y(1);
        double myPP_MSY  = pMYPPM->totCM_y(n)+pMYPPF->totCM_y(n);
        double myPP_OFL  = pMYPPM->totCM_y(1)+pMYPPF->totCM_y(1);
        cout<<"#------------------------"<<endl;
        cout<<"MYPP B0     = "<<myPP_B0  <<". B0     = "<<oflResults.B0  <<endl;
        cout<<"MYPP Bmsy   = "<<myPP_Bmsy<<". Bmsy   = "<<oflResults.Bmsy<<endl;
        cout<<"MYPP prjMMB = "<<myPP_prjB<<". prjMMB = "<<oflResults.prjB<<endl;
        cout<<"MYPP MSY    = "<<myPP_MSY <<". MSY    = "<<oflResults.MSY <<endl;
        cout<<"MYPP OFL    = "<<myPP_OFL <<". OFL    = "<<oflResults.OFL <<endl;
        cout<<"#------------------------"<<endl;
        cout<<"finished calcOFL(yr,debug,cout)"<<endl<<endl<<endl;
    }
        

