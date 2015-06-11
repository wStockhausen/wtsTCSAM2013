#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"

//**********************************************************************
//  Includes
//      ModelConfiguration
//**********************************************************************
//--------------------------------------------------------------------------------
//          ModelConfiguration
//--------------------------------------------------------------------------------
int ModelConfiguration::debug=0;
/***************************************************************
*   creation                                                   *
***************************************************************/
ModelConfiguration::ModelConfiguration(){
    asmtYr=mnYr=mxYr=0;
    inclTSD=inclGTF=inclSCF=inclRKF=inclTCF=1;//default to TRUE
    runOpMod=fitToPriors=1;//default to TRUE
    jitFrac=0.0;
}
/***************************************************************
*   destruction                                                *
***************************************************************/
ModelConfiguration::~ModelConfiguration(){
    if (debug) std::cout<<"destroying ModelConfiguration "<<this<<std::endl;
    if (debug) std::cout<<"destruction complete "<<std::endl;
}

/***************************************************************
*   function to read from file in ADMB format                  *
***************************************************************/
void ModelConfiguration::read(const adstring & fn) {
    if (debug) std::cout<<"ModelConfiguration::read(fn). Reading from '"<<fn<<"'"<<std::endl;
    cifstream strm(fn);
    read(strm);
    if (debug) std::cout<<"end ModelConfiguration::read(fn). Read from '"<<fn<<"'"<<std::endl;
}

/***************************************************************
*   function to write to file in ADMB format                   *
***************************************************************/
void ModelConfiguration::write(const adstring & fn) {
    if (debug) std::cout<<"#start ModelConfiguration::write(fn). Writing to '"<<fn<<"'"<<std::endl;
    std::ofstream strm(fn,std::ofstream::out|std::ofstream::trunc);
    write(strm); //write to file
    strm.close();
    if (debug) std::cout<<"#end ModelConfiguration::write(fn). Wrote to '"<<fn<<"'"<<std::endl;
}

/***************************************************************
*   function to read from file in ADMB format                  *
***************************************************************/
void ModelConfiguration::read(cifstream & is) {
    if (debug) std::cout<<"ModelConfiguration::read(cifstream & is)"<<std::endl;
    adstring parent = wts::getParentFolder(is);
    cout<<"Model configuration file is '"<<is.get_file_name()<<"'"<<endl;
    cout<<"Parent folder is '"<<parent<<endl;
    
    is>>cfgName;
    if (debug) std::cout<<cfgName<<std::endl;
    is>>asmtYr;//assessment year (can be < maxYr)
    is>>mnYr; //min model year
    is>>mxYr; //max model year
    is>>jitFrac; //jitter fraction
    is>>nZBs;
    if (debug){
        std::cout<<asmtYr<<tb<<"#assessment year"<<std::endl;
        std::cout<<mnYr <<tb<<"#model min year"<<std::endl;
        std::cout<<mxYr <<tb<<"#model max year"<<std::endl;
        std::cout<<jitFrac<<tb<<"#jitter fraction"<<std::endl;
        std::cout<<nZBs<<tb<<"#number of size bins"<<std::endl;
    }
    zBins.allocate(1,nZBs); 
    zBinCutPts.allocate(1,nZBs+1); 
    onesZBins.allocate(1,nZBs); onesZBins = 1.0;
    is>>zBinCutPts;
    for (int z=1;z<=nZBs;z++) zBins(z) = 0.5*(zBinCutPts(z)+zBinCutPts(z+1));
    if (debug){
        std::cout<<"#size bins (mm CW)"<<std::endl;
        std::cout<<zBins<<std::endl;
        std::cout<<"#size bin cut points (mm CW)"<<std::endl;
        std::cout<<zBinCutPts <<std::endl;
        std::cout<<"enter 1 to continue : ";
        std::cin>>debug;
        if (debug<0) exit(1);
    }
        
    adstring str;
    is>>str; inclTSD    = ModelConsts::getBooleanType(str);//include trawl survey data?
    is>>str; inclTCF    = ModelConsts::getBooleanType(str);//include directed Tanner fishery?
    is>>str; inclSCF    = ModelConsts::getBooleanType(str);//include opi fishery bycatch?
    is>>str; inclRKF    = ModelConsts::getBooleanType(str);//include red king crab fishery bycatch?
    is>>str; inclGTF    = ModelConsts::getBooleanType(str);//include groundfish trawl fishery bycatch?
    is>>str; runOpMod   = ModelConsts::getBooleanType(str);//run population model?
    is>>str; fitToPriors = ModelConsts::getBooleanType(str);//fit priors?
    
    is>>fnCtl;
    is>>fnMDS;
    fnCtl = wts::concatenateFilePaths(parent,fnCtl);
    fnMDS = wts::concatenateFilePaths(parent,fnMDS);
    
    if (debug){
        std::cout<<ModelConsts::getBooleanType(inclTSD)    <<"   #include trawl survey data?"<<std::endl;
        std::cout<<ModelConsts::getBooleanType(inclTCF)    <<"   #include directed tanner crab fishery?"<<std::endl;
        std::cout<<ModelConsts::getBooleanType(inclSCF)    <<"   #include snow crab fishery?"<<std::endl;
        std::cout<<ModelConsts::getBooleanType(inclRKF)    <<"   #include BBred king crab fishery?"<<std::endl;
        std::cout<<ModelConsts::getBooleanType(inclGTF)    <<"   #include groundfish trawl fishery?"<<std::endl;
        std::cout<<ModelConsts::getBooleanType(runOpMod)   <<"   #run operating model?"<<std::endl;
        std::cout<<ModelConsts::getBooleanType(fitToPriors)<<"   #fit to priors?"<<std::endl;
    }
        std::cout<<fnCtl<<"   #model control file"<<std::endl;
        std::cout<<fnMDS<<"   #model datasets file"<<std::endl;
    if (debug){
        std::cout<<"enter 1 to continue : ";
        std::cin>>debug;
        if (debug<0) exit(1);
    }
    
    if (debug) std::cout<<"end ModelConfiguration::read(cifstream & is)"<<std::endl;
}

/***************************************************************
*   function to write to file in ADMB format                   *
***************************************************************/
void ModelConfiguration::write(std::ostream & os) {
    if (debug) std::cout<<"#start ModelConfiguration::write(ostream)"<<std::endl;
    os<<"#######################################################"<<std::endl;
    os<<"#TCSAM2013 Model Configuration File                   #"<<std::endl;
    os<<"#######################################################"<<std::endl;
    os<<cfgName<<tb<<"#Model configuration name"<<std::endl;
    os<<asmtYr<<tb<<"#Assessment year"<<std::endl;
    os<<mnYr<<tb<<mxYr<<tb<<"#Min, max model years"<<std::endl;
    os<<jitFrac<<tb<<"#jitter fraction"<<std::endl;
    os<<nZBs<<tb<<"#Number of model size classes"<<std::endl;
    os<<"#size bin cut points"<<std::endl;
    os<<zBinCutPts <<std::endl;
    
    os<<ModelConsts::getBooleanType(inclTSD)<<tb<<"#include trawl survey data?"<<std::endl;
    os<<ModelConsts::getBooleanType(inclTCF)<<tb<<"#include directed Tanner fishery?"<<std::endl;
    os<<ModelConsts::getBooleanType(inclSCF)<<tb<<"#include opilio fishery bycatch?"<<std::endl;
    os<<ModelConsts::getBooleanType(inclRKF)<<tb<<"#include red king crab fishery bycatch?"<<std::endl;
    os<<ModelConsts::getBooleanType(inclGTF)<<tb<<"#include groudnfish trawl fishery bycatch?"<<std::endl;
    
    os<<ModelConsts::getBooleanType(runOpMod)<<tb<<"#run operating model?"<<std::endl;
    os<<ModelConsts::getBooleanType(fitToPriors)<<tb<<"#fit priors?"<<std::endl;
    
    os<<fnCtl<<tb<<"#Model parameters configuration file"<<std::endl;
    os<<fnMDS<<tb<<"#Model datasets file"<<std::endl;

    if (debug) std::cout<<"#end ModelConfiguration::write(ostream)"<<std::endl;
}

/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void ModelConfiguration::writeToR(std::ostream& os, char* nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list("<<std::endl;
    indent++;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"configName='"<<cfgName<<"', mnYr="<<mnYr<<", mxYr="<<mxYr<<", assYr="<<asmtYr<<cc<<"jitFrac="<<jitFrac<<cc<<std::endl;
        os<<"nZBins="<<nZBs<<cc;
        os<<"zBins=";      wts::writeToR(os,zBins);     os<<cc<<std::endl;
        os<<"zBinCutPts="; wts::writeToR(os,zBinCutPts);os<<cc<<std::endl;
//        writeVector(os,zBinCutPts,"zBinCutpts"); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"flags=list(";
        os<<"incFshTnr="<<inclTCF<<cc;
        os<<"incFshOpi="<<inclSCF<<cc;
        os<<"incFshRKC="<<inclRKF<<cc;
        os<<"incFshGrf="<<inclGTF<<cc;
        os<<"runOpMod="<<runOpMod<<cc;
        os<<"fitToPriors="<<fitToPriors<<"),";
        os<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"fnMPC='"<<fnCtl<<"', fnMDS='"<<fnMDS<<"')"<<std::endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
        os<<")";
}
/////////////////////////////////end ModelConfiguration/////////////////////////

