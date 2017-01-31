#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "FisheryData.hpp"
#include "ModelData.hpp"
//**********************************************************************
//  Includes
//      BioData
//      TrawlSurveyData
//      ModelDatasets
//**********************************************************************
int BioData::debug         = 0;
int TrawlSurveyData::debug = 0;
int ModelDatasets::debug   = 0;
//----------------------------------------------------------------------
//          BioData
//----------------------------------------------------------------------
/***************************************************************
*   read.                                                      *
***************************************************************/
void BioData::read(cifstream & is){
    if (debug) {
        cout<<"start BioData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading Bio Data."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    //SIZE BINS
    is>>nZBins;
    if (debug) cout<<nZBins<<tb<<"#nZBins"<<endl;
    zBins.allocate(1,nZBins);
    is>>zBins;
    if (debug) cout<<zBins<<tb<<"#zBins (mm CW)"<<endl;
    
    //WEIGHT-AT-SIZE
    {is>>unitsWatZ;
    if (debug) cout<<unitsWatZ<<tb<<"#unitsWatZ"<<endl;
    wAtZ_xmz.allocate(FEMALE,MALE,IMMATURE,MATURE,1,nZBins);
    wAtZ_xmz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    if (debug) cout<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = ModelConsts::getSexType(factors(1));
        int mat   = ModelConsts::getMaturityType(factors(2));
        if (debug) cout<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<endl;
        if (sex||mat){
            is>>wAtZ_xmz(sex,mat);
            if (debug) {
                cout<<"#wAtZ_xmz("<<factors(1)<<cc<<factors(2)<<")="<<endl;
                cout<<"#"<<zBins<<endl;
                cout<<wAtZ_xmz(sex,mat)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for wAtZ not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //PROBABILITY OF MATURING AT SIZE
    {prMature_xz.allocate(FEMALE,MALE,1,nZBins);
    prMature_xz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    if (debug) cout<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,1);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = ModelConsts::getSexType(factors(1));
        if (debug) cout<<factors(1)<<tb<<"#factors"<<endl;
        if (sex){
            is>>prMature_xz(sex);
            if (debug) {
                cout<<"#prMature_xz("<<factors(1)<<")="<<endl;
                cout<<"#"<<zBins<<endl;
                cout<<prMature_xz(sex)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for prMature_xz not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //FRACTION MATURE BY SEX, SHELL CONDITION
    {frMature_xsz.allocate(FEMALE,MALE,IMMATURE,MATURE,1,nZBins);
    frMature_xsz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    if (debug) cout<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = ModelConsts::getSexType(factors(1));
        int shl   = ModelConsts::getShellType(factors(2));
        if (debug) cout<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<endl;
        if (sex||shl){
            is>>frMature_xsz(sex,shl);
            if (debug) {
                cout<<"#frMature_xsz("<<factors(1)<<cc<<factors(2)<<")="<<endl;
                cout<<"#"<<zBins<<endl;
                cout<<frMature_xsz(sex,shl)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for frMature_xsz not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //CV for mean size at min, max size bins
    {if (debug) cout<<"#CV for min, max sizes"<<endl;
    cvMnMxZ_xc.allocate(FEMALE,MALE,1,2);
    cvMnMxZ_xc.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    if (debug) cout<<nc<<tb<<"#number of factor combinations"<<endl;
    adstring_array factors(1,1);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = ModelConsts::getSexType(factors(1));
        if (debug) cout<<factors(1)<<"#factors"<<endl;
        if (sex){
            is>>cvMnMxZ_xc(sex);
            if (debug) {
                cout<<"#cvMnMxZ_xc("<<factors(1)<<")="<<endl;
                cout<<cvMnMxZ_xc(sex)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for cvMnMxZ_xc not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    
    //MIDPOINTS OF FISHERY SEASONS BY YEAR
    if (debug) cout<<"#Fishing season midpoints"<<endl;
    is>>nyFshSeasons;
    if (debug) cout<<nyFshSeasons<<tb<<"#number of years"<<endl;
    mdptFshSeasons_yc.allocate(1,nyFshSeasons,1,2);
    is>>mdptFshSeasons_yc;
    if (debug) cout<<"#mdptFshSeasons"<<endl<<mdptFshSeasons_yc<<endl;
    
    if (debug) cout<<"end BioData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void BioData::write(ostream & os){
    if (debug) cout<<"start BioData::write(...) "<<this<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<"#Tanner crab bio data"<<endl;
    
    os<<"#-----------SIZE BINS---------------------------------#"<<endl;
    os<<nZBins<<tb<<"#number of size bins"<<endl;
    os<<"#size bins (mm CW)"<<endl<<zBins<<endl;
    
    {os<<"#-----------WEIGHT-AT-SIZE----------------------------#"<<endl;
    os<<unitsWatZ<<tb<<"#units for weight-at-size"<<endl;
    os<<MALE*MATURE<<tb<<"#number of factor combinations (sex x maturity state)"<<endl;
    adstring_array factors(1,2);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = ModelConsts::getSexType(sex);
        for (int mat=IMMATURE;mat<=MATURE;mat++){
            factors(2) = ModelConsts::getMaturityType(mat);
            os<<"#-------"<<factors(1)<<cc<<factors(2)<<endl;
            os<<factors(1)<<tb<<factors(2)<<endl;
            os<<wAtZ_xmz(sex,mat)<<endl;
        }
    }}
    
    {os<<"#-----------PROBABILITY OF MATURING-------------------#"<<endl;
    os<<MALE<<tb<<"#number of factor combinations (sex)"<<endl;
    adstring_array factors(1,1);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = ModelConsts::getSexType(sex);
        os<<"#-------"<<factors(1)<<endl;
        os<<factors(1)<<endl;
        os<<prMature_xz(sex)<<endl;
    }}
    
    {os<<"#-----------FRACTION MATURE BY SEX, SHELL CONDITION---#"<<endl;
    os<<MALE*OLD_SHELL<<tb<<"#number of factor combinations (sex x shell condition)"<<endl;
    adstring_array factors(1,2);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = ModelConsts::getSexType(sex);
        for (int shl=NEW_SHELL;shl<=OLD_SHELL;shl++){
            factors(2) = ModelConsts::getShellType(shl);
            os<<"#-------"<<factors(1)<<cc<<factors(2)<<endl;
            os<<factors(1)<<tb<<factors(2)<<endl;
            os<<frMature_xsz(sex,shl)<<endl;
        }
    }}
    
    {os<<"#-----------CV FOR MIN, MAX SIZES---------------------#"<<endl;
    os<<MALE<<tb<<"#number of factor combinations (sex)"<<endl;
    adstring_array factors(1,1);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = ModelConsts::getSexType(sex);
        os<<"#-------"<<factors(1)<<endl;
        os<<factors(1)<<endl;
        os<<cvMnMxZ_xc(sex)<<endl;
    }}
    
    {os<<"#-----------FISHERY SEASON MIDPOINTS------------------#"<<endl;
    os<<nyFshSeasons<<tb<<"#number of years"<<endl;
    os<<"#year  midpt"<<endl;
    os<<mdptFshSeasons_yc<<endl;}
    
    if (debug) cout<<"end BioData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void BioData::writeToR(ostream& os, char* nm, int indent) {
    adstring x = "x=c("+qt+STR_FEMALE+qt+cc+qt+STR_MALE+qt+")";
    adstring s = "s=c("+qt+STR_IMMATURE+qt+cc+qt+STR_MATURE+qt+")";
    adstring z = "z=c("+wts::to_qcsv(zBins)+")";
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list("<<endl;
    indent++;
        //size bins
        for (int n=0;n<indent;n++) os<<tb;
        os<<"zBins=";wts::writeToR(os,zBins); os<<","<<endl;
        
        //weight-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"wAtZ=list(units="<<qt<<unitsWatZ<<qt<<cc<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<endl;
            wts::writeToR(os,wAtZ_xmz,x,s,z); os<<endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        
        //probability of maturing-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"prMat="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            wts::writeToR(os,prMature_xz,x,z); os<<","<<endl;
        indent--;}
        
        //fraction mature-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"frMat="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            wts::writeToR(os,frMature_xsz,x,s,z); os<<","<<endl;
        indent--;}
        
        //cv for min, max sizes
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"cvZs="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring c = "cols=c('minZ','maxZ')";
            wts::writeToR(os,cvMnMxZ_xc,x,c); os<<","<<endl;
        indent--;}
        
        //fishery season midpoints
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"mdptFshSeasons="<<endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            adstring cols = "cols=c('year','midpt')";
            adstring yrs  = "y=c("+wts::to_qcsv(wts::to_ivector(column(mdptFshSeasons_yc,1)))+")";
            wts::writeToR(os,mdptFshSeasons_yc,yrs,cols); os<<endl;
        indent--;}
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
}
/////////////////////////////////end BioData/////////////////////////
//----------------------------------------------------------------------
//          TrawlSurveyData
//----------------------------------------------------------------------
/***************************************************************
*   read.                                                      *
***************************************************************/
void TrawlSurveyData::read(cifstream & is){
    if (debug) {
        cout<<"start TrawlSurveyData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading Trawl Survey Data."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    //ABUNDANCE
    is>>nyAbund;//number of years of abundance data
    if (debug) cout<<nyAbund<<tb<<"#nyAbund"<<endl;
    is>>unitsAbund;
    if (debug) cout<<unitsAbund<<tb<<"#unitsAbund"<<endl;
    inpAbund_yc.allocate(1,nyAbund,1,5);
    is>>inpAbund_yc;
    if (debug) cout<<"#abund = "<<endl<<"#year abundance cv_females cv_males"<<endl<<inpAbund_yc<<endl;
    
    yrsAbund.allocate(1,nyAbund);
    abund_y.allocate(1,nyAbund);
    abund_xy.allocate(FEMALE,MALE,1,nyAbund);
    cvsAbund_xy.allocate(FEMALE,MALE,1,nyAbund);
    
    yrsAbund = (ivector) column(inpAbund_yc,1);
    abund_y.initialize();
    for (int x=FEMALE;x<=MALE;x++) {
        abund_xy(x) = column(inpAbund_yc,1+x);
        cvsAbund_xy(x) = column(inpAbund_yc,3+x);
        
        abund_y += abund_xy(x);
    }
    
    //Mature BIOMASS
    is>>nyMatBio;//number of years of biomass data
    if (debug) cout<<nyMatBio<<tb<<"#nyMatBio"<<endl;
    is>>unitsMatBio;
    if (debug) cout<<unitsMatBio<<tb<<"#unitsMatBio"<<endl;
    inpMatBio_yc.allocate(1,nyMatBio,1,5);
    is>>inpMatBio_yc;
    if (debug) cout<<"#mature biomass = "<<endl<<"#year females males cv_females cv_males"<<endl<<inpMatBio_yc<<endl;
    
    yrsMatBio.allocate(1,nyMatBio);
    matBio_xy.allocate(FEMALE,MALE,1,nyMatBio);
    cvsMatBio_xy.allocate(FEMALE,MALE,1,nyMatBio);
    
    yrsMatBio = (ivector) column(inpMatBio_yc,1);
    for (int x=FEMALE;x<=MALE;x++) {
        matBio_xy(x)    = column(inpMatBio_yc,1+x);
        cvsMatBio_xy(x) = column(inpMatBio_yc,3+x);
    }
    
    //NUMBERS-AT-SIZE 
    is>>nyNatZ;//number of years of numbers-at-size data
    if (debug) cout<<nyNatZ<<tb<<"#nyNatZ"<<endl;
    is>>unitsNatZ;
    if (debug) cout<<unitsNatZ<<tb<<"#unitsNatZ"<<endl;
    is>>nZCutPts;
    if (debug) cout<<nZCutPts<<tb<<"#nZCutPts"<<endl;
    zCutPts.allocate(1,nZCutPts);
    is>>zCutPts;
    if (debug) cout<<zCutPts<<tb<<"#zCutPts (mm CW)"<<endl;
    zBins.allocate(1,nZCutPts-1);
    zBins = 0.5*(zCutPts(1,nZCutPts-1)+(--zCutPts(2,nZCutPts)));
    if (debug) cout<<zBins<<tb<<"#zBins"<<endl;
    
    yrsNatZ.allocate(1,nyNatZ);
    ssNatZ_xsmy.allocate(FEMALE,MALE,NEW_SHELL,OLD_SHELL,IMMATURE,MATURE,1,nyNatZ);
    nAtZ_xsmyz.allocate(FEMALE,MALE,NEW_SHELL,OLD_SHELL,IMMATURE,MATURE,1,nyNatZ,1,nZCutPts-1);
    inpNatZ_xsmyc.allocate(FEMALE,MALE,NEW_SHELL,OLD_SHELL,IMMATURE,MATURE,1,nyNatZ,1,nZCutPts-1+2);
    
    yrsNatZ.initialize();
    ssNatZ_xsmy.initialize();
    nAtZ_xsmyz.initialize();
    inpNatZ_xsmyc.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    adstring_array factors(1,3);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = ModelConsts::getSexType(factors(1));
        int shell = ModelConsts::getShellType(factors(2));
        int mat   = ModelConsts::getMaturityType(factors(3));
        if (debug) cout<<factors(1)<<tb<<factors(2)<<tb<<factors(3)<<tb<<"#factors"<<endl;
        if (sex||shell||mat){
            is>>inpNatZ_xsmyc(sex,shell,mat);
            if (debug) {
                cout<<"#inpNatZ_xsmyc("<<factors(1)<<cc<<factors(2)<<cc<<factors(3)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                cout<<inpNatZ_xsmyc(sex,shell,mat)<<endl;
            }
            yrsNatZ = (ivector) column(inpNatZ_xsmyc(sex,shell,mat),1);
            ssNatZ_xsmy(sex,shell,mat) = column(inpNatZ_xsmyc(sex,shell,mat),2);
            for (int y=1;y<=nyNatZ;y++){
                nAtZ_xsmyz(sex,shell,mat,y)  = (inpNatZ_xsmyc(sex,shell,mat,y)(3,nZCutPts-1+2)).shift(1);
            }
            if (debug) {
                cout<<"#year, sample_size, nAtZ_xsmyz("<<factors(1)<<cc<<factors(2)<<cc<<factors(3)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                for (int y=1;y<=nyNatZ;y++) cout<<yrsNatZ(y)<<tb<<ssNatZ_xsmy(sex,shell,mat,y)<<tb<<nAtZ_xsmyz(sex,shell,mat,y)<<endl;
                cout<<"Enter '1' to continue>> ";
                int resp;
                cin>>resp;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for nAtZ not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<tb<<factors(3)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }
    if (debug) cout<<"end TrawlSurveyData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void TrawlSurveyData::write(ostream & os){
    if (debug) cout<<"start TrawlSurveyData::write(...) "<<this<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<"#Tanner crab trawl survey data"<<endl;
    os<<"#-----------ABUNDANCE---------------#"<<endl;
    os<<nyAbund<<tb<<"#number of years of survey abundance data"<<endl;
    os<<unitsAbund<<tb<<"#units for crab numbers"<<endl;
    os<<"#year females males cv_females cv_males"<<endl<<inpAbund_yc<<endl;
    os<<"#-----------MATURE BIOMASS---------------#"<<endl;
    os<<nyMatBio<<tb<<"#number of years of survey mature biomass data"<<endl;
    os<<unitsMatBio<<tb<<"#units for mature biomass"<<endl;
    os<<"#year females males cv_females cv_males"<<endl<<inpMatBio_yc<<endl;
    os<<"#-----------NUMBERS-AT-SIZE DATA----------------------#"<<endl;
    os<<nyNatZ<<tb<<"#number of years of size data"<<endl;
    os<<unitsNatZ<<tb<<"#units for numbers at size of crab"<<endl;
    os<<nZCutPts<<tb<<"#number of size bin cutpoints"<<endl;
    os<<"#size bin cutpoints (mm CW)"<<endl<<zCutPts<<endl;
    os<<MALE*OLD_SHELL*MATURE<<tb<<"#number of sex x shell x maturity factor combinations to read in"<<endl;
    adstring_array factors(1,3);
    for (int sex=FEMALE;sex<=MALE;sex++){
        factors(1) = ModelConsts::getSexType(sex);
        for (int shell=NEW_SHELL;shell<=OLD_SHELL;shell++){
            factors(2) = ModelConsts::getShellType(shell);
            for (int mat=IMMATURE;mat<=MATURE;mat++){
                factors(3) = ModelConsts::getMaturityType(mat);
                os<<"#-------"<<factors(1)<<cc<<factors(2)<<cc<<factors(3)<<endl;
                os<<factors(1)<<tb<<factors(2)<<tb<<factors(3)<<endl;
                os<<"#year sample_size "<<zBins<<endl;
                os<<inpNatZ_xsmyc(sex,shell,mat)<<endl;
            }
        }
    }
    if (debug) cout<<"end TrawlSurveyData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void TrawlSurveyData::writeToR(ostream& os, char* nm, int indent) {
    if (debug) cout<<"TrawlSurveyData::writing to R"<<endl;
    adstring x = "x=c("+qt+STR_FEMALE+qt+cc+qt+STR_MALE+qt+")";
    adstring m = "m=c("+qt+STR_IMMATURE+qt +cc+ qt+STR_MATURE+qt+")";
    adstring s = "s=c("+qt+STR_IMMATURE+qt+cc+qt+STR_MATURE+qt+")";
    adstring z = "z=c("+wts::to_qcsv(zBins)+")";
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list(name="<<qt<<"trawl survey"<<qt<<cc<<endl;
        //abundance
        indent++;
        {   for (int n=0;n<indent;n++) os<<tb;
            adstring y = "yrs=c("+wts::to_qcsv(yrsAbund)+")";
            os<<"abundance=list(units="<<qt<<unitsAbund<<qt<<cc<<endl;
            indent++; 
                for (int n=0;n<indent;n++) os<<tb;
                os<<"years="; wts::writeToR(os,yrsAbund); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"cvs="; wts::writeToR(os,cvsAbund_xy,x,y); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"data="; wts::writeToR(os,abund_xy,x,y); os<<endl;
            indent--;
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        indent--;}
        
        //mature biomass
        indent++;
        {   for (int n=0;n<indent;n++) os<<tb;
            adstring y = "yrs=c("+wts::to_qcsv(yrsMatBio)+")";
            os<<"biomass=list(units="<<qt<<unitsMatBio<<qt<<cc<<endl;
            indent++; 
                for (int n=0;n<indent;n++) os<<tb;
                os<<"years="; wts::writeToR(os,yrsMatBio); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"cvs="; wts::writeToR(os,cvsMatBio_xy,x,y); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"data="; wts::writeToR(os,matBio_xy,x,y); os<<endl;
            indent--;
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        indent--;}
        
        //catch nAtZ
        {indent++;   
            for (int n=0;n<indent;n++) os<<tb;
            adstring y = "y=c("+wts::to_qcsv(yrsNatZ)+")";
            os<<"nAtZ=list(units="<<qt<<unitsNatZ<<qt<<cc<<endl; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"years="; wts::writeToR(os,yrsNatZ); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"cutpts="; wts::writeToR(os,zCutPts); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"sample.sizes="; wts::writeToR(os,ssNatZ_xsmy,x,s,m,y); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<endl;
            wts::writeToR(os,nAtZ_xsmyz,x,s,m,y,z); os<<endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<")"<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
    if (debug) cout<<"TrawlSurveyData::done writing to R"<<endl;
}
/////////////////////////////////end TrawlSurveyData/////////////////////////
/***************************************************************
*   Instantiation                                              *
***************************************************************/
ModelDatasets::ModelDatasets(ModelConfiguration* ptrMC){
    pMC=ptrMC;
    pBio=0;pTSD=0;
    pTCFR=0;pTCFD=0;pSCF=0;pRKF=0;pGTF=0;
}
/***************************************************************
*   Destruction                                                *
***************************************************************/
ModelDatasets::~ModelDatasets(){
    pMC=0;
    delete pBio;  pBio=0;
    delete pTSD;  pTSD=0;
    delete pTCFR; pTCFR=0;
    delete pTCFD; pTCFD=0;
    delete pSCF; pSCF=0;
    delete pRKF; pRKF=0;
    delete pGTF; pGTF=0;
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void ModelDatasets::read(cifstream & is){
    if (debug) {
        cout<<"start ModelDatasets::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#datasets file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    adstring parent = wts::getParentFolder(is);
    cout<<"Model datasets file is '"<<is.get_file_name()<<"'"<<endl;
    cout<<"Parent folder is '"<<parent<<endl;
    
    is>>fnBioData;
    is>>fnSurveyData_Trawl;
    is>>fnFisheryData_TCFR;
    is>>fnFisheryData_TCFD;
    is>>fnFisheryData_SnCF;
    is>>fnFisheryData_RKCF;
    is>>fnFisheryData_GrTF;

    fnBioData = wts::concatenateFilePaths(parent,fnBioData);
    fnSurveyData_Trawl = wts::concatenateFilePaths(parent,fnSurveyData_Trawl);
    fnFisheryData_TCFR = wts::concatenateFilePaths(parent,fnFisheryData_TCFR);
    fnFisheryData_TCFD = wts::concatenateFilePaths(parent,fnFisheryData_TCFD);
    fnFisheryData_SnCF = wts::concatenateFilePaths(parent,fnFisheryData_SnCF);
    fnFisheryData_RKCF = wts::concatenateFilePaths(parent,fnFisheryData_RKCF);
    fnFisheryData_GrTF = wts::concatenateFilePaths(parent,fnFisheryData_GrTF);
    
 //   if (debug){
        cout<<fnBioData         <<tb<<"#tanner crab biological data file"<<endl;
        cout<<fnSurveyData_Trawl<<tb<<"#trawl survey data file"<<endl;
        cout<<fnFisheryData_TCFR<<tb<<"#directed tanner crab fishery data (retained) file"<<endl;
        cout<<fnFisheryData_TCFD<<tb<<"#directed tanner crab fishery data (discarded) file"<<endl;
        cout<<fnFisheryData_SnCF<<tb<<"#snow crab fishery data file"<<endl;
        cout<<fnFisheryData_RKCF<<tb<<"#BB red king crab fishery data file"<<endl;
        cout<<fnFisheryData_GrTF<<tb<<"#groundfish trawl fishery data file"<<endl;
//    }
    
    //          Bio data
    {pBio = new BioData();
    cifstream strm(fnBioData,ios::in);
    strm>>(*pBio);
    if (debug){
        cout<<endl<<"#Bio Data"<<endl;
        cout<<(*pBio)<<endl;
    }}  
    
    //          Survey data
    //Groundfish trawl survey data
    {pTSD = new TrawlSurveyData();
    cifstream strm(fnSurveyData_Trawl,ios::in);
    strm>>(*pTSD);
    if (debug){
        cout<<endl<<"#Trawl Survey Data"<<endl;
        cout<<(*pTSD)<<endl;
    }}
    
    //          Fishery data
    //Directed fishery (retained catch)
    {pTCFR = new RetainedFisheryData();
    cifstream strm(fnFisheryData_TCFR,ios::in);
    strm>>(*pTCFR);
    if (debug){
        cout<<endl<<"#Directed Fishery Data (retained)"<<endl;
        cout<<(*pTCFR)<<endl;
    }}
    
    //Directed fishery (discarded catch)
    {pTCFD = new DiscardFisheryData();
    cifstream strm(fnFisheryData_TCFD,ios::in);
    strm>>(*pTCFD);
    if (debug){
        cout<<endl<<"#Directed Fishery Data (discarded)"<<endl;
        cout<<(*pTCFD)<<endl;
    }}
    
    //Snow crab fishery (discarded catch)
    {pSCF = new DiscardFisheryData();
    cifstream strm(fnFisheryData_SnCF,ios::in);
    strm>>(*pSCF);
    if (debug){
        cout<<endl<<"#Snow Crab Fishery Data"<<endl;
        cout<<(*pSCF)<<endl;
    }}
    
    //BB Red King Crab (discarded catch)
    {pRKF = new DiscardFisheryData();
    cifstream strm(fnFisheryData_RKCF,ios::in);
    strm>>(*pRKF);
    if (debug){
        cout<<endl<<"#BB Red King Crab Fishery Data (discarded)"<<endl;
        cout<<(*pRKF)<<endl;
    }}
    
    
    //Groundfish Trawl Fishery (discarded catch)
    {pGTF = new GroundfishTrawlFisheryData();
    cifstream strm(fnFisheryData_GrTF,ios::in);
    strm>>(*pGTF);
    if (debug){
        cout<<endl<<"#Groundfish Trawl Fishery (discarded)"<<endl;
        cout<<(*pGTF)<<endl;
    }}
    if (debug) cout<<"end ModelDatasets::read(...) "<<this<<endl;
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void ModelDatasets::write(ostream & os){
    if (debug) cout<<"start ModelDatasets::write(...) "<<this<<endl;
    os<<fnBioData<<tb<<"#tanner crab biological data file"<<endl;
    os<<fnSurveyData_Trawl<<tb<<"#trawl survey data file"<<endl;
    os<<fnFisheryData_TCFR<<tb<<"#directed tanner crab fishery data file"<<endl;
    os<<fnFisheryData_SnCF<<tb<<"#snow crab fishery data file"<<endl;
    os<<fnFisheryData_RKCF<<tb<<"#BB red king crab fishery data file"<<endl;
    os<<fnFisheryData_GrTF<<tb<<"#groundfish trawl fishery data file"<<endl;
    os<<(*pBio)<<endl;
    if (pTSD)  os<<(*pTSD)<<endl;
    if (pTCFR) os<<(*pTCFR)<<endl;
    if (pTCFD) os<<(*pTCFD)<<endl;
    if (pSCF) os<<(*pSCF)<<endl;
    if (pRKF) os<<(*pRKF)<<endl;
    if (pGTF) os<<(*pGTF)<<endl;
   if (debug) cout<<"end ModelDatasets::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void ModelDatasets::writeToR(ostream& os, char* nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
    os<<nm<<"=list("<<endl;
    indent++;
        //bio data
        for (int n=0;n<indent;n++) os<<tb;
        os<<"bio=list("<<endl;
        indent++;
            pBio->writeToR(os,adstring("bio"),indent); os<<endl;
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<endl;            
        
        //survey data
        for (int n=0;n<indent;n++) os<<tb;
        os<<"surveys=list("<<endl;
        indent++;
            if (pTSD) pTSD->writeToR(os,adstring("trawl"),indent); os<<endl;
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<endl;            
        
        //fishery data
        for (int n=0;n<indent;n++) os<<tb;
        os<<"fisheries=list("<<endl;
        indent++;
            if (pTCFR) pTCFR->writeToR(os,adstring("tcf.r"),indent); os<<cc<<endl;
            if (pTCFD) pTCFD->writeToR(os,adstring("tcf.d"),indent); os<<cc<<endl;
            if (pSCF) pSCF->writeToR(os,adstring("scf"),indent); os<<cc<<endl;
            if (pRKF) pRKF->writeToR(os,adstring("rkf"),indent); os<<cc<<endl;
            if (pGTF) pGTF->writeToR(os,adstring("gtf"),indent); 
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<")"<<endl;            
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
}
/////////////////////////////////end ModelDatasets/////////////////////////
