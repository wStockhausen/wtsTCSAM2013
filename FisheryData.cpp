#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "FisheryData.hpp"
//**********************************************************************
//  Includes
//      RetainedFisheryData
//      DiscardFisheryData
//      GroundfishTrawlFisheryData
//**********************************************************************
int RetainedFisheryData::debug = 0;
int DiscardFisheryData::debug  = 0;
int GroundfishTrawlFisheryData::debug  = 0;
//----------------------------------------------------------------------
//          RetainedFisheryData
//----------------------------------------------------------------------
/***************************************************************
*   read.                                                      *
***************************************************************/
void RetainedFisheryData::read(cifstream & is){
    if (debug) {
        cout<<"start RetainedFisheryData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading Directed Fishery Data."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    is>>fishery;
    if (debug) cout<<fishery<<tb<<"#fishery"<<endl;
    
    //RETAINED CATCH
    is>>nyCatch;//number of years of retained catch data
    if (debug) cout<<nyCatch<<tb<<"#nyRC"<<endl;
    is>>unitsN;
    if (debug) cout<<unitsN<<tb<<"#unitsRCN"<<endl;
    is>>unitsB;
    if (debug) cout<<unitsB<<tb<<"#unitsRCB"<<endl;
    inpCatch_yc.allocate(1,nyCatch,1,3);
    is>>inpCatch_yc;
    if (debug) cout<<"#retCatch_yc = "<<endl<<"#years  numbers  biomass"<<endl<<inpCatch_yc<<endl;
    
    yrsCatch.allocate(1,nyCatch);
    catch_ty.allocate(1,2,1,nyCatch);
    yrsCatch = (ivector) column(inpCatch_yc,1);
    for (int t=1;t<=2;t++) catch_ty(t) = column(inpCatch_yc,1+t);
    if (debug) {
        cout<<"catch_ty = "<<endl<<catch_ty<<endl;
        cout<<"Enter '1' to continue >>";
        int resp;
        cin>>resp;
    }
    
    //SIZE BINS
    is>>nZCutPts;
    if (debug) cout<<nZCutPts<<tb<<"#nZCutPts"<<endl;
    zCutPts.allocate(1,nZCutPts);
    is>>zCutPts;
    if (debug) cout<<zCutPts<<tb<<"#zCutPts (mm CW)"<<endl;
    zBins.allocate(1,nZCutPts-1);
    zBins = 0.5*(zCutPts(1,nZCutPts-1)+(--zCutPts(2,nZCutPts)));
    if (debug) cout<<zBins<<tb<<"#zBins"<<endl;
    
    //RETAINED CATCH NUMBERS-AT-SIZE 
    {is>>nyNatZ;//number of years of retained catch numbers-at-size data
    if (debug) cout<<nyNatZ<<tb<<"#nyRCNatZ"<<endl;
    is>>unitsNatZ;
    if (debug) cout<<unitsNatZ<<tb<<"#unitsRCNatZ"<<endl;
    
    yrsNatZ.allocate(1,nyNatZ);
    ssNatZ_sy.allocate(1,ALL_SHELL,1,nyNatZ);
    nAtZ_syz.allocate(1,ALL_SHELL,1,nyNatZ,1,nZCutPts-1);
    inpNatZ_syc.allocate(1,ALL_SHELL,1,nyNatZ,1,nZCutPts-1+2);
    
    yrsNatZ.initialize();
    ssNatZ_sy.initialize();
    nAtZ_syz.initialize();
    inpNatZ_syc.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    adstring_array factors(1,1);
    for (int i=0;i<nc;i++){
        is>>factors;
        int shell = ModelConsts::getShellType(factors(1));
        if (debug) cout<<factors(1)<<tb<<"#factors"<<endl;
        if (shell){
            is>>inpNatZ_syc(shell);
            if (debug) {
                cout<<"#inpNatZ_syc("<<factors(1)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                cout<<inpNatZ_syc(shell)<<endl;
            }
            for (int y=1;y<=nyNatZ;y++){
                yrsNatZ(y) = inpNatZ_syc(shell,y,1);
                ssNatZ_sy(shell,y) = inpNatZ_syc(shell,y,2);
                nAtZ_syz(shell,y) = (inpNatZ_syc(shell,y)(3,nZCutPts-1+2)).shift(1);
            }
            if (debug) {
                cout<<"#year, sample_size, nAtZ_syz("<<factors(1)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                for (int y=1;y<=nyNatZ;y++) cout<<yrsNatZ(y)<<tb<<ssNatZ_sy(shell,y)<<tb<<nAtZ_syz(shell,y)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for nAtZ_syc not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    if (debug) cout<<"end RetainedFisheryData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void RetainedFisheryData::write(std::ostream & os){
    if (debug) cout<<"start RetainedFisheryData::write(...) "<<this<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<fishery<<endl;
    
    os<<"#-----------RETAINED CATCH (NUMBERS/BIOMASS)---------------#"<<endl;
    os<<nyCatch<<tb<<"#number of years of retained catch data"<<endl;
    os<<unitsN<<tb<<"#units for numbers retained"<<endl;
    os<<unitsB<<tb<<"#units for biomass retained"<<endl;
    os<<"#year numbers biomass"<<endl<<inpCatch_yc<<endl;
    
    os<<"#-----------SIZE BINS---------------#"<<endl;
    os<<nZCutPts<<tb<<"#number of size bin cutpoints"<<endl;
    os<<"#size bin cutpoints (mm CW)"<<endl<<zCutPts<<endl;
    
    {os<<"#-----------RETAINED CATCH NUMBERS-AT-SIZE DATA----------#"<<endl;
    os<<nyNatZ<<tb<<"#number of years of retained catch numbers-at-size data"<<endl;
    os<<unitsNatZ<<tb<<"#units for numbers at retained catch numbers-at-size"<<endl;
    os<<ALL_SHELL<<tb<<"#number of shell factor combinations to read in"<<endl;
    adstring_array factors(1,1);
    for (int shell=NEW_SHELL;shell<=ALL_SHELL;shell++){
        factors(1) = ModelConsts::getShellType(shell);
        os<<"#-------"<<factors(1)<<endl;
        os<<factors(1)<<endl;
        os<<"#year sample_size "<<zBins<<endl;
        os<<inpNatZ_syc(shell)<<endl;
    }}
    
    if (debug) cout<<"end RetainedFisheryData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void RetainedFisheryData::writeToR(std::ostream& os, char* nm, int indent) {
    if (debug) cout<<"RetainedFisheryData::writing to R"<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list(name="<<qt<<fishery<<qt<<cc<<endl;
        //retained catch
        indent++;
        {   for (int n=0;n<indent;n++) os<<tb;
            adstring y  = wts::to_qcsv(yrsCatch);
            os<<"catch=list("<<endl;
            indent++; 
                for (int n=0;n<indent;n++) os<<tb;
                os<<"years="; wts::writeToR(os,yrsCatch); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"numbers=list(units="<<qt<<unitsN<<qt<<cc<<"data="; wts::writeToR(os,catch_ty(1),y); os<<"),"<<endl;
                for (int n=0;n<indent;n++) os<<tb; 
                os<<"biomass=list(units="<<qt<<unitsB<<qt<<cc<<"data="; wts::writeToR(os,catch_ty(2),y); os<<")"<<endl;
            indent--;
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        indent--;}
        
        //retained catch NatZ
        {indent++;   
            for (int n=0;n<indent;n++) os<<tb;
            adstring s = "s=c("+qt+STR_NEW_SHELL+qt +cc+ qt+STR_OLD_SHELL+qt +cc+ qt+STR_ALL_SHELL+qt+")";
            adstring y = "y=c("+wts::to_qcsv(yrsNatZ)+")";
            adstring z = "z=c("+wts::to_qcsv(zBins)+")";        
            os<<"nAtZ=list(units="<<qt<<unitsNatZ<<qt<<cc<<endl; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"years="; wts::writeToR(os,yrsNatZ); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"cutpts="; wts::writeToR(os,zCutPts); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"sample.sizes="; wts::writeToR(os,ssNatZ_sy,s,y); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<endl;
            wts::writeToR(os,nAtZ_syz,s,y,z); os<<endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<")"<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
    if (debug) cout<<"RetainedFisheryData::done writing to R"<<endl;
}
/////////////////////////////////end RetainedFisheryData/////////////////////////
//----------------------------------------------------------------------
//          DiscardFisheryData
//----------------------------------------------------------------------
/***************************************************************
*   read.                                                      *
***************************************************************/
void DiscardFisheryData::read(cifstream & is){
    if (debug) {
        cout<<"start DiscardFisheryData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading DiscardFisheryData."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    is>>fishery;
    if (debug) cout<<fishery<<tb<<"#fishery"<<endl;
    
    //DISCARDED CATCH
    is>>nyCatch;//number of years of discarded catch data
    if (debug) cout<<nyCatch<<tb<<"#nyCatch"<<endl;
    is>>unitsCatch;
    if (debug) cout<<unitsCatch<<tb<<"#unitsCatch"<<endl;
    inpCatch_yc.allocate(1,nyCatch,1,3);
    is>>inpCatch_yc;
    if (debug) cout<<"#dscCatch_yc = "<<endl<<"#year  females  male"<<endl<<inpCatch_yc<<endl;
    
    yrsCatch.allocate(1,nyCatch);
    catch_xy.allocate(1,nSEXS,1,nyCatch);
    yrsCatch = (ivector) column(inpCatch_yc,1);
    for (int x=FEMALE;x<=MALE;x++) catch_xy(x) = column(inpCatch_yc,1+x);
    
    //EFFORT DATA
    is>>nyEff;//number of years of effort data
    if (debug) cout<<nyEff<<tb<<"#nyEff"<<endl;
    is>>unitsPLs;
    if (debug) cout<<unitsPLs<<tb<<"#unitsPLs"<<endl;
    effort_yc.allocate(1,nyEff,1,2);
    is>>effort_yc;
    if (debug) cout<<"#effort_yc = "<<endl<<"#year potlifts"<<endl<<effort_yc<<endl;
    
    yrsEffort.allocate(1,nyEff);
    effort_y.allocate(1,nyEff);
    yrsEffort = (ivector) column(effort_yc,1);
    effort_y  = column(effort_yc,2);
    
    //SIZE BINS
    is>>nZCutPts;
    if (debug) cout<<nZCutPts<<tb<<"#nZCutPts"<<endl;
    zCutPts.allocate(1,nZCutPts);
    is>>zCutPts;
    if (debug) cout<<zCutPts<<tb<<"#zCutPts (mm CW)"<<endl;
    zBins.allocate(1,nZCutPts-1);
    zBins = 0.5*(zCutPts(1,nZCutPts-1)+(--zCutPts(2,nZCutPts)));
    if (debug) cout<<zBins<<tb<<"#zBins"<<endl;
    
    //CATCH NUMBERS-AT-SIZE 
    {is>>nyNatZ;//number of years of retained catch numbers-at-size data
    if (debug) cout<<nyNatZ<<tb<<"#nyNatZ"<<endl;
    is>>unitsNatZ;
    if (debug) cout<<unitsNatZ<<tb<<"#unitsNatZ"<<endl;
    
    yrsNatZ.allocate(1,nyNatZ);
    ssNatZ_xsy.allocate(1,ALL_SEXES,1,ALL_SHELL,1,nyNatZ);
    nAtZ_xsyz.allocate(1,ALL_SEXES,1,ALL_SHELL,1,nyNatZ,1,nZCutPts-1);
    inpNatZ_xsyc.allocate(1,ALL_SEXES,1,ALL_SHELL,1,nyNatZ,1,nZCutPts-1+2);
    
    yrsNatZ.initialize();
    ssNatZ_xsy.initialize();
    nAtZ_xsyz.initialize();
    inpNatZ_xsyc.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = ModelConsts::getSexType(factors(1));
        int shell = ModelConsts::getShellType(factors(2));
        if (debug) cout<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<endl;
        if (shell){
            is>>inpNatZ_xsyc(sex,shell);
            if (debug) {
                cout<<"#inpNatZ_xsyc("<<factors(1)<<cc<<factors(2)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                cout<<inpNatZ_xsyc(sex,shell)<<endl;
            }
            yrsNatZ = (ivector) column(inpNatZ_xsyc(sex,shell),1);
            ssNatZ_xsy(sex,shell) = column(inpNatZ_xsyc(sex,shell),2);
            for (int y=1;y<=nyNatZ;y++){
                nAtZ_xsyz(sex,shell,y)  = (inpNatZ_xsyc(sex,shell,y)(3,nZCutPts-1+2)).shift(1);
            }
            if (debug) {
                cout<<"#year, sample_size, nAtZ_xsyz("<<factors(1)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                for (int y=1;y<=nyNatZ;y++) cout<<yrsNatZ(y)<<tb<<ssNatZ_xsy(sex,shell,y)<<tb<<nAtZ_xsyz(sex,shell,y)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for nAtZ_xsyc not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    if (debug) cout<<"end DiscardFisheryData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void DiscardFisheryData::write(std::ostream & os){
    if (debug) cout<<"start DiscardFisheryData::write(...) "<<this<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<fishery<<endl;
    
    os<<"#-----------CATCH (BIOMASS)----------------------#"<<endl;
    os<<nyCatch<<tb<<"#number of years of discarded catch data"<<endl;
    os<<unitsCatch<<tb<<"#units for biomass"<<endl;
    os<<"#year  females  males"<<endl<<inpCatch_yc<<endl;
    
    os<<"#-----------EFFORT---------------#"<<endl;
    os<<nyEff<<tb<<"#number of years of effort data"<<endl;
    os<<unitsPLs<<tb<<"#units for pot lifts"<<endl;
    os<<"#year   potlifts"<<endl<<effort_yc<<endl;

    os<<"#-----------SIZE BINS---------------#"<<endl;
    os<<nZCutPts<<tb<<"#number of size bin cutpoints"<<endl;
    os<<"#size bin cutpoints (mm CW)"<<endl<<zCutPts<<endl;
    
    {os<<"#-----------CATCH NUMBERS-AT-SIZE DATA----------#"<<endl;
    os<<nyNatZ<<tb<<"#number of years of catch numbers-at-size data"<<endl;
    os<<unitsNatZ<<tb<<"#units for catch numbers-at-size"<<endl;
    os<<ALL_SEXES*ALL_SHELL<<tb<<"#number of factor combinations to read in"<<endl;
    adstring_array factors(1,2);
    for (int sex=FEMALE;sex<=ALL_SEXES;sex++){
        for (int shell=NEW_SHELL;shell<=ALL_SHELL;shell++){
            factors(1) = ModelConsts::getSexType(sex);
            factors(2) = ModelConsts::getShellType(shell);
            os<<"#-----"<<endl;
            os<<factors(1)<<tb<<factors(2)<<endl;
            os<<"#year sample_size "<<zBins<<endl;
            os<<inpNatZ_xsyc(sex,shell)<<endl;
        }
    }}
    if (debug) cout<<"end DiscardFisheryData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void DiscardFisheryData::writeToR(std::ostream& os, char* nm, int indent) {
    if (debug) cout<<"DiscardFisheryData::writing to R"<<endl;
    adstring x = "x=c("+qt+STR_FEMALE+qt+cc+qt+STR_MALE+qt+")";
    adstring m = "m=c("+qt+STR_IMMATURE+qt +cc+ qt+STR_MATURE+qt+")";
    adstring s = "s=c("+qt+STR_IMMATURE+qt+cc+qt+STR_MATURE+qt+")";
    adstring z = "z=c("+wts::to_qcsv(zBins)+")";
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list(name="<<qt<<fishery<<qt<<cc<<endl;
        //discarded catch
        indent++;
        {   for (int n=0;n<indent;n++) os<<tb;
            adstring y  = "y=c("+wts::to_qcsv(yrsCatch)+")";
            os<<"catch=list("<<endl;
            indent++; 
                for (int n=0;n<indent;n++) os<<tb;
                os<<"years="; wts::writeToR(os,yrsCatch); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"biomass=list(units="<<qt<<unitsCatch<<qt<<cc<<"data="; wts::writeToR(os,catch_xy,x,y); os<<")"<<endl;
            indent--;
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        indent--;}
        
        //effort
        indent++;
        {   for (int n=0;n<indent;n++) os<<tb;
            adstring y  = "y=c("+wts::to_qcsv(yrsCatch)+")";
            os<<"effort=list("<<endl;
            indent++; 
                for (int n=0;n<indent;n++) os<<tb;
                os<<"units="<<qt<<unitsPLs<<qt<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"years="; wts::writeToR(os,yrsEffort); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"data="; wts::writeToR(os,effort_y,y); os<<endl;
            indent--;
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        indent--;}
        
        //discarded catch nAtZ
        {indent++;   
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = "x=c("+qt+STR_FEMALE+qt +cc+ qt+STR_MALE+qt +cc+ qt+STR_ALL_SEXES+qt+")";
            adstring s = "s=c("+qt+STR_NEW_SHELL+qt +cc+ qt+STR_OLD_SHELL+qt +cc+ qt+STR_ALL_SHELL+qt+")";
            adstring y = "y=c("+wts::to_qcsv(yrsNatZ)+")";
            os<<"nAtZ=list(units="<<qt<<unitsNatZ<<qt<<cc<<endl; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"years="; wts::writeToR(os,yrsNatZ); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"cutpts="; wts::writeToR(os,zCutPts); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"sample.sizes="; wts::writeToR(os,ssNatZ_xsy,x,s,y); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<endl;
            wts::writeToR(os,nAtZ_xsyz,x,s,y,z); os<<endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<")"<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
    if (debug) cout<<"DiscardFisheryData::done writing to R"<<endl;
}
/////////////////////////////////end DiscardFisheryData/////////////////////////
//----------------------------------------------------------------------
//          GroundfishTrawlFisheryData
//----------------------------------------------------------------------
/***************************************************************
*   read.                                                      *
***************************************************************/
void GroundfishTrawlFisheryData::read(cifstream & is){
    if (debug) {
        cout<<"start GroundfishTrawlFisheryData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading GroundfishTrawlFisheryData."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    is>>fishery;
    if (debug) cout<<fishery<<tb<<"#fishery"<<endl;
    
    //DISCARDED CATCH
    is>>nyCatch;//number of years of discarded catch data
    if (debug) cout<<nyCatch<<tb<<"#nyCatch"<<endl;
    is>>unitsCatch;
    if (debug) cout<<unitsCatch<<tb<<"#unitsCatch"<<endl;
    inpCatch_yc.allocate(1,nyCatch,1,2);//year, total discard biomass
    is>>inpCatch_yc;
    if (debug) cout<<"#dscCatch_yc = "<<endl<<"#year  total biomass"<<endl<<inpCatch_yc<<endl;
    
    yrsCatch.allocate(1,nyCatch);
    catch_y.allocate(1,nyCatch);
    yrsCatch = (ivector) column(inpCatch_yc,1);
    catch_y  = column(inpCatch_yc,2);//total discard biomass
    
    //SIZE BINS
    is>>nZCutPts;
    if (debug) cout<<nZCutPts<<tb<<"#nZCutPts"<<endl;
    zCutPts.allocate(1,nZCutPts);
    is>>zCutPts;
    if (debug) cout<<zCutPts<<tb<<"#zCutPts (mm CW)"<<endl;
    zBins.allocate(1,nZCutPts-1);
    zBins = 0.5*(zCutPts(1,nZCutPts-1)+(--zCutPts(2,nZCutPts)));
    if (debug) cout<<zBins<<tb<<"#zBins"<<endl;
    
    //CATCH NUMBERS-AT-SIZE 
    {is>>nyNatZ;//number of years of retained catch numbers-at-size data
    if (debug) cout<<nyNatZ<<tb<<"#nyNatZ"<<endl;
    is>>unitsNatZ;
    if (debug) cout<<unitsNatZ<<tb<<"#unitsNatZ"<<endl;
    
    yrsNatZ.allocate(1,nyNatZ);
    ssNatZ_xsy.allocate(1,ALL_SEXES,1,ALL_SHELL,1,nyNatZ);
    nAtZ_xsyz.allocate(1,ALL_SEXES,1,ALL_SHELL,1,nyNatZ,1,nZCutPts-1);
    inpNatZ_xsyc.allocate(1,ALL_SEXES,1,ALL_SHELL,1,nyNatZ,1,nZCutPts-1+2);
    
    yrsNatZ.initialize();
    ssNatZ_xsy.initialize();
    nAtZ_xsyz.initialize();
    inpNatZ_xsyc.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = ModelConsts::getSexType(factors(1));
        int shell = ModelConsts::getShellType(factors(2));
        if (debug) cout<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<endl;
        if (shell){
            is>>inpNatZ_xsyc(sex,shell);
            if (debug) {
                cout<<"#inpNatZ_xsyc("<<factors(1)<<cc<<factors(2)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                cout<<inpNatZ_xsyc(sex,shell)<<endl;
            }
            for (int y=1;y<=nyNatZ;y++){
                yrsNatZ(y) = inpNatZ_xsyc(sex,shell,y,1);
                ssNatZ_xsy(sex,shell,y) = inpNatZ_xsyc(sex,shell,y,2);
                nAtZ_xsyz(sex,shell,y) = (inpNatZ_xsyc(sex,shell,y)(3,nZCutPts-1+2)).shift(1);
            }
            if (debug) {
                cout<<"#year, sample_size, nAtZ_xsyz("<<factors(1)<<")="<<endl;
                cout<<"#year sample_size "<<zBins<<endl;
                for (int y=1;y<=nyNatZ;y++) cout<<yrsNatZ(y)<<tb<<ssNatZ_xsy(sex,shell,y)<<tb<<nAtZ_xsyz(sex,shell,y)<<endl;
            }
        } else {
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            cout<<"Reading file name "<<is.get_file_name()<<endl;
            cout<<"Factors for nAtZ_xsyc not recognized!"<<endl;
            cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<endl;
            cout<<"Terminating..."<<endl;
            cout<<"-----------------------------------------------------------"<<endl;
            exit(-1);
        }
    }}
    if (debug) cout<<"end GroundfishTrawlFisheryData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void GroundfishTrawlFisheryData::write(std::ostream & os){
    if (debug) cout<<"start GroundfishTrawlFisheryData::write(...) "<<this<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<"#-----------------------------------------------------#"<<endl;
    os<<fishery<<endl;
    
    os<<"#-----------CATCH (BIOMASS)----------------------#"<<endl;
    os<<nyCatch<<tb<<"#number of years of discarded catch data"<<endl;
    os<<unitsCatch<<tb<<"#units for biomass"<<endl;
    os<<"#year  total biomass"<<endl<<inpCatch_yc<<endl;
    
    os<<"#-----------SIZE BINS---------------#"<<endl;
    os<<nZCutPts<<tb<<"#number of size bin cutpoints"<<endl;
    os<<"#size bin cutpoints (mm CW)"<<endl<<zCutPts<<endl;
    
    {os<<"#-----------CATCH NUMBERS-AT-SIZE DATA----------#"<<endl;
    os<<nyNatZ<<tb<<"#number of years of catch numbers-at-size data"<<endl;
    os<<unitsNatZ<<tb<<"#units for catch numbers-at-size"<<endl;
    os<<MALE*ALL_SHELL<<tb<<"#number of factor combinations to read in"<<endl;
    adstring_array factors(1,2);
    for (int sex=FEMALE;sex<=MALE;sex++){
        for (int shell=NEW_SHELL;shell<=ALL_SHELL;shell++){
            factors(1) = ModelConsts::getSexType(sex);
            factors(2) = ModelConsts::getShellType(shell);
            os<<"#-----"<<endl;
            os<<factors(1)<<tb<<factors(2)<<endl;
            os<<"#year sample_size "<<zBins<<endl;
            os<<inpNatZ_xsyc(sex,shell)<<endl;
        }
    }}
    if (debug) cout<<"end GroundfishTrawlFisheryData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void GroundfishTrawlFisheryData::writeToR(std::ostream& os, char* nm, int indent) {
    if (debug) cout<<"GroundfishTrawlFisheryData::writing to R"<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list(name="<<qt<<fishery<<qt<<cc<<endl;
        //discarded catch
        indent++;
        {   for (int n=0;n<indent;n++) os<<tb;
            adstring y  = wts::to_qcsv(yrsCatch);
            os<<"catch=list("<<endl;
            indent++; 
                for (int n=0;n<indent;n++) os<<tb;
                os<<"years="; wts::writeToR(os,yrsCatch); os<<cc<<endl;
                for (int n=0;n<indent;n++) os<<tb;
                os<<"biomass=list(units="<<qt<<unitsCatch<<qt<<cc<<"data="; wts::writeToR(os,catch_y,y); os<<")"<<endl;
            indent--;
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<endl;
        indent--;}
        
        //discarded catch nAtZ
        {indent++;   
            for (int n=0;n<indent;n++) os<<tb;
            adstring x = "x=c("+qt+STR_FEMALE+qt +cc+ qt+STR_MALE+qt +cc+ qt+STR_ALL_SEXES+qt+")";
            adstring s = "s=c("+qt+STR_NEW_SHELL+qt +cc+ qt+STR_OLD_SHELL+qt +cc+ qt+STR_ALL_SHELL+qt+")";
            adstring y = "y=c("+wts::to_qcsv(yrsNatZ)+")";
            adstring z = "z=c("+wts::to_qcsv(zBins)+")";        
            os<<"nAtZ=list(units="<<qt<<unitsNatZ<<qt<<cc<<endl; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"years="; wts::writeToR(os,yrsNatZ); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"cutpts="; wts::writeToR(os,zCutPts); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"sample.sizes="; wts::writeToR(os,ssNatZ_xsy,x,s,y); os<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<endl;
            wts::writeToR(os,nAtZ_xsyz,x,s,y,z); os<<endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<")"<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
    if (debug) cout<<"GroundfishTrawlFisheryData::done writing to R"<<endl;
}
/////////////////////////////////end GroundfishTrawlFisheryData/////////////////////////

