/* 
 * File:   ModelData.hpp
 * Author: william.stockhausen
 *
 * Created on May 14, 2013, 6:52 AM
 */

#ifndef MODELDATA_HPP
#define MODELDATA_HPP

//**********************************************************************
//  Includes
//      BioData
//      TrawlSurveyData
//      ModelDatasets
//**********************************************************************
class ModelConfiguration; //forward definition
class RetainedFisheryData;
class DiscardFisheryData;
class GroundfishTrawlFisheryData;
//--------------------------------------------------------------------------------
//          BioData
//  
//--------------------------------------------------------------------------------
    class BioData {
    public:
        static int debug;
    private:
    public:
        int nZBins;           //number of size bin cut pts
        dvector zBins;        //size bins
        adstring unitsWatZ;   //units for weight-at-size
        d3_array wAtZ_xmz;    //weight at size
        dmatrix prMature_xz;  //probability of maturity at size by sex
        d3_array frMature_xsz;//fraction mature at size by sex, shell condition
        dmatrix cvMnMxZ_xc;    //cv of min size, max size by sex
        
        int nyFshSeasons;     //number of years for fishery season midpoints
        dmatrix mdptFshSeasons_yc;//midpoint of fishery season by year
    public:
        BioData(){}
        ~BioData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, char* nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, BioData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   BioData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          TrawlSurveyData
//  
//--------------------------------------------------------------------------------
    class TrawlSurveyData {
    public:
        static int debug;
    private:
        dmatrix inpAbund_yc;     //input abundance data (year,year+female abundance+male abundance+cv_females+cv_males)
        dmatrix inpMatBio_yc;    //input mature biomass data (year,year+female biomass+male biomass+cv_females+cv_males)
        d5_array inpNatZ_xsmyc;  //input numbers-at-size data (sex,shell,maturity,year,year+sample_size+nAtZ)
    public:
        int nyAbund;          //number of years of abundance data
        adstring unitsAbund;  //units for abundance data
        ivector yrsAbund;     //years for abundance data
        dvector abund_y;      //abundance by year
        dmatrix abund_xy;     //abundance by sex, year
        dmatrix cvsAbund_xy;  //cvs by sex, year
        
        int nyMatBio;          //number of years of mature biomass data
        adstring unitsMatBio;  //units for mature biomass data
        ivector yrsMatBio;     //years for mature biomass data
        dmatrix matBio_xy;     //mature biomass by sex, year
        dmatrix cvsMatBio_xy;  //cvs by sex, year
        
        int nZCutPts;         //number of size bin cut pts
        dvector zCutPts;      //cut points for size bins
        dvector zBins;        //size bins
        
        int nyNatZ;           //number of years of numbers-at-size data
        adstring unitsNatZ;   //units for numbers-at-size data
        ivector  yrsNatZ;     //years of size frequency data
        d4_array ssNatZ_xsmy; //sample sizes for size frequency data
        d5_array nAtZ_xsmyz;  //input numbers-at-size data (sex,shell,maturity,year,size)
    public:
        TrawlSurveyData(){}
        ~TrawlSurveyData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, char* nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, TrawlSurveyData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   TrawlSurveyData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//         ModelDatasets
//  
//--------------------------------------------------------------------------------
    class ModelDatasets {
    public:
        static int debug;
    public:
        ModelConfiguration* pMC;//pointer to ModelConfiguration object
        adstring fnBioData;        //bio data
        adstring fnSurveyData_Trawl;//trawl survey data
        adstring fnFisheryData_TCFR;//directed tanner crab fishery data (retained)
        adstring fnFisheryData_TCFD;//directed tanner crab fishery data (discarded)
        adstring fnFisheryData_SnCF;//snow crab fishery bycatch data
        adstring fnFisheryData_RKCF;//BB red king crab bycatch data
        adstring fnFisheryData_GrTF;//groundfish trawl bycatch data
        
        BioData*               pBio;//pointer to bio dataset object
        TrawlSurveyData*       pTSD;//pointer to trawl survey dataset object
        RetainedFisheryData*   pTCFR;//pointer to directed tanner crab fishery retained catch data
        DiscardFisheryData*    pTCFD;//pointer to directed tanner crab fishery discarded catch data
        DiscardFisheryData*    pSCF;//pointer to snow crab fishery discard catch data
        DiscardFisheryData*    pRKF;//pointer to BB red king crab discard catch data
        GroundfishTrawlFisheryData* pGTF;//pointer to groudnfish trawl fishery discard data
    public:
        ModelDatasets(ModelConfiguration* ptrMC);
        ~ModelDatasets();
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, char* nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, ModelDatasets & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   ModelDatasets & obj){obj.write(os); return os;}
    };
#endif  /* MODELDATA_HPP */

