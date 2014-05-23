/* 
 * File:   FisheryData.hpp
 * Author: william.stockhausen
 *
 * Created on May 17, 2013, 3:29 PM
 */

#ifndef FISHERYDATA_HPP
#define	FISHERYDATA_HPP

//**********************************************************************
//  Includes
//      RetainedFisheryData
//      DiscardFisheryData
//**********************************************************************
//--------------------------------------------------------------------------------
//          RetainedFisheryData
//  
//--------------------------------------------------------------------------------
    class RetainedFisheryData {
    public:
        static int debug;
    private:
        dmatrix inpCatch_yc;//input retained catch data (year,year+ret catch numbers+ret catch biomass)
        d3_array inpNatZ_syc; //input retained catch (males only) numbers-at-size data (shell,year,year+sample_size+nAtZ)
    public:
        adstring fishery;  //fishery name
        
        int nyCatch;          //number of years of retained catch data
        adstring unitsN; //units for retained catch numbers
        adstring unitsB; //units for retained catch biomass
        ivector yrsCatch;  //years of catch data
        dmatrix catch_ty;   //retained catch data (year,number + biomass)
        
        int nZCutPts;         //number of size bin cut pts
        dvector zCutPts;      //cut points for size bins
        dvector zBins;        //size bins
        
        int nyNatZ;           //number of years of retained catch (males only) numbers-at-size data
        adstring unitsNatZ;   //units for retained catch (males only) numbers-at-size data
        ivector  yrsNatZ;     //years of size frequency data
        dmatrix  ssNatZ_sy;   //sample sizes for size frequency data
        d3_array nAtZ_syz;    //size frequency data
    public:
        RetainedFisheryData(){}
        ~RetainedFisheryData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, char* nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, RetainedFisheryData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   RetainedFisheryData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          DiscardFisheryData
//  
//--------------------------------------------------------------------------------
    class DiscardFisheryData {
    public:
        static int debug;
    private:
        dmatrix inpCatch_yc;  //input catch data (year,year+female biomass+male biomass)
        d4_array inpNatZ_xsyc;//input numbers-at-size data (sex,shell,year,year+sample_size+nAtZ)
    public:
        adstring fishery;     //fishery name
        
        int nyCatch;          //number of years of discarded catch data
        adstring unitsCatch;  //units for discarded catch biomass
        ivector yrsCatch;    //years of catch data
        dmatrix catch_xy;    //catch data (year, female + male biomass)
        
        int nyEff;           //number of years of effort data
        adstring unitsPLs;   //units for potlifts
        dmatrix  effort_yc;  //input effort data (year,year+potlifts)
        dvector yrsEffort;   //years w/ effort data
        dvector effort_y;    //effort data
        
        int nZCutPts;         //number of size bin cut pts
        dvector zCutPts;      //cut points for size bins
        dvector zBins;        //size bins
        
        int nyNatZ;           //number of years of numbers-at-size data
        adstring unitsNatZ;   //units for numbers-at-size data
        ivector  yrsNatZ;     //years of size frequency data
        d3_array ssNatZ_xsy;  //sample sizes for size frequency data
        d4_array nAtZ_xsyz;   //size frequency data
    public:
        DiscardFisheryData(){}
        ~DiscardFisheryData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, char* nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, DiscardFisheryData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   DiscardFisheryData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          GroundfishTrawlFisheryData
//  
//--------------------------------------------------------------------------------
    class GroundfishTrawlFisheryData {
    public:
        static int debug;
    private:
        dmatrix inpCatch_yc;//input catch data (year,year+total biomass)
        d4_array inpNatZ_xsyc; //input numbers-at-size data (sex,shell,year,year+sample_size+nAtZ)
    public:
        adstring fishery;     //fishery name
        
        int nyCatch;        //number of years of discarded catch data
        adstring unitsCatch;//units for discarded catch biomass
        ivector yrsCatch;   //years of catch data
        dvector catch_y;   //catch data (total discard biomass by year)
        
        int nZCutPts;         //number of size bin cut pts
        dvector zCutPts;      //cut points for size bins
        dvector zBins;        //size bins
        
        int nyNatZ;            //number of years of numbers-at-size data
        adstring unitsNatZ;    //units for numbers-at-size data
        ivector  yrsNatZ;     //years of size frequency data
        d3_array ssNatZ_xsy;  //sample sizes for size frequency data
        d4_array nAtZ_xsyz;   //size frequency data
    public:
        GroundfishTrawlFisheryData(){}
        ~GroundfishTrawlFisheryData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(ostream & os); //write object to file in ADMB format
        void writeToR(ostream& os, char* nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, GroundfishTrawlFisheryData & obj){obj.read(is); return is;}
        friend ostream&   operator <<(ostream & os,   GroundfishTrawlFisheryData & obj){obj.write(os); return os;}
    };

#endif	/* FISHERYDATA_HPP */

