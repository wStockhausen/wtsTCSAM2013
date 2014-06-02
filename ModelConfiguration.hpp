/* 
 * File:   ModelConfiguration.hpp
 * Author: william.stockhausen
 *
 * Created on March 12, 2013, 8:59 AM
 */

#ifndef MODELCONFIGURATION_HPP
    #define	MODELCONFIGURATION_HPP

//--------------------------------------------------------------------------------
//          ModelConfiguration
//--------------------------------------------------------------------------------
    class ModelConfiguration {
    public:
        static int debug;  //flag to print debug info
    public:
        adstring cfgName;//model configuration name
        int asmtYr; //assessment year
        int mnYr;  //min model year
        int mxYr;  //max model year
        int nZBins;//number of size bins 
        
        double jitFrac;//jitter fraction

        int inclTSD;  //flag to include trawl survey data
        int inclTCF;  //flag to include directed tanner fishery
        int inclSCF;  //flag to include opilio bycatch fishery
        int inclRKF;  //flag to include red king crab bycatch fishery
        int inclGTF;  //flag to include groundfish trawl bycatch fishery
        int runOpMod;    //flag to run operating model
        int fitToPriors; //flag to fit to priors
        
        dvector zBins;     //size bin midpoints (CW in mm)
        dvector zBinCutPts;//size bin cutpoints (CW in mm)
        dvector onesZBins; //vector of 1's

        adstring fnMPC; //model parameters configuration file name
        adstring fnMDS; //model datasets file name

    public:
        ModelConfiguration();
        ~ModelConfiguration();
        ModelConfiguration& operator =(const ModelConfiguration & rhs);
        
        int isModelYear(int yr){if ((mnYr<=yr)&&(yr<=mxYr)) return 1; return 0;}
        void read(const adstring & fn);   //read file in ADMB format
        void write(const adstring & fn);  //write object to file in ADMB format
        void read(cifstream & is);        //read file in ADMB format
        void write(std::ostream & os);         //write object to file in ADMB format
        void writeToR(std::ostream& os, char* nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, ModelConfiguration & obj){obj.read(is);return is;}
        friend std::ostream&   operator <<(std::ostream & os,   ModelConfiguration & obj){obj.write(os);;return os;}
    };
    
    
#endif	/* MODELCONFIGURATION_HPP */

