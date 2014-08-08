#ifndef EUTelFitTuple_h
#define EUTelFitTuple_h 1

#include "marlin/Processor.h"

#include "EUTelAlignmentConstant.h"
#include "EUTelReferenceHit.h"
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>


// system includes <>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TVectorT.h>
#include <gsl/gsl_matrix_double.h>

namespace eutelescope {
  class EUTelAPIXTbTrackTuple : public marlin::Processor {
    
  public:
    virtual Processor*  newProcessor() { return new EUTelAPIXTbTrackTuple; }
    
    EUTelAPIXTbTrackTuple() ;
    virtual void init() ;
    virtual void processRunHeader( LCRunHeader* run ) ;
    virtual void processEvent( LCEvent * evt ) ;
    virtual void check( LCEvent * /*evt*/ ){;} ;
    virtual void end() ;
    
    
    void readAlignment(LCEvent * event);
    void invertAlignment(EUTelAlignmentConstant * alignment);
    void invertGear();
    void reverseAlign(double& x, double& y, double& z, int iden, double nomZPos);
    gsl_matrix* invertLU(int dim, gsl_matrix* m);
    void getDUTRot(EUTelAlignmentConstant * alignment);


 
    
       
  protected:
    //TbTrack additions
    void prepareTree();
    void clear();

    bool readZsHits(std::string colName, LCEvent* event);
    bool readTracks(LCEvent* event);
    bool readHits( std::string hitColName, LCEvent* event );

    std::string _inputTrackColName;
    std::string _inputTrackerHitColName;
    std::string _inputTelPulseCollectionName;
    std::string _inputDutPulseCollectionName;
    std::string _telZsColName;
    std::string _dutZsColName;

    std::string _path2file;

    std::vector<int> _DUTIDs;
 
    // Internal processor variables
    // ----------------------------
    int _nRun ;
    int _nEvt ;
    int _runNr;
    int _evtNr;

    bool _isFirstEvent;
    
    TFile* _file;
    
    TTree* _eutracks;
    
    int _nTrackParams;
    std::vector<double> *_xPos;
    std::vector<double> *_yPos;
    std::vector<double> *_dxdz;
    std::vector<double> *_dydz;
    std::vector<int>    *_trackIden;
    std::vector<int>    *_trackNum;
    std::vector<double> *_chi2;
    std::vector<double> *_ndof;    

    TTree* _zstree;
    int _nPixHits;
    std::vector<int> *p_col;
    std::vector<int> *p_row;
    std::vector<int> *p_tot;
    std::vector<int> *p_iden;
    std::vector<int> *p_lv1;
    std::vector<int> *p_chip;

    TTree* _euhits;
    int _nHits;
    std::vector<double> *_hitXPos;
    std::vector<double> *_hitYPos;
    std::vector<double> *_hitZPos;
    std::vector<int>    *_hitSensorId;
    
    TTree* _rottree;
    std::vector<int> *_rotDUTId;  
    std::vector<double> *_alpha;
    std::vector<double> *_beta; 
    std::vector<double> *_gamma;  
    std::vector<double> *_rotZY;
    std::vector<double> *_rotZX; 
    std::vector<double> *_rotXY;
    std::vector<double> *_rotZYerr;
    std::vector<double> *_rotZXerr; 
    std::vector<double> *_rotXYerr;     
    
    TVectorD* _versionVec;
    
    
    
    bool _foundAllign;
    bool _doScales;
    std::vector<std::string> _alignColNames;
    std::map<int, bool> _rotationstored;
    int _countrotstored;
    std::string _referenceHitCollectionName;
    LCCollectionVec* _referenceHitVec; 
    
    std::map<int, std::vector< std::vector<double> > > _alignShift;
    std::map<int, std::vector< std::vector<double> > > _alignRotations;
    //std::map< int, std::vector< std::vector<double> > >
    std::map<int, std::vector<gsl_matrix*> > _alignRot;
    std::map<int, std::vector<double> > _gearShift;
    std::map<int, std::vector<double> > _gearOffset;
    std::map<int, std::vector<double> > _gearSize;
    std::map<int, gsl_matrix* > _gearRot;
    std::map<int, std::pair<double, double> > _gearPitch;
    std::map<int, std::vector<double> > _gearEulerRot;

    gear::SiPlanesParameters * _siPlanesParameters;
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;
  };

  //! A global instance of the processor.
  EUTelAPIXTbTrackTuple aEUTelAPIXTbTrackTuple;
}
#endif
