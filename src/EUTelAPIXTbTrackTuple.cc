// eutelescope inlcudes
#include "EUTelAPIXTbTrackTuple.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDDecoder.h>


//TbTrack include
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <TVector3.h>
#include <algorithm>

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

EUTelAPIXTbTrackTuple::EUTelAPIXTbTrackTuple()
: Processor("EUTelAPIXTbTrackTuple"),
  _inputTrackColName (""),
  _inputTrackerHitColName (""),
  _inputTelPulseCollectionName(""),
  _inputDutPulseCollectionName(""),
  _telZsColName(""),
  _dutZsColName(""),
  _path2file(""),
  _DUTIDs(std::vector<int>()),
  _nRun (0),
  _nEvt (0),
  _runNr(0),
  _evtNr(0),
  _isFirstEvent(false),
  _file(NULL),
  _eutracks(NULL),
  _nTrackParams(0),
  _xPos(NULL),
  _yPos(NULL),
  _dxdz(NULL),
  _dydz(NULL),
  _trackIden(NULL),
  _trackNum(NULL),
  _chi2(NULL),
  _ndof(NULL),
  _zstree(NULL),
  _nPixHits(0),
  p_col(NULL),
  p_row(NULL),
  p_tot(NULL),
  p_iden(NULL),
  p_lv1(NULL),
  _euhits(NULL),
  _nHits(0),
  _hitXPos(NULL),
  _hitYPos(NULL),
  _hitZPos(NULL),
  _hitSensorId(NULL)
 {
  //processor description
  _description = "Prepare tbtrack style n-tuple with track fit results" ;


  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input Track collection",
		  		_inputTrackColName, std::string("fittracks"));

  	   
  registerInputCollection( LCIO::TRACKERHIT, "InputTrackerHitCollectionName", "Name of the plane-wide hit-data hit collection"  ,
		_inputTrackerHitColName, std::string("fitpoints") );
 
  registerProcessorParameter ("DutZsColName", "DUT zero surpressed data colection name",
 		     _dutZsColName, std::string("zsdata_apix"));
 
  registerProcessorParameter ("OutputPath", "Path/File where root-file should be stored",
			      _path2file, std::string("NTuple.root"));
  
  registerProcessorParameter ("DUTIDs", "Int std::vector containing the IDs of the DUTs",
		  		_DUTIDs, std::vector<int>());
  
  registerOptionalParameter("ReferenceCollection","This is the name of the reference it collection (init at 0,0,0)",
                           _referenceHitCollectionName, static_cast< string > ( "reference" ) );
  
  registerOptionalParameter("AlignmentCollectionNames", "Names of alignment collections, should be in same order as application", _alignColNames, std::vector<std::string>());
			   
			
 
}


void EUTelAPIXTbTrackTuple::init()
{
	// usually a good idea to
	printParameters();

	_isFirstEvent = true; 
 
	_nRun = 0;
	_nEvt = 0;
	
	 _foundAllign = false;

  	_rotationstored.clear();
  	_countrotstored = 0;

	// Prepare TTree/TFiles
	prepareTree();
	
#ifndef USE_GEAR
	streamlog_out ( ERROR4 ) <<  "Marlin was not built with GEAR support." << endl <<  "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;
	exit(-1);
#else
	// check if the GEAR manager pointer is not null!
	if ( Global::GEAR == 0x0 ) 
	{
		streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << std::endl;
		exit(-1);
	}

	_siPlanesParameters  = const_cast<SiPlanesParameters* > ( &(Global::GEAR->getSiPlanesParameters()) );
	_siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout()) );

	invertGear();
	
#endif
}

void EUTelAPIXTbTrackTuple::processRunHeader( LCRunHeader* runHeader) 
{
	std::auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
	eutelHeader->addProcessor( type() );
	_nRun++;
	
	// Decode and print out Run Header information - just a check
	_runNr = runHeader->getRunNumber();
}

void EUTelAPIXTbTrackTuple::processEvent( LCEvent * event )
{
	_nEvt ++;
	_evtNr = event->getEventNumber();
	EUTelEventImpl* euEvent = static_cast<EUTelEventImpl*> ( event );

	if( euEvent->getEventType() == kEORE )
       	{
		streamlog_out( DEBUG5) << "EORE found: nothing else to do." << std::endl;
    		return;
	}
  	if( not _foundAllign ) { readAlignment(event); } //Sets _foundAlign to true if found.
  	if( not _foundAllign ) {
    		streamlog_out  ( ERROR5 ) << "Have not found the needed alignment collections, will skip this event ( " 
			     << event->getEventNumber() << " )." << endl; 
    		return;
  		}
	//Clear all event info containers
	clear();

	//try to read in hits (e.g. fitted hits in local frame)	
  	if( !readHits( _inputTrackerHitColName, event )) 
	{ 
		return; 
	}
  	
	//read in raw data	
	if(!readZsHits( _dutZsColName , event )) 
	{ 
		return; 
	}
  
  	//read in tracks
	if(!readTracks(event))
       	{
		return;
	}
 
        //fill the trees	
	_zstree->Fill();
	_eutracks->Fill();
	_euhits->Fill();

	_isFirstEvent = false;
}

void EUTelAPIXTbTrackTuple::end()
{
	//Maybe some stats output?
	_file->Write();
}

//Read in TrackerHit(Impl) to later dump them
bool EUTelAPIXTbTrackTuple::readHits( std::string hitColName, LCEvent* event )
{
	LCCollection* hitCollection = NULL;
  
	try
       	{
		hitCollection = event->getCollection( hitColName ); 
  	} 
	catch(lcio::DataNotAvailableException& e) 
	{
		streamlog_out( DEBUG2 ) << "Hit collection " << hitColName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
    		return false;
  	}
  	
	int nHit = hitCollection->getNumberOfElements();
	_nHits = nHit;
 
  	for(int ihit=0; ihit< nHit ; ihit++)
       	{
    		TrackerHitImpl* meshit = dynamic_cast<TrackerHitImpl*>( hitCollection->getElementAt(ihit) ) ;
    		const double* pos = meshit->getPosition();	
    		LCObjectVec clusterVec = (meshit->getRawHits());

		UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
    		int sensorID = hitDecoder(meshit)["sensorID"];

		//Only dump DUT hits
		if( std::find( _DUTIDs.begin(), _DUTIDs.end(), sensorID) == _DUTIDs.end() )
		{
			continue;
		}

    		double x = pos[0];
    		double y = pos[1];
    		double z = pos[2];
    
    		_hitXPos->push_back(x);
    		_hitYPos->push_back(y);
    		_hitZPos->push_back(z);
    		_hitSensorId->push_back(sensorID);
	}

	return true;
}

//Read in TrackerHit to later dump
bool EUTelAPIXTbTrackTuple::readTracks(LCEvent* event)
{
	LCCollection* trackCol = NULL;

	try
	{
    		trackCol = event->getCollection( _inputTrackColName ) ;
  	}
	catch(lcio::DataNotAvailableException& e)
	{
		streamlog_out( DEBUG2 ) << "Track collection " << _inputTrackColName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
		return false;
  	}

	// setup cellIdDecoder to decode the hit properties
	UTIL::CellIDDecoder<TrackerHitImpl>  hitCellDecoder( EUTELESCOPE::HITENCODING );

	int nTrackParams=0;
	for(int itrack=0; itrack< trackCol->getNumberOfElements(); itrack++)
       	{
		lcio::Track* fittrack = dynamic_cast<lcio::Track*>( trackCol->getElementAt(itrack) ) ;
		
		std::vector<EVENT::TrackerHit*> trackhits = fittrack->getTrackerHits();
		double chi2 = fittrack->getChi2();
		double ndof = fittrack->getNdf();
		double dxdz = fittrack->getOmega();
		double dydz = fittrack->getPhi();

 		//Get the (fitted) hits belonging to this track, they are in local frame!
		for(unsigned int ihit=0; ihit< trackhits.size(); ihit++)
	       	{
      			TrackerHitImpl* fittedHit = dynamic_cast<TrackerHitImpl*>( trackhits.at(ihit) );
      			const double* pos = fittedHit->getPosition();
      			if( (hitCellDecoder(fittedHit)["properties"] & kFittedHit) == 0)
		       	{
			       	continue;
		       	}

     			UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
      			int sensorID = hitDecoder(fittedHit)["sensorID"];
	  		double nominalZpos = _siPlanesLayerLayout->getSensitivePositionZ(ihit);

			//Dump the (fitted) hits for the DUTs
			if( std::find( _DUTIDs.begin(), _DUTIDs.end(), sensorID) == _DUTIDs.end() )
			{
				continue;
			}

      			nTrackParams++;
      			
			double x = pos[0];
      			double y = pos[1];
      			double z = pos[2]; //not used!
			reverseAlign(x,y,z,sensorID, nominalZpos);
			//eutrack tree
      			_xPos->push_back(x);
      			_yPos->push_back(y);
      			_dxdz->push_back(dxdz);
      			_dydz->push_back(dydz);
      			_trackIden->push_back(sensorID);
      			_trackNum->push_back(itrack);
      			_chi2->push_back(chi2);
      			_ndof->push_back(ndof);
    		}
  	}

	_nTrackParams = nTrackParams;
	return true;
}

//Read in raw (zs) TrackerData(Impl) to later dump
bool EUTelAPIXTbTrackTuple::readZsHits( std::string colName, LCEvent* event)
{
	LCCollectionVec* zsInputCollectionVec = NULL;
 
	try
	{
		zsInputCollectionVec = dynamic_cast<LCCollectionVec*>( event->getCollection(colName) );
  	}
       	catch(DataNotAvailableException& e)
	{
		streamlog_out( DEBUG2 ) << "Raw ZS data collection " << colName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
    		return false;
  	}
	
	UTIL::CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
	for( unsigned int plane = 0; plane < zsInputCollectionVec->size(); plane++ ) 
	{
		TrackerDataImpl* zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( plane ) );
		SparsePixelType type = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
    		int sensorID = cellDecoder( zsData )["sensorID"];
    
		if (type == kEUTelGenericSparsePixel  ) 
		{
			std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > apixData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ));
      			EUTelGenericSparsePixel apixPixel;
      				
			for( unsigned int iHit = 0; iHit < apixData->size(); iHit++ ) 
			{
				apixData->getSparsePixelAt( iHit, &apixPixel);
				_nPixHits++;
				p_iden->push_back( sensorID );
				p_row->push_back( apixPixel.getYCoord() );
				p_col->push_back( apixPixel.getXCoord() );
				p_tot->push_back( static_cast< int >(apixPixel.getSignal()) );
				p_lv1->push_back( static_cast< int >(apixPixel.getTime()) );
      			}
    		}
		else
		{
			//PANIC
		}
  	}
  	return true;
}

void EUTelAPIXTbTrackTuple::clear()
{
  /* Clear zsdata */
  p_col->clear();
  p_row->clear();
  p_tot->clear();
  p_iden->clear();
  p_lv1->clear();
  _nPixHits = 0;
  /* Clear hittrack */
  _xPos->clear();
  _yPos->clear();
  _dxdz->clear();
  _dydz->clear();
  _trackNum->clear();
  _trackIden->clear();
  _chi2->clear();
  _ndof->clear();
  //Clear hits
  _hitXPos->clear();
  _hitYPos->clear();
  _hitZPos->clear();
  _hitSensorId->clear();
 }

void EUTelAPIXTbTrackTuple::prepareTree()
{
	_file = new TFile(_path2file.c_str(),"RECREATE");

	_xPos = new std::vector<double>();	     
	_yPos = new std::vector<double>();	   
	_dxdz = new std::vector<double>();	   
	_dydz = new std::vector<double>();	   
	_trackIden  = new std::vector<int>();
	_trackNum = new std::vector<int>();
	_chi2 = new std::vector<double>();	   
	_ndof = new std::vector<double>();    

	p_col = new std::vector<int>();
	p_row = new std::vector<int>();
	p_tot = new std::vector<int>();
	p_iden = new std::vector<int>();
	p_lv1 = new std::vector<int>();

	_hitXPos = new std::vector<double>();
	_hitYPos = new std::vector<double>();
	_hitZPos = new std::vector<double>();
	_hitSensorId  = new std::vector<int>();

	_rotDUTId = new vector<int>();
	_rotXY = new vector<double>();
	_rotZX = new vector<double>();
	_rotZY = new vector<double>();
	
	_euhits = new TTree("fitpoints","fitpoints");
	_euhits->Branch("nHits", &_nHits);
	_euhits->Branch("xPos", &_hitXPos);
	_euhits->Branch("yPos", &_hitYPos);
	_euhits->Branch("zPos", &_hitZPos);
	_euhits->Branch("sensorId", &_hitSensorId);

	_zstree = new TTree("rawdata", "rawdata");
	_zstree->Branch("nPixHits", &_nPixHits);
	_zstree->Branch("euEvt",    &_nEvt);
	_zstree->Branch("col",      &p_col);
	_zstree->Branch("row",      &p_row);
	_zstree->Branch("tot",      &p_tot);
	_zstree->Branch("lv1",      &p_lv1);
	_zstree->Branch("iden",     &p_iden);

	//Tree for storing all track param info
	_eutracks = new TTree("tracks", "tracks");
	_eutracks->Branch("nTrackParams", &_nTrackParams);
	_eutracks->Branch("euEvt", &_nEvt);
	_eutracks->Branch("xPos", &_xPos);
	_eutracks->Branch("yPos", &_yPos);
	_eutracks->Branch("dxdz", &_dxdz);
	_eutracks->Branch("dydz", &_dydz);
	_eutracks->Branch("trackNum", &_trackNum);
	_eutracks->Branch("iden", &_trackIden);
	_eutracks->Branch("chi2", &_chi2);
	_eutracks->Branch("ndof", &_ndof);
	
	
	//TTree for DUT rot info
	_rottree = new TTree("DUTrotation", "DUTrotation");
	_rottree->Branch("ID", &_rotDUTId);
	_rottree->Branch("Alpha", &_alpha);
	_rottree->Branch("Beta", &_beta);
	_rottree->Branch("Gamma", &_gamma);
	_rottree->Branch("RotXY", &_rotXY);
	_rottree->Branch("RotZX", &_rotZX);
	_rottree->Branch("RotZY", &_rotZY);
	_rottree->Branch("RotXYErr", &_rotXYerr);
	_rottree->Branch("RotZXErr", &_rotZXerr);
	_rottree->Branch("RotZYErr", &_rotZYerr);

	_euhits->AddFriend(_zstree);
	_euhits->AddFriend(_eutracks);

	_versionVec = new TVectorD(1);
	_versionVec[0] = 1.1;
        _versionVec->Write("ver");
}



gsl_matrix* EUTelAPIXTbTrackTuple::invertLU(int dim, gsl_matrix* matrix){
  gsl_permutation* perm = gsl_permutation_alloc(dim);
  int s = 0;
  gsl_matrix * inverse = gsl_matrix_alloc(dim,dim);
  gsl_linalg_LU_decomp(matrix, perm, &s);
  gsl_linalg_LU_invert(matrix, perm, inverse);
  gsl_permutation_free(perm);
  //gsl_matrix_fprintf(stdout,inverse,"%f");
  return (inverse);
}



void EUTelAPIXTbTrackTuple::invertAlignment(EUTelAlignmentConstant * alignment){
  int iden = alignment->getSensorID();
  
  //Rotations
  gsl_matrix * alignM = gsl_matrix_alloc(3,3);
  //Fill alignment matrix
  if(_doScales){
    gsl_matrix_set( alignM, 0, 0, (1.0 + alignment->getAlpha() ));
    gsl_matrix_set( alignM, 1, 1, (1.0 + alignment->getBeta()  ));
    gsl_matrix_set( alignM, 2, 2, 1.0 );
    gsl_matrix_set( alignM, 0, 1, alignment->getGamma() );
    gsl_matrix_set( alignM, 1, 0, -1 * alignment->getGamma() );
    gsl_matrix_set( alignM, 0, 2, 0.0 );
    gsl_matrix_set( alignM, 2, 0, 0.0 );
    gsl_matrix_set( alignM, 1, 2, 0.0 );
    gsl_matrix_set( alignM, 2, 1, 0.0 );
  } else {
    gsl_matrix_set( alignM, 0, 0, 1.0);
    gsl_matrix_set( alignM, 1, 1, 1.0);
    gsl_matrix_set( alignM, 2, 2, 1.0 );
    gsl_matrix_set( alignM, 0, 1, -1* alignment->getGamma() ); // Rubinsky 09-10-2011
    gsl_matrix_set( alignM, 1, 0,     alignment->getGamma());
    gsl_matrix_set( alignM, 0, 2,     alignment->getBeta() );
    gsl_matrix_set( alignM, 2, 0, -1* alignment->getBeta());
    gsl_matrix_set( alignM, 1, 2, -1* alignment->getAlpha() );
    gsl_matrix_set( alignM, 2, 1,     alignment->getAlpha());
  }
  message<DEBUG5> ( log() << "Inverting alignment matrix for iden" << iden  ) ;
  gsl_matrix * inverse = invertLU(3, alignM);
  gsl_matrix_free(alignM);
  _alignRot[iden].push_back(inverse);

  //Shifts
  std::vector<double> shifts;
  shifts.push_back( alignment->getXOffset()  );
  shifts.push_back( alignment->getYOffset() );
  shifts.push_back( alignment->getZOffset() );
  _alignShift[iden].push_back(shifts);
  //Rotations
  std::vector<double> rotations;
  rotations.push_back( alignment->getAlpha()  );
  rotations.push_back( alignment->getBeta() );
  rotations.push_back( alignment->getGamma() );
  _alignRotations[iden].push_back(rotations);

  //check if alignment for plane has already been stored
  if(_rotationstored[iden] == false){
    getDUTRot(alignment);
    _rotationstored[iden] = true;
  }

  streamlog_out( DEBUG5 ) << "Iden: " << iden << endl
			 << "X-shift: "<< alignment->getXOffset() << endl
			 << "Y-shift: "<< alignment->getYOffset() << endl
			 << "Z-shift: "<< alignment->getZOffset() << endl
			 << "alpha  : "<< alignment->getAlpha() << endl
			 << "beta   : "<< alignment->getBeta() << endl
			 << "gamma  : "<< alignment->getGamma()<< endl;
}





void EUTelAPIXTbTrackTuple::invertGear(){
  for ( int layerIndex = 0 ; layerIndex < _siPlanesParameters->getSiPlanesNumber() ; ++layerIndex ) {
    int iden = _siPlanesLayerLayout->getID( layerIndex );
    //Rotations
    gsl_matrix * gearM = gsl_matrix_alloc(2,2);
    gsl_matrix_set( gearM, 0, 0, _siPlanesLayerLayout->getSensitiveRotation1(layerIndex));
    gsl_matrix_set( gearM, 0, 1, _siPlanesLayerLayout->getSensitiveRotation2(layerIndex));
    gsl_matrix_set( gearM, 1, 0, _siPlanesLayerLayout->getSensitiveRotation3(layerIndex));
    gsl_matrix_set( gearM, 1, 1, _siPlanesLayerLayout->getSensitiveRotation4(layerIndex));

    //Ugly as sin, but following the suggestions of EUTelHitMaker
    double xSign(0.0), ySign(0.0);
    if( (gsl_matrix_get(gearM, 0 ,0) < -0.7) or (gsl_matrix_get(gearM,0,1) < -0.7)){
      xSign = -1.0;
    } else if( (gsl_matrix_get(gearM, 0 ,0) > 0.7) or (gsl_matrix_get(gearM,0,1) > 0.7)){
      xSign = 1.0;
    }
    if( (gsl_matrix_get(gearM, 1 ,0) < -0.7) or (gsl_matrix_get(gearM,1,1) < -0.7)){
      ySign = -1.0;
    } else if( (gsl_matrix_get(gearM, 1 ,0) > 0.7) or (gsl_matrix_get(gearM,1,1) > 0.7)){
      ySign = 1.0;
    }
  
    message<DEBUG5> ( log() << "Inverting gear matrix for iden" << iden  ) ;
    gsl_matrix * inverse = invertLU(2, gearM);
    gsl_matrix_free(gearM);
    _gearRot[iden] = inverse;
    
    //shifts
    vector<double> shifts;
    shifts.push_back(-1 * ( _siPlanesLayerLayout->getSensitivePositionX(layerIndex) 
			    - xSign * 0.5 * _siPlanesLayerLayout->getSensitiveSizeX(layerIndex)));
    shifts.push_back(-1 * ( _siPlanesLayerLayout->getSensitivePositionY(layerIndex)
			    - ySign * 0.5 * _siPlanesLayerLayout->getSensitiveSizeY(layerIndex)));
    shifts.push_back(-1 * _siPlanesLayerLayout->getSensitivePositionZ(layerIndex));
    _gearShift[iden] = shifts;

    //offset
    vector<double> offset;
    offset.push_back( _siPlanesLayerLayout->getSensitivePositionX(layerIndex) );
    offset.push_back( _siPlanesLayerLayout->getSensitivePositionY(layerIndex) );
    offset.push_back( _siPlanesLayerLayout->getSensitivePositionZ(layerIndex) );
    _gearOffset[iden] = offset;

    //sizes
    vector<double> size;
    size.push_back(-1 * (  - xSign * 0.5 * _siPlanesLayerLayout->getSensitiveSizeX(layerIndex)));
    size.push_back(-1 * (  - ySign * 0.5 * _siPlanesLayerLayout->getSensitiveSizeY(layerIndex)));
    _gearSize[iden] = size;


    //Pitches
    _gearPitch[iden] = make_pair( _siPlanesLayerLayout->getSensitivePitchX(layerIndex),
				  _siPlanesLayerLayout->getSensitivePitchY(layerIndex));

    std::vector<double> gearRotations(3, 0.0);
    //Convert angles to radians
    double conv = 3.1415926/180.0;
    gearRotations.at(0) = conv * _siPlanesLayerLayout->getLayerRotationXY(layerIndex);
    gearRotations.at(1) = conv * _siPlanesLayerLayout->getLayerRotationZX(layerIndex);
    gearRotations.at(2) = conv * _siPlanesLayerLayout->getLayerRotationZY(layerIndex);
    streamlog_out ( DEBUG5 )  << "Plane iden: " << iden << endl;
    streamlog_out ( DEBUG5 )  << "gearRotationsXY = " << gearRotations.at(0) << endl
                             << "gearRotationsZX = " << gearRotations.at(1) << endl
                             << "gearRotationsZY = " << gearRotations.at(2) << endl;
    _gearEulerRot[iden] = gearRotations;
  }
}



void EUTelAPIXTbTrackTuple::getDUTRot(EUTelAlignmentConstant * alignment){
  int iden = alignment->getSensorID();

  for ( int layerIndex = 0 ; layerIndex < _siPlanesParameters->getSiPlanesNumber() ; ++layerIndex ) {
    int idencheck = _siPlanesLayerLayout->getID( layerIndex );

    if (idencheck == iden){
      //get rotations from gearfile
      double rotXY = _siPlanesLayerLayout->getLayerRotationXY(layerIndex);
      double rotZX = _siPlanesLayerLayout->getLayerRotationZX(layerIndex);
      double rotZY = _siPlanesLayerLayout->getLayerRotationZY(layerIndex);

      _rotDUTId->push_back(iden);

      //corrections from alignment
      _alpha->push_back(alignment->getAlpha());
      _beta->push_back(alignment->getBeta());
      _gamma->push_back(alignment->getGamma());

      //rotations from gearfile
      _rotXY->push_back(rotXY);
      _rotZX->push_back(rotZX);
      _rotZY->push_back(rotZY);

      //errors from alignment
      _rotXYerr->push_back(alignment->getAlphaError());
      _rotZXerr->push_back(alignment->getBetaError());
      _rotZYerr->push_back(alignment->getGammaError());
    }
  }

  //count stored planes
  _countrotstored++;

  //message<MESSAGE5> ( log() << "Planes stored " << _countrotstored );
  //message<MESSAGE5> ( log() << "# Planes " << _siPlanesParameters->getSiPlanesNumber() );

  //fill tree omly once
  //first plane has no alignment
  if ((_countrotstored + 1) == _siPlanesParameters->getSiPlanesNumber()){
    _rottree->Fill();
  }
}



void EUTelAPIXTbTrackTuple::readAlignment(LCEvent* event){
  //Check for alignment collections, if found invert and store, if not keep _foundAlign false so it will look in next event.
  _foundAllign = true;

  LCCollectionVec * alignmentCollectionVec;
  streamlog_out ( DEBUG5 )  << "Trying to read " << _alignColNames.size() << " collections" << endl;
  for(size_t ii = 0; ii < _alignColNames.size(); ii++){
    size_t index = _alignColNames.size() - ii -1;
    streamlog_out ( DEBUG5 )  << "Trying to read alignment collection " << _alignColNames.at(index) << endl;
    try{
      alignmentCollectionVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _alignColNames.at(index)));
      for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) {
	 streamlog_out ( DEBUG5 ) << "Inverting plane " << iPos << endl;
	invertAlignment( static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) ) );
      }
    } catch( DataNotAvailableException& e) {
      streamlog_out ( ERROR5 ) << "Could not find  alignment collections on event " << event->getEventNumber() << "." << endl;
      _foundAllign = false;
    }
  }
}



void EUTelAPIXTbTrackTuple::reverseAlign(double& x, double& y, double &z, int iden, double nomZpos){
//printf("reverseAlign-- %5d -- x:%5.2f %5.2f %5.2f \n", iden, x,y,z);
  // Apply alignment translations 
  double xTemp(0.0),yTemp(0.0), zTemp(0.0);

  if( _referenceHitVec != 0  )
  {
    int iloc = -1;
    for(int ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
    {
      EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
      if( iden == refhit->getSensorID()) iloc = ii;
    }
   
    if( iloc < 0 ) streamlog_out( WARNING ) << "Uknown sensor found at :" << iden << endl;
    EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(iloc) ) ;
    if( refhit == 0 ) streamlog_out( WARNING ) << "refhit is empty :" << refhit << endl;
    double  xrefhit = refhit->getXOffset();
    double  yrefhit = refhit->getYOffset();
    double  zrefhit = refhit->getZOffset();

    double temp_x = x;
    double temp_y = y;
    double temp_z = z;

    if( abs(xrefhit) < 1e-03 && abs(xrefhit) < 1e-03 )
    {
       for( size_t trans = 0; trans < _alignShift[iden].size(); trans++)
       { 
         xrefhit -= _alignShift[iden].at(trans).at(0);
         yrefhit -= _alignShift[iden].at(trans).at(1);
         zrefhit -= _alignShift[iden].at(trans).at(2);
       }
    }   
//    printf("APIXTbTrackTuple: original refhit[%5.3f %5.3f %7.3f] local[%5.3f %5.3f %7.3f] unaligned[%5.3f %5.3f %7.3f]\n", 
//                                     xrefhit, yrefhit, zrefhit, 
//                                     xTemp, yTemp, zTemp, 
//                                     temp_x, temp_y, temp_z);


// alignment(s) part
    for( size_t trans = 0; trans < _alignShift[iden].size(); trans++)
    { 
      xrefhit += _alignShift[iden].at(trans).at(0);
      yrefhit += _alignShift[iden].at(trans).at(1);
      zrefhit += _alignShift[iden].at(trans).at(2);

//      arefhit = _RotatedVector[0];
//      brefhit = _RotatedVector[1];
//      grefhit = _RotatedVector[2];
      
      //if( _alignShift.find(iden) != _alignShift.end() ){
      double x_rot = temp_x - xrefhit;
      double y_rot = temp_y - yrefhit;
      double z_rot = temp_z - zrefhit;

      //Apply alignment rotations
      TVector3 iCenterOfSensorFrame(x_rot, y_rot, z_rot);
      iCenterOfSensorFrame.RotateZ( _alignRotations[iden].at(trans).at(2) );
      iCenterOfSensorFrame.RotateY( _alignRotations[iden].at(trans).at(1) );
      iCenterOfSensorFrame.RotateX( _alignRotations[iden].at(trans).at(0) );
      xTemp = iCenterOfSensorFrame[0];
      yTemp = iCenterOfSensorFrame[1];
      zTemp = iCenterOfSensorFrame[2];
 
//      gsl_matrix* m = _alignRot[iden].at(trans);
//      xTemp = x_rot * gsl_matrix_get(m,0,0) + y_rot * gsl_matrix_get(m,0,1) + z_rot * gsl_matrix_get(m,0,2); 
//      yTemp = x_rot * gsl_matrix_get(m,1,0) + y_rot * gsl_matrix_get(m,1,1) + z_rot * gsl_matrix_get(m,1,2); 
//      zTemp = x_rot * gsl_matrix_get(m,2,0) + y_rot * gsl_matrix_get(m,2,1) + z_rot * gsl_matrix_get(m,2,2); 
      temp_x = xTemp + xrefhit + _alignShift[iden].at(trans).at(0);
      temp_y = yTemp + yrefhit + _alignShift[iden].at(trans).at(1);
      temp_z = zTemp + zrefhit + _alignShift[iden].at(trans).at(2);

//      printf("APIXTbTrackTuple: align[%1d] refhit[%5.3f %5.3f %7.3f] local[%5.3f %5.3f %7.3f] unaligned[%5.3f %5.3f %7.3f]\n", trans, 
//                                     xrefhit, yrefhit, zrefhit, 
//                                     xTemp, yTemp, zTemp, 
//                                     temp_x, temp_y, temp_z);

      //Do not need z's from here, I guess
    }

    x = temp_x;
    y = temp_y;
    z = temp_z;


// gear part
    if(_gearOffset.find(iden) != _gearOffset.end())
    {
      x -= _gearOffset[iden].at(0);
      y -= _gearOffset[iden].at(1);
      z -= _gearOffset[iden].at(2);
    }
 
    TVector3 RotatedSensorHit( x, y, z);

    std::vector<double>& rots = _gearEulerRot[iden];
    if( TMath::Abs(rots.at(0)) > 1e-6 ){
      RotatedSensorHit.RotateZ( -1.0 * rots.at(0) ); // in XY
    }
    if( TMath::Abs(rots.at(1)) > 1e-6 ){
      RotatedSensorHit.RotateY( -1.0 * rots.at(1) ); // in ZX 
    }
    if( TMath::Abs(rots.at(2)) > 1e-6 ){
      RotatedSensorHit.RotateX( -1.0 * rots.at(2) ); // in ZY
    }
    x = RotatedSensorHit.X();
    y = RotatedSensorHit.Y();
    z = RotatedSensorHit.Z();
//printf("unrotated x:%5.2f y:%5.2f z:%5.2f\n", x,y,z);
 
    if(_gearSize.find(iden) != _gearSize.end())
    {
      x += _gearSize[iden].at(0);
      y += _gearSize[iden].at(1);
    }
//printf("size corrected x:%5.2f y:%5.2f \n", x,y);
 
    //Apply gear rot
    if(_gearRot.find(iden) != _gearRot.end()){
      gsl_matrix* m = _gearRot[iden];
      xTemp = x * gsl_matrix_get(m,0,0) + y * gsl_matrix_get(m,0,1);
      yTemp = x * gsl_matrix_get(m,1,0) + y * gsl_matrix_get(m,1,1);
      x = xTemp; y = yTemp;
    }

 
  }else{
 
   for( size_t trans = 0; trans < _alignShift[iden].size(); trans++){
    //if( _alignShift.find(iden) != _alignShift.end() ){
    x += _alignShift[iden].at(trans).at(0);
    y += _alignShift[iden].at(trans).at(1);
    z += _alignShift[iden].at(trans).at(2);
    double zShift(z);
    if(not _doScales){ zShift = z - nomZpos;}
    //Apply alignment rotations
    gsl_matrix* m = _alignRot[iden].at(trans);
    xTemp = x * gsl_matrix_get(m,0,0) + y * gsl_matrix_get(m,0,1) + zShift * gsl_matrix_get(m,0,2); 
    yTemp = x * gsl_matrix_get(m,1,0) + y * gsl_matrix_get(m,1,1) + zShift * gsl_matrix_get(m,1,2); 
    zTemp = x * gsl_matrix_get(m,2,0) + y * gsl_matrix_get(m,2,1) + z * gsl_matrix_get(m,2,2); 
    x = xTemp; y = yTemp; z = zTemp;
    //Do not need z's from here, I guess
   }

   zTemp = z - nomZpos;
    TVector3 RotatedSensorHit( x, y, zTemp);
    std::vector<double>& rots = _gearEulerRot[iden];
    if( TMath::Abs(rots.at(0)) > 1e-6 ){
      RotatedSensorHit.RotateZ( -1.0 * rots.at(0) ); // in XY
    }
    if( TMath::Abs(rots.at(1)) > 1e-6 ){
      RotatedSensorHit.RotateY( -1.0 * rots.at(1) ); // in ZX 
    }
    if( TMath::Abs(rots.at(2)) > 1e-6 ){
      RotatedSensorHit.RotateX( -1.0 * rots.at(2) ); // in ZY
    }
    x = RotatedSensorHit.X();
    y = RotatedSensorHit.Y();
 
    //Apply gear local trans
    if(_gearShift.find(iden) != _gearShift.end())
    {
      x += _gearShift[iden].at(0);
      y += _gearShift[iden].at(1);
    }
    //Apply gear rot
    if(_gearRot.find(iden) != _gearRot.end()){
      gsl_matrix* m = _gearRot[iden];
      xTemp = x * gsl_matrix_get(m,0,0) + y * gsl_matrix_get(m,0,1);
      yTemp = x * gsl_matrix_get(m,1,0) + y * gsl_matrix_get(m,1,1);
      x = xTemp; y = yTemp;
    }

    // Change reference from sensor corner to center of pixel 0,0
    x -= 0.5 * _gearPitch[iden].first;
    y -= 0.5 * _gearPitch[iden].second;
  }

//  printf("iden %5d at x=%5.2f y=%5.2f z=%5.2f \n", iden, x, y, z);
}











