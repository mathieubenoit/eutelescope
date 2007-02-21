// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: clustertest.cc,v 1.2 2007-02-21 10:55:13 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#include "lcio.h"
#include  "EUTelRunHeaderImpl.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackerRawDataImpl.h"
#include "UTIL/CellIDEncoder.h"
#include "UTIL/LCTime.h"

#include <vector>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace eutelescope;

const int nDetector  = 5;
const int nPedeEvent = 100;
const int nDataEvent = 100;

const int noiseDivider       = 3;
const int maxClusterPerEvent = 10;

const int   xCluSize  = 5;
const int   yCluSize  = 5;
const short clusterSignal[xCluSize * yCluSize] =  {0, 1, 2, 1, 0,
						   1, 5, 8, 6, 2,
						   3, 9, 13,7, 2,
						   2, 5, 8, 4, 1,
						   0, 1, 2, 1, 0};


const int xNPixel = 512;
const int yNPixel = 512;


void getXYFromIndex(int index, int& x, int& y);
int  getIndexFromXY(int x, int y);
void usage();
ofstream logfile;

int main(int argc, char ** argv) {
  
  bool doPede = false;
  bool doData = false;

  if (argc == 1) {
    usage();
    return 0;
  } else {
  
    for (int arg = 1; arg < argc; arg++) {
      if (strcmp(argv[arg],"pede") == 0) {
	doPede = true;
      }
      if (strcmp(argv[arg],"data") == 0) {
	doData = true;
      }
    }

  }

  vector<int> minX(nDetector,0);
  vector<int> minY(nDetector,0);
  vector<int> maxX(nDetector,xNPixel - 1);
  vector<int> maxY(nDetector,yNPixel - 1);

  if (doPede) {

    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
    try {
      lcWriter->open("pede_input.slcio", LCIO::WRITE_NEW);
    } catch (IOException& e) {
      cerr << e.what() << endl;
      return 0;
    }
  
    EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl(); 
    runHeader->setRunNumber(0);
    runHeader->setDetectorName("test");
    runHeader->setHeaderVersion(0.0011);
    runHeader->setDataType(EUTELESCOPE::CONVDATA);
    runHeader->setDateTime();
    runHeader->setDAQHWName(EUTELESCOPE::SUCIMAIMAGER);
    runHeader->setDAQHWVersion(0.0001);
    runHeader->setDAQSWName(EUTELESCOPE::SUCIMAIMAGER);
    runHeader->setDAQSWVersion(0.0001);  
    runHeader->setNoOfEvent(nPedeEvent);
    runHeader->setNoOfDetector(nDetector);
    runHeader->setMinX(minX);
    runHeader->setMaxX(maxX);
    runHeader->setMinY(minY);
    runHeader->setMaxY(maxY);

    lcWriter->writeRunHeader(runHeader);
    delete runHeader;


    for (int iEvent = 0; iEvent < nPedeEvent; iEvent++) {
      if ( iEvent % 10 == 0 )
	cout << "Pedestal on event " << iEvent << endl;

      LCEventImpl * event = new LCEventImpl;
      event->setDetectorName("test");
      
      LCTime * now = new LCTime;
      event->setTimeStamp(now->timeStamp());
      delete now;

      event->setEventNumber(iEvent);
      LCCollectionVec * rawData = new LCCollectionVec(LCIO::TRACKERRAWDATA);

      for (int iDetector = 0; iDetector < nDetector ; iDetector++) {

	short  baseSignal = 1;
	
	TrackerRawDataImpl * rawMatrix = new TrackerRawDataImpl;
	CellIDEncoder<TrackerRawDataImpl> idEncoder("sensorID:5,xMin:12,xMax:12,yMin:12,yMax:12", rawData);
	idEncoder["sensorID"] = iDetector;
	idEncoder["xMin"]     = 0;
	idEncoder["xMax"]     = xNPixel - 1;
	idEncoder["yMin"]     = 0;
	idEncoder["yMax"]     = yNPixel - 1;
	idEncoder.setCellID(rawMatrix);
	
	for (int yPixel = 0; yPixel < yNPixel; yPixel++) {
	  for (int xPixel = 0; xPixel < xNPixel; xPixel++) {
	    rawMatrix->adcValues().push_back(baseSignal + (short) (rand()/(RAND_MAX / noiseDivider)));
	    ++baseSignal;
	  }
	}
	rawData->push_back(rawMatrix);
      }
      event->addCollection(rawData, "rawdata");
      lcWriter->writeEvent(event);
      delete event;
    }

    lcWriter->close();

  }


  if (doData) {

    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
    try {
      lcWriter->open("data_input.slcio", LCIO::WRITE_NEW);
    } catch (IOException& e) {
      cerr << e.what() << endl;
      return 0;
    }

    EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl(); 
    runHeader->setRunNumber(0);
    runHeader->setDetectorName("test");
    runHeader->setHeaderVersion(0.0011);
    runHeader->setDataType(EUTELESCOPE::CONVDATA);
    runHeader->setDateTime();
    runHeader->setDAQHWName(EUTELESCOPE::SUCIMAIMAGER);
    runHeader->setDAQHWVersion(0.0001);
    runHeader->setDAQSWName(EUTELESCOPE::SUCIMAIMAGER);
    runHeader->setDAQSWVersion(0.0001);  
    runHeader->setNoOfEvent(nDataEvent);
    runHeader->setNoOfDetector(nDetector);
    runHeader->setMinX(minX);
    runHeader->setMaxX(maxX);
    runHeader->setMinY(minY);
    runHeader->setMaxY(maxY);

    lcWriter->writeRunHeader(runHeader);
    delete runHeader;

    short * matrix = new short[xNPixel * yNPixel];

    logfile.open("clustering.log");

    for (int iEvent = 0; iEvent < nDataEvent; iEvent++) {
      if ( iEvent % 10 == 0) 
	cout << "Data on event " << iEvent << endl;

      logfile << "Event " << iEvent << endl;
      LCEventImpl * event = new LCEventImpl;
      event->setEventNumber(iEvent);
      event->setDetectorName("test");
      
      LCTime * now = new LCTime;
      event->setTimeStamp(now->timeStamp());
      delete now;

      LCCollectionVec * rawData = new LCCollectionVec(LCIO::TRACKERRAWDATA);

      for (int iDetector = 0; iDetector < nDetector; iDetector++) {
	logfile << "Working on detector " << iDetector << endl;
	
	short  baseSignal = 1;
	
	TrackerRawDataImpl * rawMatrix = new TrackerRawDataImpl;
	CellIDEncoder<TrackerRawDataImpl> idEncoder("sensorID:5,xMin:12,xMax:12,yMin:12,yMax:12", rawData);
	idEncoder["sensorID"] = iDetector;
	idEncoder["xMin"]     = 0;
	idEncoder["xMax"]     = xNPixel - 1;
	idEncoder["yMin"]     = 0;
	idEncoder["yMax"]     = yNPixel - 1;
	idEncoder.setCellID(rawMatrix);

	int iPixel = 0;
	for (int yPixel = 0; yPixel < yNPixel; yPixel++) {
	  for (int xPixel = 0; xPixel < xNPixel; xPixel++) {
	    matrix[iPixel] = baseSignal +  + (short) (rand()/(RAND_MAX / noiseDivider));
	    ++iPixel;
	    ++baseSignal;
	  }
	}
	
	// set the number of cluster 
	int clusterPerEvent = rand() / (RAND_MAX / maxClusterPerEvent);
	logfile << "  injected " << clusterPerEvent << " clusters " << endl; 
	for (int iSeed = 0; iSeed < clusterPerEvent; iSeed++) {
	  int index = rand() / (RAND_MAX / (xNPixel * yNPixel));
	  int xSeed, ySeed;
	  getXYFromIndex(index, xSeed, ySeed);
	  logfile << "iSeed " << iSeed << " xSeed " << xSeed << " ySeed " << ySeed << endl;
	  int iCluPos = 0;
	  for (int yPixel = ySeed - (yCluSize / 2); yPixel <=  ySeed + (yCluSize / 2); yPixel++) {
	    for (int xPixel = xSeed - (xCluSize / 2); xPixel <=  xSeed + (xCluSize / 2); xPixel++) {
	      if ( ( xPixel >= 0 ) && ( xPixel < xNPixel ) &&
		   ( yPixel >= 0 ) && ( yPixel < yNPixel ) ) {
		index = getIndexFromXY(xPixel, yPixel);
		matrix[index] += clusterSignal[iCluPos];
		//		logfile << "x = " << xPixel << " y = " << yPixel << " s = " << matrix[index] << endl;
	      }
	      ++iCluPos;
	    }
	  }
	}
      
	iPixel = 0;
	for (int yPixel = 0; yPixel < yNPixel; yPixel++) {
	  for (int xPixel = 0; xPixel < xNPixel; xPixel++) {
	    rawMatrix->adcValues().push_back(matrix[iPixel]);
	    ++iPixel;
	  }
	}
	rawData->push_back(rawMatrix);
      }
      event->addCollection(rawData,"rawdata");
      lcWriter->writeEvent(event);
      delete event;
    }
    
    logfile.close();
    lcWriter->close();
    delete [] matrix;
    
  }
  return 0;
}
  
void getXYFromIndex(int index, int& x, int& y) {

  y = (index / xNPixel);
  x = index - (y * xNPixel);
}

int getIndexFromXY(int x, int y) {
  return x + y * xNPixel;
}

void usage() {
  cout << "./clustertest pede       to produce only the pedestal file " << endl
       << "./clustertest data       to produce only the data file     " << endl
       << "./clustertest pede data  to produce both! " << endl;
}
