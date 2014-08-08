// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id$

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelTimepix3SparsePixel.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace eutelescope;


EUTelTimepix3SparsePixel::EUTelTimepix3SparsePixel()
{
 _pixel_data.x=0;
 _pixel_data.y=0;
 _pixel_data.tot=0;
 _pixel_data.ts=0;

 _trigger_data.int_nr=0;
 _trigger_data.tlu_nr=0;
 _trigger_data.ts=0;

  _noOfElements = 7;
  _type = kEUTelTimepix3SparsePixel;
}

EUTelTimepix3SparsePixel::EUTelTimepix3SparsePixel( unsigned char xCoord,
													unsigned char yCoord,
													unsigned short signal,
													uint64_t hit_timestamp,
													unsigned short int_number,
													unsigned short TLU_number,
													uint64_t tlu_timestamp)
{

	_pixel_data.x=xCoord;
	_pixel_data.y=yCoord;
	_pixel_data.tot=signal;
	_pixel_data.ts=hit_timestamp;

	_trigger_data.int_nr=int_number;
	_trigger_data.tlu_nr=TLU_number;
	_trigger_data.ts=tlu_timestamp;

	_noOfElements = 7;
  _type = kEUTelTimepix3SparsePixel;
}


unsigned int EUTelTimepix3SparsePixel::getNoOfElements() const {
  return _noOfElements;
}

SparsePixelType EUTelTimepix3SparsePixel::getSparsePixelType() const {
  return _type;
}

void EUTelTimepix3SparsePixel::print(std::ostream& os) const {
  int bigWidth = 50;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
  os << endl;
  int width = 20;
  os << setw(width) << setiosflags(ios::left) << "Type: "     << _type << endl
     << setw(width) << "Elements: " << _noOfElements << endl
     << setw(width) << "x coord: "  << _pixel_data.x << endl
     << setw(width) << "y coord: "  << _pixel_data.y << endl
     << setw(width) << "signal: "  << _pixel_data.tot<< endl;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
}
