/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELTIMEPIX3SPARSEPIXEL_H
#define EUTELTIMEPIX3SPARSEPIXEL_H

// personal includes ".h"
#include "EUTelBaseSparsePixel.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>

namespace eutelescope {

// Structure to store pixel info
struct PIXEL
	{
		unsigned char  x, y;
		unsigned short tot;
		uint64_t ts;
	};

// Structure to store trigger info
struct TRIGGER
	{
		unsigned short int_nr, tlu_nr;
		uint64_t ts;
	};


  //! Helper class for simple sparsified pixel
  /*! This class contains only the pixel coordinates and signal as
   *  integer numbers.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */

  class EUTelTimepix3SparsePixel  :  public EUTelBaseSparsePixel  {

  public:

    //! Default constructor with all arguments
	EUTelTimepix3SparsePixel::EUTelTimepix3SparsePixel( unsigned char xCoord,
														unsigned char yCoord,
														unsigned short signal,
														uint64_t hit_timestamp,
														unsigned short int_number,
														unsigned short TLU_number,
														uint64_t tlu_timestamp);
    //! Default constructor with no args
    EUTelTimepix3SparsePixel();

    //! Default destructor
    virtual ~EUTelTimepix3SparsePixel() { ; }

        //! Get the number of elements in the data structure
    /*! This method returns the number of elements the sparse pixel
     *  contains.
     *
     *  @return The number of elements in the data structure
     */
    virtual unsigned int getNoOfElements() const;

  //! Get the sparse pixel type using the enumerator
    /*! This methods returns the sparse pixel type using the
     *  enumerator defined in EUTELESCOPE.h
     *
     *  @return The sparse pixel type using the enumerator
     */
    virtual SparsePixelType getSparsePixelType() const;

    //! Print method
    /*! This method is used to print out the contents of the sparse
     *  pixel
     *
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const ;

    //! Setter for x coordinate
    void setXCoord(unsigned char xCoord) { _pixel_data.x = xCoord ; }

    //! Setter for y coordinate
    void setYCoord(unsigned char yCoord) { _pixel_data.y = yCoord ; }

    //! Setter for the signal
    void setSignal(unsigned short signal) { _pixel_data.tot = signal ; }

    //! Getter for the x coordinate
    inline unsigned short getXCoord() const { return (unsigned short)_pixel_data.x ; }

    //! Getter for the y coordinate
    inline unsigned short getYCoord() const { return (unsigned short)_pixel_data.y ; }

    //! Getter for the signal
    inline unsigned short getSignal() const { return  _pixel_data.tot ; }


  private:
    PIXEL _pixel_data;
    TRIGGER _trigger_data;
  };
}

#endif
