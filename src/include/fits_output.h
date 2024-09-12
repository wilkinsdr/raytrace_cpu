/*
 * fits_output.h
 *
 *  Created on: 17 Oct 2013
 *      Author: drw
 */

#ifndef FITS_OUTPUT_H_
#define FITS_OUTPUT_H_

#include "fitsio.h"
#include <string>

#define COL_STRING "8A"
#define COL_STRING16 "16A"
#define COL_STRING32 "32A"
#define COL_STRING64 "64A"
#define COL_BOOL = "1L"
#define COL_INT "1I"
#define COL_UINT "1U"
#define COL_INT32 "1J"
#define COL_UINT32 "1V"
#define COL_INT64 "1K"
#define COL_FLOAT "1E"
#define COL_FLOAT64 "1D"
#define COL_COMPLEX "1C"
#define COL_COMPLEX64 "1M"

class FITSOutputException : public exception
{
	string what_msg;
public:
	FITSOutputException(string msg, int status = 0) : what_msg("FITSOutput ERROR : " + msg)
	{
		if(status > 0) fits_report_error(stderr, status);
	}
	virtual ~FITSOutputException() throw()
	{   }

	virtual const char* what() const throw()
	{
		return what_msg.c_str();
	}
};

template <typename T>
class FITSOutput
{
private:
	fitsfile* fptr;
	int status;
	bool open;

	int nextcol, tabrows;

	void print_error(int status)
	{
	    //
	    // Print out cfitsio error messages
	    //
	    if (status)
	    {
	    	fits_report_error(stderr, status);
	    }
	    return;
	}

public:
	FITSOutput(char* filename, bool clobber = true) : status(0)
	{
		cout << "Opening FITS file for output: " << filename;
		if(clobber) cout << " (will overwrite if file exists)" << endl;

		if(clobber) remove(filename);

		if( fits_create_file(&fptr, filename, &status) )
		{
			open = false;
			throw FITSOutputException( string("Could not open file"), status);
		}
		open = true;
	}

	FITSOutput(string filename, bool clobber = true) : status(0)
	{
		cout << "Opening FITS file for output: " << filename;
		if(clobber) cout << " (will overwrite if file exists)" << endl;

		if(clobber) remove(filename.c_str());

		if( fits_create_file(&fptr, filename.c_str(), &status) )
		{
			open = false;
			throw FITSOutputException( string("Could not open file"), status);
		}
		open = true;
	}

	void close( )
	{
		if(open) fits_close_file(fptr, &status);

		open = false;
	}

	~FITSOutput( )
	{
		close( );
	}

	void create_primary( )
	{
		//
		// Create an empty image extension to be used as the primary header in a multi-extension file
		//
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}

		long naxes[2];
		naxes[0]=0;
		naxes[1]=0;

		cout << "Creating empty primary extension in FITS file" << endl;

		if( fits_create_img(fptr, BYTE_IMG, 0, naxes, &status) ) print_error(status);
	}


	// ------ Image output functions ------

	void write_image_array(double* frame, int Nx, int Ny)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}

		long naxes[2];
		naxes[0]=static_cast<long>(Nx);
		naxes[1]=static_cast<long>(Ny);

		long startpixel[2];
		startpixel[0] = 1;
		startpixel[1] = 1;

		cout << "Adding " << Nx << 'x' << Ny << " image extension to FITS file" << endl;

		if( fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status) ) print_error(status);
		if( fits_write_pix(fptr, TDOUBLE, startpixel, Nx*Ny, frame, &status) ) print_error(status);

	}

	void write_image(T** data, int Nx, int Ny, bool transpose = false, bool flip_x = false, bool flip_y = false)
	{
		//
		// writes a 2 dimentional array data[x][y] to a FITS image extension
		//
		// X runs from left to right, Y from bottom to top (unless flipped)
		// In transpose mode, X runs from bottom to top, Y from left to right (unless flipped)
		// Note that flip_x and flip_y correspond to the original image array, not the transposed output
		//
		double* frame;

		frame = new double[Nx*Ny];

		// flatten the 2D image into a 1D array for cfitsio
		// note the reversed ordering of x/y axes between input and output
		for(int j=0; j<Nx; j++)
			for(int k=0; k<Ny; k++)
			{
				const int x = (flip_x) ? (Nx - 1 - j) : j;
				const int y = (flip_y) ? (Ny - 1 - k) : k;

				if(transpose)
					frame[j*Ny + k] = data[x][y];
				else
					frame[k*Nx + j]	= data[x][y];

			}

		if (transpose) write_image_array(frame, Ny, Nx);
		else write_image_array(frame, Nx, Ny);

		delete[] frame;
	}

	void write_me_data_cube(T*** data, int Nframes, int Nx, int Ny)
	{
		//
		// Write a multi-extension datacube (series of image frames) to a FITS file
		//
		double* frame;

		cout << "Writing Frame..." << endl;
		for(int i=0; i<Nframes; i++)
		{
			frame = new double[Ny*Nx];

			for(int j=0; j<Ny; j++)
				for(int k=0; k<Nx; k++)
					frame[(Ny - 1 - j)*Nx + (Nx - 1 - k)] = data[i][j][k];	// (flipping the y-axis)

			write_image_array(frame, Nx, Ny);

			delete[] frame;
		}
	}


	// ------ Table output functions ------

	void create_table(const char* extname, int cols, long rows, char *col_names[], char *col_format[], char *col_units[])
	{
		//
		// create a new binary table extension
		//
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}

		cout << "Adding binary table with " << cols << " columns and " << rows << " rows to FITS file" << endl;

		if ( fits_create_tbl(fptr, BINARY_TBL, rows, cols, col_names, col_format, col_units, extname, &status) )
			throw FITSOutputException( string("Could not create table extension"), status );

		// save these for convenience
		nextcol = 1;
		tabrows = rows;
	}

	void write_table_column(double* data, long rows = -1, int col = -1, long firstrow = 1, long firstelem = 1)
	{
		//
		// write a column to the last created table extension
		//
		cout << "Writing column with " << rows << " rows to FITS table extension" << endl;

		if(rows < 0) rows = tabrows;
		if(col < 0) col = nextcol;

		if( fits_write_col(fptr, TDOUBLE, col, firstrow, firstelem, rows, data, &status) )
			throw FITSOutputException( string("Could not add column to table extension"), status );

		nextcol = col + 1;
	}

	void write_table_column(int* data, long rows = -1, int col = -1, long firstrow = 1, long firstelem = 1)
	{
		//
		// write a column to the last created table extension
		//
		cout << "Writing column with " << rows << " rows to FITS table extension" << endl;

		if(rows < 0) rows = tabrows;
		if(col < 0) col = nextcol;

		if( fits_write_col(fptr, TINT, col, firstrow, firstelem, rows, data, &status) )
			throw FITSOutputException( string("Could not add column to table extension"), status);

		nextcol = col + 1;
	}


	// ------ Keywords ------

	//
	// The following (overloaded) function writes a keyword or a comment to the header of the previously created extension
	// N.B. Must create an extension first. For the primary header, use close_primary() to make an empty extension
	//
	void write_keyword(const char* keyname, const char* comment, int value)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_key(fptr, TINT, keyname, &value, comment, &status) ) print_error(status);
	}
	void write_keyword(const char* keyname, const char* comment, long value)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_key(fptr, TLONG, keyname, &value, comment, &status) ) print_error(status);
	}
	void write_keyword(const char* keyname, const char* comment, double value)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_key(fptr, TDOUBLE, keyname, &value, comment, &status) ) print_error(status);
	}
	void write_keyword(const char* keyname, const char* comment, const char* value)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_key(fptr, TSTRING, keyname, (char*)value, comment, &status) ) print_error(status);
	}
	void write_keyword(const char* keyname, const char* comment, char* value)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_key(fptr, TSTRING, keyname, value, comment, &status) ) print_error(status);
	}
	void write_keyword(const char* keyname, const char* comment, string value)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_key(fptr, TSTRING, keyname, (char*)value.c_str(), comment, &status) ) print_error(status);
	}
	void write_keyword(const char* keyname, const char* comment, bool value)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		int logical_value = (value) ? 1 : 0;
		if( fits_write_key(fptr, TLOGICAL, keyname, &logical_value, comment, &status) ) print_error(status);
	}

	void write_comment(const char* comment)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_comment(fptr, (char*)comment, &status) ) print_error(status);
	}

	void write_comment(char* comment)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		if( fits_write_comment(fptr, comment, &status) ) print_error(status);
	}

	void set_ext_name(char* extname)
	{
		if(!open)
		{
			throw FITSOutputException( string("File is not open") );
		}
		char keyname[] = "EXTNAME";
		char comment[] = "Name of this extension";
		if( fits_write_key(fptr, TSTRING, keyname, extname, comment, &status) ) print_error(status);
	}
};


#endif /* FITS_OUTPUT_H_ */
