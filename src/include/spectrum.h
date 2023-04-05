/*
 * spectrum.h
 *
 *  Created on: 17 Oct 2013
 *      Author: drw
 */

#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include <fstream>
using namespace std;

template <typename T>
class Spectrum
{
public:
	T *energy, *counts;
	int bins;

private:
	int get_size(const char* filename)
	{
		ifstream spec_file;
		T en, c;

		int count = 0;

		spec_file.open(filename);

		while( !spec_file.eof() )
		{
		  spec_file >> en >> c;
		  count++;
		}

		count--;

		spec_file.close();

		return count;
	}

	int get_qdp_size(const char* filename)
	{
		ifstream spec_file;

		T en, err, c;

		int count = 0;

		spec_file.open(filename);

		while( !spec_file.eof() )
		{
		  spec_file >> en >> err >> c;

		  // throw it away if this line is not in the correct format
		  if(spec_file.good())
		  {
			  count++;
		  }
		  else
		  {
			  spec_file.clear();
		  	  spec_file.ignore(10,'\n');
		  }
		}

		count--;

		spec_file.close();

		return count;
	}

	void read_spectrum(const char* filename)
	{
		//
		// read in emission lines from a plain text file
		//
		// arguments:
		//   filename    char*     name of file to read emission lines from
		//   n_bins      int&      stores number of lines
		//   line_en     gpu_float*   stores the line energy
		//   line_strength gpu_float* stores the relative photon counts in each line
		//
		ifstream lines_file;

		int count = 0;

		lines_file.open(filename);

		while( !lines_file.eof() )
		{
			if(count >= bins) break;

			lines_file >> energy[count] >> counts[count];
			count++;
		}

		lines_file.close();
	}

	void read_qdp_spectrum(const char* filename)
	{
		//
		// read in emission lines from a QDP text file
		//
		// arguments:
		//   filename    char*     name of file to read emission lines from
		//   n_bins      int&      stores number of lines
		//   line_en     gpu_float*   stores the line energy
		//   line_strength gpu_float* stores the relative photon counts in each line
		//
		ifstream spec_file;

		T en, err, cts;

		int count = 0;

		spec_file.open(filename);

		while( !spec_file.eof() )
		{
		  spec_file >> en >> err >> cts;

		  // throw it away if this line is not in the correct format
		  if(spec_file.good())
		  {
			  if(count >= bins) break;

			  energy[count] = en;
			  counts[count] = cts;
			  count++;
		  }
		  else
		  {
			  spec_file.clear();
		  	  spec_file.ignore(10,'\n');
		  }
		}

		spec_file.close();
	}
public:
	Spectrum(char* filename, bool qdp = false)
	{
		cout << "Reading spectrum: " << filename << " (" << ( (qdp) ? "QDP" : "plain text" ) << ')' << endl;

		bins = (qdp) ? get_qdp_size(filename) : get_size(filename);

		energy = new T[bins];
		counts = new T[bins];

		if(qdp)
			read_qdp_spectrum(filename);
		else
			read_spectrum(filename);

		cout << "Spectrum has " << bins << " bins " << energy[0] << ':' << energy[bins-1] << endl;
	}

	~Spectrum( )
	{
		delete[] energy;
		delete[] counts;
	}
};


#endif /* SPECTRUM_H_ */
