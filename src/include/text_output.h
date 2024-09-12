/*
 * text_output.h
 *
 * Class for outputting data into columnar text files
 * Use the stream operator << to output a column of data to the file
 *
 *  Created on: 22 Mar 2016
 *      Author: drw
 */

#ifndef TEXT_OUTPUT_H_
#define TEXT_OUTPUT_H_

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

class TextOutput
{
private:
	ofstream outfile;
	bool open;

	int width;

public:
	TextOutput(const char* filename, bool append = false, int precision = 8, int colwidth = 20, ios_base::fmtflags format = ios::scientific)
		: width(colwidth)
	{
		cout << "Writing results to text file: " << filename << endl;
		outfile.open(filename, (append) ? ios::app : ios::out);

		if(outfile)
		{
			open = true;
		}
		else
		{
			cerr << "************" << endl;
			cerr << "TextOutput ERROR: Could not open output file " << filename << endl;
			cerr << "************" << endl;
			open = false;
			return;
		}

		outfile.setf(format);
		outfile.precision(precision);
	}

	TextOutput(string filename, bool append = false, int precision = 8, int colwidth = 20, ios_base::fmtflags format = ios::scientific)
			: width(colwidth)
	{
		cout << "Writing results to text file: " << filename << endl;
		outfile.open(filename.c_str(), (append) ? ios::app : ios::out);

		if(outfile)
		{
			open = true;
		}
		else
		{
			cerr << "************" << endl;
			cerr << "TextOutput ERROR: Could not open output file " << filename << endl;
			cerr << "************" << endl;
			open = false;
			return;
		}

		outfile.setf(format);
		outfile.precision(precision);
	}

	void close( )
	{
		outfile.close();
		open = false;
	}

	~TextOutput( )
	{
        close();
	}

	void set_format(ios_base::fmtflags format = ios::scientific)
	{
		outfile.setf(format);
	}

	void set_precision(int precision = 8)
	{
		outfile.precision(precision);
	}

	void newline(int n = 1)
	{
		for(int i=0; i<n; i++)
			outfile << endl;
	}

	// insertion operator for manipulators
	TextOutput& operator << ( ostream& (*m)(ostream&) )
	{
	    if(open)
	    	outfile << *m;
	    else
	    	cerr << "TextOutput ERROR: File is not open" << endl;

	    return *this;
	}

	// insertion operators for variables
	template<typename T>
	friend TextOutput& operator << (TextOutput& t, T const& data);

};

template<typename T>
TextOutput& operator << (TextOutput& t, T const& data)
{
	if(t.open)
		t.outfile << setw(t.width) << data;
	else
		cerr << "TextOutput ERROR: File is not open" << endl;
	return t;
}


#endif /* TEXT_OUTPUT_H_ */
