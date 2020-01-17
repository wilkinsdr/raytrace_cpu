//
// Created by drw on 15/01/2020.
//

#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <string.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
using namespace std;

class ProgressBar
{
private:
	long end;
	int length;
	char* label;
	int digits;

public:
	ProgressBar(long end, char* init_label = "", int init_length = 0) : end(end), length(init_length), label(init_label)
	{
		if(length == 0)
		{
			struct winsize w;
			ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
			int label_length = strlen(label);
			digits = 0;
			int temp_end = end;
			while(temp_end)
			{
				temp_end /= 10;
				++digits;
			}
			length = w.ws_col - label_length - (2*digits) - 5;
		}
		cout << "\e[?25l";
	}

	void show(long prog)
	{
		ostringstream outstream;

		outstream << "\e[?25l";
		outstream << "\r" << label << ' ' << setw(digits) << prog << '/' << end << " [";
		for(int step=1; step<=length; step++)
		{
			outstream << ((((float)step/length) <= ((float)prog/end)) ? "\u2588" : "-");
		}
		outstream << ']';
		outstream << "\e[?25h";

		cout << outstream.str();
	}

	void show_complete()
	{
		ostringstream outstream;

		outstream << "\e[?25l";
		outstream << "\r" << label << ' ' << setw(digits) << end << '/' << end << " [";
		for(int step=1; step<=length; step++)
		{
			outstream << "\u2588";
		}
		outstream << ']';
		outstream << "\e[?25h";

		cout << outstream.str();
	}

	void done()
	{
		cout << endl;
	}
};

#endif //PROGRESS_BAR_H
