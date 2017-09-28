//
// Created by Dan Wilkins on 7/6/17.
//

#ifndef RAYTRACE_PAR_FILE_H_H
#define RAYTRACE_PAR_FILE_H_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <exception>
#include <typeinfo>
using namespace std;

#define COMMENT_CHAR '#'
#define SEP_CHAR '='

class ParameterException : public exception
{
	string what_msg;
public:
	ParameterException(string msg) : what_msg("ParameterFile ERROR : " + msg)
	{

	}

	virtual ~ParameterException() throw()
	{   }

	virtual const char* what() const throw()
	{
		return what_msg.c_str();
	}
};

class ParameterFile
{
private:
	map<string, string> pars;

	void remove_comment(string& line)
	{
		if (line.find(COMMENT_CHAR) != line.npos)
			line.erase(line.find(COMMENT_CHAR));
	}

	bool line_is_blank(string const& line)
	{
		return (line.find_first_not_of(' ') == line.npos);
	}

	void trim_whitespace(string& str)
	{
		str.erase(0, str.find_first_not_of("\t "));
		str.erase(str.find_last_not_of("\t ") + 1);
	}

	bool check_line(string& line)
	{
		string temp = line;

		trim_whitespace(temp);
		size_t sep_pos = line.find(SEP_CHAR);
		if(sep_pos == 0 || sep_pos == temp.npos)
			return false;

		return true;
	}

	string parse_key(string const& line)
	{
		size_t sep_pos = line.find(SEP_CHAR);
		string key = line.substr(0, sep_pos);


		// remove leading and trailing whitespace
		key.erase(0, key.find_first_not_of("\t "));
		key.erase(key.find_last_not_of("\t ") + 1);

		return key;
	}
	string parse_value(string const& line)
	{
		size_t sep_pos = line.find(SEP_CHAR);
		string value = line.substr(sep_pos + 1);

		// remove leading and trailing whitespace
		value.erase(0, value.find_first_not_of("\t "));
		value.erase(value.find_last_not_of("\t ") + 1);

		return value;
	}

	void parse_key_value(string const& line)
	{
		string key, value;
		key = parse_key(line);
		value = parse_value(line);

		if(!key_exists(key))
			pars.insert(pair<string, string>(key, value));
		else
			cerr << "ParameterFile ERROR: Duplicate definition of " << key << endl;
	}

	void parse_file(string const& filename)
	{
		ifstream par_file(filename.c_str());
		string line;

		if(!par_file.is_open())
			throw ParameterException( string("ParameterFile ERROR: Could not open file ") + filename);

		while(getline(par_file, line))
		{
			remove_comment(line);
			if(line_is_blank(line)) continue;
			parse_key_value(line);
		}

		par_file.close();
	}

public:
	ParameterFile(string filename)
	{
		parse_file(filename);
	}

	bool key_exists(string const& key) const
	{
		return pars.find(key) != pars.end();
	}

	template <typename T>
	T get_parameter(string const& key) const
	{
		if(!key_exists(key))
			throw ParameterException( string("ParameterFile ERROR: ") + key + " not found in parameter file");

		istringstream istr(pars.find(key)->second);
		T parsed_val;
		if (!(istr >> parsed_val))
		{
			throw ParameterException( string("Could not parse value of ") + key + " (expected type " + typeid(T).name() + ")" );
		}
		return parsed_val;
	}
	template <typename T>
	T get_parameter(string const& key, T default_value) const
	{
		if(!key_exists(key))
			return default_value;

		istringstream istr(pars.find(key)->second);
		T parsed_val;
		if (!(istr >> parsed_val))
		{
			throw ParameterException( string("Could not parse value of ") + key + " (expected type " + typeid(T).name() + ")" );
		}
		return parsed_val;
	}

	template <typename T>
	void get_parameter_array(string const& key, T* arr, int n_elem) const
	{
		if(!key_exists(key))
			throw ParameterException( string("ParameterFile ERROR: ") + key + " not found in parameter file");

		istringstream istr(pars.find(key)->second);

		int read_values;
		for(int i=0; i<n_elem; i++)
		{
			if (!(istr >> arr[i]))
			{
				ostringstream errormsg;
				errormsg << "Could not parse array " << key << " (expected " << n_elem << " values of type " << typeid(T).name() << ")";
				throw ParameterException( errormsg.str() );
			}
		}
	}

	char* get_string_parameter(string const& key) const
	{
		if(!key_exists(key))
			throw ParameterException( string("ParameterFile ERROR: ") + key + " not found in parameter file");

		string value = pars.find(key)->second;
		return (char*)value.c_str();
	}

	template <typename T>
	static T string_to_T(std::string const &str)
	{
		istringstream istr(str);
		T parsed_val;
		if (!(istr >> parsed_val))
		{
			cerr << "ParameterFile ERROR : Could not parse value";
			return 0;
		}
		return parsed_val;
	}
};

#endif //RAYTRACE_PAR_FILE_H_H
