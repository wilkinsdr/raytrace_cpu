//
// Created by Dan Wilkins on 7/6/17.
//

#ifndef RAYTRACE_PAR_ARGS_H_H
#define RAYTRACE_PAR_ARGS_H_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <exception>
#include <typeinfo>
#include <vector>
using namespace std;

#define SEP_CHAR '='
#define KEY_PREFIX "--"

class ArgumentException : public exception
{
	string what_msg;
public:
	ArgumentException(string msg) : what_msg("Argument ERROR : " + msg)
	{

	}
	virtual ~ArgumentException() throw()
	{   }

	virtual const char* what() const throw()
	{
		return what_msg.c_str();
	}
};


class ParameterArgs
{
private:
    map <string, string> pars;
    vector <string> positional_args;

    void trim_whitespace(string &str)
    {
        str.erase(0, str.find_first_not_of("\t "));
        str.erase(str.find_last_not_of("\t ") + 1);
    }

    bool is_key(string &line)
    {
        string temp = line;

        trim_whitespace(temp);
        size_t prefix_pos = temp.find(KEY_PREFIX);
        if(prefix_pos == temp.npos)
            return false;

        return true;
    }

    bool is_key_pair(string &line)
    {
        string temp = line;

        trim_whitespace(temp);
        size_t sep_pos = temp.find(SEP_CHAR);
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
			throw ArgumentException( string("Duplicate definition of ") + key);
	}

	void parse_argv(int argc, char** argv)
    {
        for(int i = 1; i < argc; i++)
        {
            string line(argv[i]);
            if(is_key(line))
                parse_key_value(line);
            else
                positional_args.push_back(line);
        }
    }

public:
	ParameterArgs(int argc, char** argv)
	{
		parse_argv(argc, argv);
	}

	bool key_exists(string const& key) const
	{
		return pars.find(key) != pars.end();
	}

	template <typename T>
	T get_parameter(string const& key) const
	{
		if(!key_exists(key))
			throw ArgumentException( string("Required parameter ") + key + " not specified");

		istringstream istr(pars.find(key)->second);
		T parsed_val;
		if (!(istr >> parsed_val))
		{
			throw ArgumentException( string("Could not parse value of ") + key + " (expected type " + typeid(T).name() + ")" );
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
			throw ArgumentException( string("Could not parse value of ") + key + " (expected type " + typeid(T).name() + ")" );
		}
		return parsed_val;
	}

	template <typename T>
	void get_parameter_array(string const& key, T* arr, int n_elem) const
	{
		if(!key_exists(key))
			throw ArgumentException( key + " not found in parameter file");

		istringstream istr(pars.find(key)->second);

//		int read_values;
		for(int i=0; i<n_elem; i++)
		{
			if (!(istr >> arr[i]))
			{
				ostringstream errormsg;
				errormsg << "Could not parse array " << key << " (expected " << n_elem << " values of type " << typeid(T).name() << ")";
				throw ArgumentException( errormsg.str() );
			}
		}
	}

	string get_string_parameter(string const& key) const
	{
		if(!key_exists(key))
			throw ArgumentException( key + " not found in parameter file");

		string value = pars.find(key)->second;
		return value;
	}

	template <typename T>
	static T string_to_T(std::string const &str)
    {
        istringstream istr(str);
        T parsed_val;
        if(!(istr >> parsed_val))
        {
            cerr << "ParameterFile ERROR : Could not parse value";
            return 0;
        }
        return parsed_val;
    }

    int num_positional()
    {
        return positional_args.size();
    }

    string operator[](int i)
    {
        if(i < 0 || i >= positional_args.size())
        {
            ostringstream errormsg;
            errormsg << "Positional argument " << i + 1 << " not supplied";
            throw ArgumentException(errormsg.str());
        }
        return positional_args[i];
    }
};

#endif //RAYTRACE_PAR_ARGS_H_H
