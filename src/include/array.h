/*
 * array.h
 *
 * Template classes to manage allocations of multi-dimensional arrays
 * in contiguous memory
 *
 *  Created on: 28 May 2018
 *      Author: drw
 */

#ifndef ARRAY_H_
#define ARRAY_H_

#include <fstream>

#include "H5Cpp.h"
using namespace H5;
#include <typeinfo>

template <typename T>
class Array2D
{
public:
	T** ptr;
	//T* pool;
	int num_x, num_y;
public:
	Array2D(int Nx, int Ny) : num_x(Nx), num_y(Ny)
	{
		ptr = new T*[Nx];
		T* pool = new T[Nx*Ny];
		for(int i=0; i<Nx; ++i, pool+=Ny)
			ptr[i] = pool;

		zero();
	}

	~Array2D()
	{
		delete[] ptr[0];
		delete[] ptr;
	}

	T* operator [] (int i)
	{
		return ptr[i];
	}

	operator T**()
	{
		return ptr;
	}

	operator T*()
	{
		return ptr[0];
	}

	void zero()
	{
		for(int i=0; i<num_x*num_y; i++)
			ptr[0][i] = 0;
	}

	void write(ofstream* outfile)
	{
		outfile->write(reinterpret_cast<char*> (ptr[0]), num_x * num_y * sizeof(T));
	}

	void read(ifstream* infile)
	{
		infile->read(reinterpret_cast<char*> (ptr[0]), num_x * num_y * sizeof(T));
	}

	void write_hdf(Group* container, const char* name)
	{
		PredType type = PredType::NATIVE_DOUBLE;
		if(typeid(T) == typeid(double)) type = PredType::NATIVE_DOUBLE;
		else if(typeid(T) == typeid(float)) type = PredType::NATIVE_FLOAT;
		else if(typeid(T) == typeid(long)) type = PredType::NATIVE_LONG;
		else if(typeid(T) == typeid(int)) type = PredType::NATIVE_INT;

		hsize_t arr_dims[] = {num_x, num_y};
		DataSpace arr_dataspace(2, arr_dims);
		DataSet arr_dataset = container->createDataSet(name, type, arr_dataspace);
		arr_dataset.write(ptr[0], type);
	}

	template<typename T2>
	void operator /= (Array2D<T2>& other)
	{
		if(num_x != other.num_x || num_y != other.num_y)
		{
			cerr << "Array2D ERROR: Cannot divide arrays with different dimensions";
			return;
		}
		for(int i=0; i<num_x*num_y; i++)
			ptr[0][i] /= other.ptr[0][i];
	}

    void operator /= (float other)
    {
        for(int i=0; i<num_x*num_y; i++)
            ptr[0][i] /= other;
    }
    void operator /= (double other)
    {
        for(int i=0; i<num_x*num_y; i++)
            ptr[0][i] /= other;
    }
    void operator /= (int other)
    {
        for(int i=0; i<num_x*num_y; i++)
            ptr[0][i] /= other;
    }
    void operator /= (long other)
    {
        for(int i=0; i<num_x*num_y; i++)
            ptr[0][i] /= other;
    }

	template<typename T2>
	void operator += (Array2D<T2>& other)
	{
		if(num_x != other.num_x || num_y != other.num_y)
		{
			cerr << "Array2D ERROR: Cannot add arrays with different dimensions";
			return;
		}
		for(int i=0; i<num_x*num_y; i++)
			ptr[0][i] += other.ptr[0][i];
	}
};

template <typename T>
class Array3D
{
public:
	T*** ptr;
	T* pool;
	int num_x, num_y, num_z;
public:
	Array3D(int Nx, int Ny, int Nz) : num_x(Nx), num_y(Ny), num_z(Nz)
	{
		ptr = new T**[Nx];
		pool = new T[Nx*Ny*Nz];
		for(int i=0; i<Nx; ++i)
		{
			ptr[i] = new T*[Ny];
			for(int j=0; j<Ny; j++)
			{
				ptr[i][j] = pool + i*(Ny*Nz) + j*(Nz);
			}
		}
	}

	~Array3D()
	{
		delete[] ptr[0];
		delete[] ptr;
	}

	T** operator [] (int i)
	{
		return ptr[i];
	}

	operator T**()
	{
		return ptr[0];
	}

	operator T*()
	{
		return ptr[0][0];
	}

	void zero()
	{
		for(int i=0; i<num_x*num_y*num_z; i++)
			pool[i] = 0;
	}

	void write(ofstream* outfile)
	{
		outfile->write(reinterpret_cast<char*> (ptr[0][0]), num_x * num_y * num_z * sizeof(T));
	}

	void read(ifstream* infile)
	{
		infile->read(reinterpret_cast<char*> (ptr[0][0]), num_x * num_y * num_z * sizeof(T));
	}

	template<typename T2>
	void operator /= (Array3D<T2>& other)
	{
		if(num_x != other.num_x || num_y != other.num_y || num_z != other.num_z)
		{
			cerr << "Array2D ERROR: Cannot divide arrays with different dimensions";
			return;
		}
		for(int i=0; i<num_x*num_y*num_z; i++)
			pool[i] /= other.pool[i];
	}

    void operator /= (float other)
    {
        for(int i=0; i<num_x*num_y*num_z; i++)
            pool[i] /= other;
    }
    void operator /= (double other)
    {
        for(int i=0; i<num_x*num_y*num_z; i++)
            pool[i] /= other;
    }
    void operator /= (int other)
    {
        for(int i=0; i<num_x*num_y*num_z; i++)
            pool[i] /= other;
    }
    void operator /= (long other)
    {
        for(int i=0; i<num_x*num_y*num_z; i++)
            pool[i] /= other;
    }

	void write_hdf(Group* container, const char* name)
	{
		PredType type = PredType::NATIVE_DOUBLE;
		if(typeid(T) == typeid(double)) type = PredType::NATIVE_DOUBLE;
		else if(typeid(T) == typeid(float)) type = PredType::NATIVE_FLOAT;
		else if(typeid(T) == typeid(long)) type = PredType::NATIVE_LONG;
		else if(typeid(T) == typeid(int)) type = PredType::NATIVE_INT;

		hsize_t arr_dims[] = {num_x, num_y, num_z};
		DataSpace arr_dataspace(3, arr_dims);
		DataSet arr_dataset = container->createDataSet(name, type, arr_dataspace);
		arr_dataset.write(ptr[0][0], type);
	}
};

#endif /* ARRAY_H_ */
