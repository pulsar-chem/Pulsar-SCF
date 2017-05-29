#include "pulsar_scf/ga_builders/GATensor.hpp"
#include <ga.h>
#include <vector>

namespace pulsar_scf {

//Global Arrays doesn't use size_t instead they use int (modify if they add
//support for long int or even better unsigned long int)
using GAint=int;
using shape_t=GATensor::shape_t;

//Details that don't need to see the light of day
namespace detail_ {

//Makes an appropriate unit stride array for a rank "rank" tensor
inline std::vector<GAint> unit_stride(size_t rank){
    return std::vector<GAint>(rank ? rank-1 : 0,1);
}

//Converts size_t arrays to GAint arrays
inline std::vector<GAint> to_GAint(size_t len, const size_t* in){
    return std::vector<GAint>(in,in+len);
}

//Returns the number of elements in a shape
inline size_t nelems(const shape_t& shape){
    size_t total=1;
    for(size_t i=0;i<shape.size();++i)
        total*=shape[i].second-shape[i].first;
    return total;
}

//Returns the number of elements in a GA shape
inline size_t nelems(size_t rank, const size_t* lo,const size_t* hi){
    size_t total=1;
    for(size_t i=0;i<rank;++i)total*=(hi[i]-lo[i]+1);
    return total;
}

//Code factorization for getting shapes (lo/hi use GA range)
shape_t get_shape(int handle,size_t rank, int proc,GAint* lo,GAint* hi)
{
    NGA_Distribution(handle,proc,lo,hi);
    shape_t shape(rank);
    for(size_t i=0;i<rank;++i)
    {
        if(lo[i]==-1 || hi[i]==-2)//No data
            return shape_t();
        //Make the range [lo,hi) to conform to usual C++
        shape[i]=std::make_pair(lo[i],hi[i]+1);
    }
    return shape;
}

//Code factorization for converting shape to GA range
std::pair<std::vector<size_t>,std::vector<size_t>> shape_2_GA(const shape_t& shape)
{

    const size_t rank=shape.size();
    std::vector<size_t> lo(rank),hi(rank);
    for(size_t i=0;i<rank;++i)
    {
        lo[i]=shape[i].first;
        hi[i]=shape[i].second-1;
    }
    return std::make_pair(lo,hi);
}

//Wrapper around code for making initial shape
inline shape_t get_shape(size_t rank,int handle){
        int proc=GA_Nodeid();
        std::vector<int> lo(rank),hi(rank);
        return detail_::get_shape(handle,rank,proc,lo.data(),hi.data());
}

}//End namespace detail_

//This CTor uses the fact that members are initialized in declaration order
GATensor::GATensor(size_t rank,
                   const size_t* dims,
                   const size_t* chunks,
                   const char * name):
    handle_(GA_Create_handle()),
    rank_(rank),
    dims_(dims,dims+rank),
    shape_(),
    nelems_(detail_::nelems(shape_))
{
    //RAII, i.e. don't make a GATensor if you don't want the memory allocated
    auto temp=detail_::to_GAint(rank_,dims_.data());
    GA_Set_data(handle_,rank_,temp.data(),C_DBL);
    if(chunks)
    {
        temp=detail_::to_GAint(rank,chunks);
        GA_Set_chunk(handle_,temp.data());
    }
    if(name){
        std::vector<char> new_name(name,name+strlen(name));
        GA_Set_array_name(handle_,new_name.data());
    }
    GA_Allocate(handle_);
    shape_=detail_::get_shape(rank_,handle_);
}


GATensor::GATensor(const GATensor& other):
    handle_(GA_Create_handle()),
    rank_(other.rank_),
    dims_(other.dims_),
    shape_(other.shape_),
    nelems_(other.nelems_)
{
    auto temp=detail_::to_GAint(rank_,dims_.data());
    GA_Set_data(handle_,rank_,temp.data(),C_DBL);
    GA_Allocate(handle_);
    GA_Copy(other.handle_,handle_);
}

void GATensor::swap(GATensor& other)
{
    std::swap(handle_,other.handle_);
    std::swap(rank_,other.rank_);
    std::swap(dims_,other.dims_);
    std::swap(shape_,other.shape_);
    std::swap(nelems_,other.nelems_);
}

GATensor::~GATensor(){
    //The other half of RAII, don't let the object go away if you still want it
    GA_Destroy(handle_);
}

bool GATensor::my_block(const shape_t& shape)const{
    for(size_t i=0;i<shape.size();++i)
        if(shape[i].first<shape_[i].first||shape[i].second>shape_[i].second)
            return false;
    return true;
}


void GATensor::set_value(const size_t* lo,const size_t* hi,
                         const double* value)const{

    //TODO: figure out if we can just const_cast,
    //i.e. is GA going to modify value in any shape or form or take ownership?
    //TODO: avoid copy if local
    buffer_t buffer(value,value+detail_::nelems(rank_,lo,hi));
    std::vector<GAint> _hi(detail_::to_GAint(rank_,hi)),
                       _lo(detail_::to_GAint(rank_,lo)),
                      strides(detail_::unit_stride(rank_));
    NGA_Put(handle_,_lo.data(),_hi.data(),buffer.data(),strides.data());
}


GATensor::buffer_t GATensor::get_value(const size_t* low,const size_t* high)const{
    std::vector<double> buffer(detail_::nelems(rank_,low,high));
    std::vector<GAint> strides(detail_::unit_stride(rank_)),
                       _high(detail_::to_GAint(rank_,high)),
                       _low(detail_::to_GAint(rank_,low));
    //This does not allocate memory, but takes existing memory
    NGA_Get(handle_,_low.data(),_high.data(),buffer.data(),strides.data());
    return buffer;
}

void GATensor::set_value(const shape_t& shape,const double* buffer)
{
    std::vector<size_t> lo,hi;
    std::tie(lo,hi) = detail_::shape_2_GA(shape);
    set_value(lo.data(),hi.data(),buffer);
}

GATensor::buffer_t GATensor::my_data()const{
    if(!shape_.size())return buffer_t();//No data
    return get_value(shape_);
}

GATensor::buffer_t GATensor::get_value(const shape_t &shape)const
{
    std::vector<size_t> lo,hi;
    std::tie(lo,hi)=detail_::shape_2_GA(shape_);
    return get_value(lo.data(),hi.data());
}

void GATensor::print_out()const{
    GA_Print(handle_);
}

void GATensor::fill(double value)const{
    GA_Fill(handle_,&value);
}

GATensor symmetrize(const GATensor &tensor)
{
    GATensor rv(tensor);
    GA_Symmetrize(rv.handle());
    return rv;
}

}//End namespace pulsar_scf
