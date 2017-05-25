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

}//End namespace detail_

//This CTor uses the fact that members are initialized in declaration order
GATensor::GATensor(size_t rank,
                   const size_t* dims,
                   const size_t* chunks,
                   const char * name):
    handle_(NGA_Create(
        C_DBL,
        rank,
        detail_::to_GAint(rank,dims).data(),
        const_cast<char*>(name),
        chunks? detail_::to_GAint(rank,chunks).data() : nullptr
    )),
    rank_(rank),
    dims_(dims,dims+rank),
    shape_([&](){
        int proc=GA_Nodeid();
        std::vector<int> lo(rank_),hi(rank_);
        return detail_::get_shape(handle_,rank_,proc,lo.data(),hi.data());
    }()),
    nelems_(detail_::nelems(shape_))
{
}

void GATensor::set_value(const size_t* lo,
                         const size_t* hi,
                         const double* value)const{
    std::vector<GAint> _hi(detail_::to_GAint(rank_,hi)),
                       _lo(detail_::to_GAint(rank_,lo)),
                      strides(detail_::unit_stride(rank_));
    NGA_Put(
        handle_,
        _lo.data(),
        _hi.data(),
        const_cast<double*>(value),
        strides.data()
    );
}


GATensor::buffer_t GATensor::get_value(const size_t* low,const size_t* high)const{
    size_t total=1;
    for(size_t i=0;i<rank_;++i)
        total*=high[i]-low[i]+1;
    auto buffer=std::make_shared<std::vector<double>>(total);
    std::vector<GAint> strides(detail_::unit_stride(rank_)),
                       _high(detail_::to_GAint(rank_,high)),
                       _low(detail_::to_GAint(rank_,low));
    NGA_Get(
        handle_,
        _low.data(),
        _high.data(),
        buffer->data(),
        strides.data()
    );
    return buffer;
}

void GATensor::set_value(const shape_t& shape,const double* buffer)
{
    std::vector<size_t> lo,hi;
    std::tie(lo,hi) = detail_::shape_2_GA(shape);
    set_value(lo.data(),hi.data(),const_cast<double *>(buffer));
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

GATensor symmetrize(const GATensor& tensor)
{
    GA_Symmetrize(tensor.handle());
    return tensor;
}

}//End namespace pulsar_scf
