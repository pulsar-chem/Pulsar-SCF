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
    nelems_(0),
    transpose_(false)
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
    nelems_=detail_::nelems(shape_);
}

GATensor::GATensor(GATensor&& other):
    handle_(std::move(other.handle_)),
    rank_(std::move(other.rank_)),
    dims_(std::move(other.dims_)),
    shape_(std::move(other.shape_)),
    nelems_(std::move(other.nelems_)),
    transpose_(std::move(other.transpose_))
{
    //We own the memory now
    other.handle_=0;
}

GATensor::GATensor(const GATensor& other):
    handle_(GA_Create_handle()),
    rank_(other.rank_),
    dims_(other.dims_),
    shape_(other.shape_),
    nelems_(other.nelems_),
    transpose_(other.transpose_)
{
    auto temp=detail_::to_GAint(rank_,dims_.data());
    GA_Set_data(handle_,rank_,temp.data(),C_DBL);
    GA_Allocate(handle_);
    GA_Copy(other.handle_,handle_);
}

GATensor::GATensor():
    handle_(0),
    rank_(0),
    dims_(),
    shape_(),
    nelems_(0),
    transpose_(false)
{

}

GATensor::~GATensor(){
    //The other half of RAII, don't let the object go away if you still want it
    if(handle_)
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

GATensor vec2matrix(const GATensor &tensor)
{
    const size_t len=tensor.dims()[0];
    GATensor rv(std::array<size_t,2>({len,len}),0.0);
    auto data=tensor.my_data();
    auto rvdata=rv.my_data();
    //TODO: case when all data isn't local
    for(size_t i=0;i<len;++i)
        rvdata[i*len+i]=data[i];
    rv.set_value(rv.my_shape(),rvdata.data());
    return rv;

}

//TODO: Figure out how to get this PEIGS library and then actually call GA's eigensolver
std::pair<GATensor,GATensor> GeneralEigenSolver(const GATensor& tensor,
                                                const GATensor& metric)
{
    if(tensor.dims().size()!=2)
        throw pulsar::PulsarException("Not sure how to diagonalize a vector or proper tensor");
    const size_t len=tensor.dims()[0];
    GATensor values(std::array<size_t,1>({len})),
             vectors(std::array<size_t,2>({tensor.dims()[0],tensor.dims()[1]}));
    //std::vector<double> scratch(len,0.0);
    //Assume it's all local
    auto data=tensor.my_data(),metric_data=metric.my_data();
    Eigen::Map<Eigen::MatrixXd> H(data.data(),len,len);
    Eigen::Map<Eigen::MatrixXd> S(metric_data.data(),len,len);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver(H,S);
    auto P = eig_solver.eigenvectors().transpose();
    //Eigenvalues are in increasing order
    auto E = eig_solver.eigenvalues();
    values.set_value(values.my_shape(),E.data());
    //GA_Diag(tensor.handle(),metric.handle(),vectors.handle(),scratch.data());
    vectors.set_value(vectors.my_shape(),P.data());
    return std::make_pair(values,vectors);
}

std::pair<GATensor,GATensor> EigenSolver(const GATensor& tensor)
{
    if(tensor.dims().size()!=2)
        throw pulsar::PulsarException("Not sure how to diagonalize a vector or proper tensor");
    const size_t len=tensor.dims()[0];
    GATensor values(std::array<size_t,1>({len})),
              vectors(std::array<size_t,2>({len,len}));
    //std::vector<double> scratch(len,0.0);
    //GA_Diag_std(tensor.handle(),vectors.handle(),scratch.data());
    auto data=tensor.my_data();
    Eigen::Map<Eigen::MatrixXd> H(data.data(),len,len);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver(H);
    Eigen::MatrixXd P = eig_solver.eigenvectors().transpose();
    //Eigenvalues are in increasing order
    auto E = eig_solver.eigenvalues();
    values.set_value(values.my_shape(),E.data());
    vectors.set_value(vectors.my_shape(),P.data());
    //GA_Diag(tensor.handle(),metric.handle(),vectors.handle(),scratch.data());
    return std::make_pair(values,vectors);
}

GATensor GApow(const GATensor& tensor,double power)
{
    GATensor rv(tensor);
    auto data=rv.my_data();
    for(double& x :data)x=std::pow(x,power);
    rv.set_value(rv.my_shape(),data.data());
    return rv;
}

GATensor gemm(double alpha,const GATensor& A,const GATensor&B)
{
    std::array<size_t,2> dims({A.dims()[A.is_transposed()],
                               B.dims()[!B.is_transposed()]});
    GATensor C(dims,0.0);
    return gemm(alpha,A,B,1.0,C);
}

GATensor gemm(double alpha,const GATensor& A,const GATensor&B, double beta,const GATensor& C)
{
    if(A.rank()==0||A.rank()>2||
       B.rank()==0||B.rank()>2||
       C.rank()==0||C.rank()>2)
        throw pulsar::PulsarException("NYI gemm with rank!= 1 or 2");
    //TODO: Skip copy if A and B are not modified by GA
    GATensor rv(C),_A(A),_B(B);
    char ta=A.is_transposed()?'T':'N';
    char tb=B.is_transposed()?'T':'N';

    //A and B can be vectors (so only grab 1st dimension)
    //Vectors are column vectors
    GAint m=A.dims()[0];//Right if A is not transposed
    GAint n=B.dims()[0];//Right if B is transposed
    GAint k=A.dims()[0];//Right if A is transposed
    if(ta=='T')m=(A.rank()==1? 1 : A.dims()[1]);
    if(tb=='N')n=(B.rank()==1? 1 : B.dims()[1]);
    if(ta=='N')k=(A.rank()==1? 1 : A.dims()[1]);
    GA_Dgemm(ta,tb,m,n,k,alpha,_A.handle(),_B.handle(),beta,rv.handle());
    return rv;
}

void GAaccumulate(GATensor& B,double alpha,const GATensor& A)
{
    double one=1.0;
    GA_Add(&one,B.handle(),&alpha,A.handle(),B.handle());
}

GATensor GAAdd(double alpha,const GATensor& A, double beta, const GATensor& B)
{
    GATensor C(A.dims());
    GA_Add(&alpha,A.handle(),&beta,B.handle(),C.handle());
    return C;
}


}//End namespace pulsar_scf
