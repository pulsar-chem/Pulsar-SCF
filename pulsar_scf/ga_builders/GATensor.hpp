#pragma once
#include<array>
#include<vector>
#include<utility>
#include<pulsar/math/TensorImpl.hpp>
#include<pulsar/math/EigenImpl.hpp>
#include<pulsar/math/Cast.hpp>
#include<bphash/types/All.hpp>

namespace pulsar_scf {

/** \brief A C++ wrapper to the Global Arrays library
 *
 *  Global Arrays, has C++ bindings, but they are just very thin wrappers around
 *  the C bindings.  This class is my take on a more proper Tensor class powered
 *  by Global Arrays.  In addition to serving as a more C++-like API, this
 *  class also buries the Global Array dependency in the pulsar_scf library to
 *  avoid propagating it.
 *
 *  A couple of things to note:
 *
 *  You can get the handle to the tensor so you can call GA yourself.  This is
 *  strongly discouraged.
 *
 *  As far as I can tell blocks in GA are always contiguous as far as indices go
 *  that is to say that the elements in the flattened buffer are something like
 *  (0,1), (0,2), ... and not (0,1), (0,11), (0,21)... that is the fast index has
 *  period 1 and not say 10.  Regardless of whether or not GA supports block
 *  cyclic or not, I assume the blocks have period 1.
 *
 *  In C++ it is canonical to specify ranges [begin,end), that is ranges are
 *  inclusive on the starting element and exclusive on the ending element.  When
 *  this class asks for a range I expect a C++-like range.  Conversion to GA's
 *  range standard, i.e. [begin,end], happens under the hood.
 *
 *  At the moment the shape of an instance can't be reshaped.  If you want to
 *  change its shape you'll have to make a new instance.
 *
 *  This class uses RAII, what this means is when you make an instance it
 *  obtains the memory it will need.  When the instance goes out of scope it
 *  releases the memory.  In other words, you're not responsible for the memory
 *  managed by this class.
 *
 */
class GATensor{
public:
    ///The type of a dimension's limits
    using limit_t = std::pair<size_t,size_t>;

    ///The type of a tensor's shape
    using shape_t = std::vector<limit_t>;

    ///The type of the buffer
    using buffer_t = std::vector<double>;

    ///The type of an N-dimensional index
    template<size_t N>
    using index_t=std::array<size_t,N>;

    ///Makes a GATensor from an existing handle
    GATensor(int handle):handle_(handle){}

    ///Makes a GATensor with dimensions given by \p dims and blocked via \p chunks
    template<size_t N>
    GATensor(const index_t<N>& dims,const index_t<N>& chunks,
             const char * name=nullptr):
        GATensor(N,dims.data(),chunks.size()?chunks.data():nullptr,name)
    {}

    ///Makes a GATensor with dimensions \p dims, blocks determined automatically
    template<size_t N>
    explicit GATensor(const index_t<N>& dims,const char * name=nullptr):
        GATensor(dims,index_t<N>(),name){}

    ///Makes a GATensor with dimension \p dims, set to value
    template<size_t N>
    GATensor(const index_t<N>& dims,double value,const char* name=nullptr):
        GATensor(dims,name)
    {
        fill(value);
    }

    ///Deep copies a GATensor
    GATensor(const GATensor& other);

    ///Deep copies a GATensor
    const GATensor& operator=(GATensor other){
        swap(other);
        return *this;
    }

    ///Swaps the contents of this GATensor with other
    void swap(GATensor& other);

    ///Takes ownership of other
    GATensor(GATensor&& other){swap(other);}

    ///Takes ownership of other
    const GATensor& operator=(GATensor&& other)
    {
        swap(other);
        return *this;
    }

    ///Deallocates the tensor
    ~GATensor();

    ///Sets a single value
    template<size_t N>
    void set_value(const index_t<N>& idx, double value)const{
        set_value(idx.data(),idx.data(),&value);
    }

    ///Sets a block from a buffer, uses normal C++ [lo,hi) range
    void set_value(const shape_t& shape,const double* buffer);

    ///Sets all elements of the tensor to a value
    void fill(double value)const;

    ///Returns a single value
    template<size_t N>
    double get_value(const index_t<N>& idx)const{
        return get_value(idx.data(),idx.data())[0];
    }

    ///Returns a block, not necessarily the local one
    buffer_t get_value(const shape_t& shape)const;

    ///Returns the local block
    buffer_t my_data()const;

    ///Returns the tensor's (local) shape
    shape_t my_shape()const{return shape_;}

    ///Function to check if a block is local
    bool my_block(const shape_t& shape)const;

    ///Returns the tensor's full shape (the low end is all 0's and implied)
    const std::vector<size_t>& dims()const{return dims_;}

    ///Prints the tensor for debuggin'
    void print_out() const;

    ///Returns the tensor's handle
    int handle()const{return handle_;}

private:
    ///The handle assigned to this tensor
    int handle_;

    ///The rank of the tensor
    size_t rank_;

    ///The dimensions of the tensor
    std::vector<size_t> dims_;

    ///The local shape of the tensor
    shape_t shape_;

    ///The number of local elements
    size_t nelems_;

    ///The main constructor
    GATensor(size_t rank,
             const size_t* dims,
             const size_t* chunks,
             const char * name);

    ///The main set function, uses GA's inclusive range
    void set_value(const size_t* low,
                   const size_t* hi,
                   const double* value)const;

    ///The main get function, uses GA's inclusive range
    buffer_t get_value(const size_t* low,const size_t* high)const;

    DECLARE_SERIALIZATION_FRIENDS
    BPHASH_DECLARE_HASHING_FRIENDS

    ///Makes an usable tensor (for serialization only!!!)
    GATensor()=default;

    ///Serializes the tensor (sort of...)
    template<class Archive>
    void save(Archive & archive) const
    {
        archive(handle_,rank_,dims_,shape_,nelems_);
    }

    ///Unserializes the tensor (sort of...)
    template<class Archive>
    void load(Archive & archive)
    {
        archive(handle_,rank_,dims_,shape_,nelems_);
    }

    ///Hashes the tensor
    void hash(bphash::Hasher & h) const
    {
        h(handle_,rank_,dims_,shape_,nelems_);
    }

};

//@{
///Free functions for common operations

/** \relates GATensor
 *  \brief Returns a symmetrized version of \p tensor
 */
GATensor symmetrize(const GATensor& tensor);
//@}



//Register the tensor type with Pulsar
template<size_t N>
class GATensorImpl : public pulsar::TensorImpl<N,double>{
private:
    std::shared_ptr<GATensor> tensor_;///<The actual tensor

    DECLARE_SERIALIZATION_FRIENDS
    BPHASH_DECLARE_HASHING_FRIENDS

    /*! \brief For serialization only
     *
     * \warning NOT FOR USE OUTSIDE OF SERIALIZATION
     * \todo Replace if cereal fixes this
     */
    GATensorImpl() = delete;

    template<class Archive>
    void save(Archive & archive) const
    {
       archive(tensor_);
    }

    template<class Archive>
    void load(Archive & archive)
    {
        archive(tensor_);
    }

    void hash(bphash::Hasher & h) const
    {
        h(*tensor_);
    }

public:
    using tensor_type=GATensor;
    using shared_tensor=std::shared_ptr<tensor_type>;

    ///Aliases this GA tensor to \p other
    GATensorImpl(shared_tensor other):
        tensor_(other)
    {}

    ///Aliases underlying GA tensor
    GATensorImpl(const tensor_type& other)
        : tensor_(std::make_shared<tensor_type>(other)) { }

    ///Moves GA tensor from \p mat
    GATensorImpl(tensor_type && other)
        : tensor_(std::make_shared<tensor_type>(std::move(other))) { }

    ///True if the tensors point to the same instance, TODO: check value
    bool operator==(const GATensorImpl& rhs)const
    {
        return tensor_==rhs.tensor_;
    }

    ///True if the underlying tensor are not the same instance, TODO: check_value
    bool operator!=(const GATensorImpl& rhs)const
    {
        return !((*this)==rhs);
    }

    bphash::HashValue my_hash(void) const
    {
        return bphash::make_hash(bphash::HashType::Hash128, *this);
    }

    ///\copydoc TensorImpl::sizes
    virtual std::array<size_t, N> sizes(void) const
    {
        auto dims = tensor_->dims();
        std::array<size_t,N> rv;
        for(size_t i=0;i<N;++i)rv[i]=dims[i];
        return rv;
    }

    ///\copydoc TensorImpl::get_value
    virtual double get_value(std::array<size_t,N> idx) const
    {
        return tensor_->get_value(idx);
    }

    ///\copydoc TensorImpl::set_Value
    virtual void set_value(std::array<size_t,N> idx, double val)
    {
        tensor_->set_value(idx, val);
    }

    ///Allows you to get the actual matrix (in constant form)
    std::shared_ptr<const tensor_type> get_tensor(void) const
    {
        return tensor_;
    }
};

//Allow us to cast from Pulsar's stub pointer back to a GATensor
template<size_t N>
std::shared_ptr<const GATensor>
convert_to_GA(const pulsar::TensorImpl<N,double>& ten)
{
    auto test = dynamic_cast<const GATensorImpl<N>*>(&ten);
    if(test)
        return test->get_tensor();
    else{
        //Try Eigen Matrix
        auto test2 = dynamic_cast<const pulsar::EigenMatrixImpl*>(&ten);
        if(test2)
        {
            const auto& eigen_mat = *test2->get_matrix();
            std::array<size_t,2> dims({
                pulsar::numeric_cast<size_t>(eigen_mat.rows()),
                pulsar::numeric_cast<size_t>(eigen_mat.cols())});
            GATensor::shape_t shape({{0,dims[0]},{0,dims[1]}});
            GATensor rv(dims);
            double *data=new double[dims[0]*dims[1]];
            std::copy(eigen_mat.data(),eigen_mat.data()+dims[0]*dims[1],data);
            rv.set_value(shape,data);
            return std::make_shared<GATensor>(rv);
        }
        //Try Eigen Tensor
        auto test3 = dynamic_cast<const pulsar::EigenTensorImpl<N>*>(&ten);
        if(test3)
        {
            /*const auto& eigen_tensor = *test3->get_tensor();
            std::array<size_t,N> dims({eigen_mat.rows(),eigen_mat.cols()});
            shape_t shape({{0,dims[0]},{0,dims[1]}});
            GATensor rv(dims);
            rv.set_value(shape,eigen_mat.data());
            return make_shared<GATensorImpl<2>(rv);*/
        }
    }
    throw pulsar::PulsarException("Don't know this conversion");
}

}//End namespace pulsar_scf
