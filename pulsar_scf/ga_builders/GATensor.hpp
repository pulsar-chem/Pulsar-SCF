#pragma once
#include<array>
#include<vector>
#include<utility>
#include<pulsar/math/TensorImpl.hpp>
#include<pulsar/math/EigenImpl.hpp>
#include<pulsar/math/Cast.hpp>

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
 *  GA's C API is not very explicit about who owns memory and when.  I've
 *  modeled this class somewhat after the std::vector.  The short and sweet
 *  explanation of this statement is this class owns the memory associated with
 *  it.
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
 *  change its shape you'll have to make a new instance.  This means the
 *
 */
class GATensor{
public:
    ///The type of a dimension's limits
    using limit_t = std::pair<size_t,size_t>;

    ///The type of a tensor's shape
    using shape_t = std::vector<limit_t>;

    ///The type of the buffer
    using buffer_t = std::shared_ptr<std::vector<double>>;

    ///Makes a GATensor from an existing handle
    ///(call set_value to give this class the memory)
    GATensor(int handle):handle_(handle){}

    ///Makes a GATensor with dimensions given by \p dims and blocked via \p chunks
    template<size_t N>
    GATensor(const std::array<size_t,N>& dims,const std::array<size_t,N>& chunks,
             const char * name=""):
        GATensor(N,dims.data(),chunks.size()?chunks.data():nullptr,name)
    {}

    ///Makes a GATensor with dimensions \p dims, blocks determined automatically
    template<size_t N>
    explicit GATensor(const std::array<size_t,N>& dims,const char * name=""):
        GATensor(dims,std::array<size_t,N>(),name){}

    ///Makes a GATensor with dimension \p dims, set to value
    template<size_t N>
    GATensor(const std::array<size_t,N>& dims,double value,const char* name=""):
        GATensor(dims,name)
    {
        fill(value);
    }

    ///Sets a single value
    template<size_t N>
    void set_value(const std::array<size_t,N>& idx, double value)const{
        set_value(idx.data(),idx.data(),&value);
    }

    ///Sets a block from a buffer, uses normal C++ [lo,hi) range
    void set_value(const shape_t& shape,const double* buffer);

    ///Sets all elements of the tensor to a value
    void fill(double value)const;

    ///Returns a single value
    template<size_t N>
    double get_value(const std::array<size_t,N>& idx)const{
        return (*get_value(idx.data(),idx.data()))[0];
    }

    ///Returns a block, not necessarily the local one
    buffer_t get_value(const shape_t& shape)const;

    ///Returns the local block
    buffer_t my_data()const;

    ///Returns the tensor's (local) shape
    shape_t my_shape()const{return shape_;}

    ///Function to check if an index is local
    template<size_t N>
    bool my_idx(const std::array<size_t,N>& idx)const{
        for(size_t i=0;i<N;++i)
            if(idx[i]<shape_[i].first || idx[i]>=shape_[i].second)
                return false;
        return true;
    }

    ///Returns the tensor's full shape
    const std::vector<size_t>& dims()const{return dims_;}

    ///Prints the tensor for debuggin'
    void print_out() const;

    ///Returns the tensor's handle
    int handle()const{return handle_;}

private:
    ///The handle assigned to this tensor
    int handle_;

    ///The rank of the tensor (it's not a non-type template parameter to
    ///hide GA from the outside world)
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
    buffer_t get_value(const size_t* low,
                                  const size_t* high)const;

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

    BPHASH_DECLARE_HASHING_FRIENDS

    void hash(bphash::Hasher & h) const
    {
        //This is a very bad hash...
        h(tensor_->handle());
    }
public:
    using tensor_type=GATensor;
    using shared_tensor=std::shared_ptr<tensor_type>;
    /*! \brief For serialization only
     *
     * \warning NOT FOR USE OUTSIDE OF SERIALIZATION
     * \todo Replace if cereal fixes this
     */
    GATensorImpl() = default;

    ///Aliases this GA tensor to \p other
    GATensorImpl(const shared_tensor& other):
        tensor_(other)
    {}

    ///Aliases underlying GA tensor
    GATensorImpl(const tensor_type& other)
        : tensor_(std::make_shared<tensor_type>(other)) { }

    ///Moves GA tensor from \p mat
    GATensorImpl(tensor_type && other)
        : tensor_(std::make_shared<tensor_type>(std::move(other))) { }

    ///True if the underlying Eigen matrices are the same
    bool operator==(const GATensorImpl& rhs)const
    {
        return *tensor_==*rhs.tensor_;
    }

    ///True if the underlying Eigen matrices are not the same
    bool operator!=(const GATensorImpl& rhs)const
    {
        return !((*this)==rhs);
    }

    /*! \brief Obtain a hash of the data
         *
         * Details depend on what kind of data is stored.
         * See the hashing functions of the stored type
         * for details.
    */
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

//Allow us to cast from Pulsar's stub pointer back
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
            rv.set_value(shape,eigen_mat.data());
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
