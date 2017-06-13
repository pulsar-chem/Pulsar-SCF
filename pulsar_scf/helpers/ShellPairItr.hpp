#pragma once
#include<pulsar/system/BasisSet.hpp>
#include<array>

namespace pulsar_scf {

/** \brief In order to disguise the details of iterating over shell pairs one
 *         uses this iterator.
 *
 *  There are a lot of times when one wants to iterate over shell pairs, this is
 *  done with two nested for loops; however, particularly when the matrices are
 *  symmetric the logic becomes quite nasty and easy to mess up.  This class
 *  encapsulates said logic keeping your code clean.
 *
 *  This iterator is sort of an iterator within an iterator.  The outer iterator
 *  loops over the shell pairs.  Then for each shell pair it is possible to get
 *  a pair of iterators over the basis functions in that shell pair.  Typically
 *  this leads to code like:
 *
 *   \code{.cpp}
 *   ShellPairItr itr(bs);//Bs is the basis set we are iterating over
 *   while(itr)
 *   {
 *     //Do stuff with only the value of the shells
 *
 *     for(const auto x: itr)
 *     {
 *        //Do stuff with the values of the basis functions within the shell
 *     }
 *      ++itr;
 *   }
 *
 *   \endcode
 *
 *   \TODO: Possible to generalize this triplet and quartet into one class?
 *   \TODO: Allow for two basis sets and specialize for case when they are not
 *          the same
 */
class ShellPairItr{
    using array_t=std::array<size_t,2>;
    ///The starting index for each shell (last element is the number of basis functions)
    const std::vector<size_t> off_;
    ///The total number of shells
    const size_t nshells_;
    ///List of all shell pairs
    const std::vector<array_t> shell_pairs_;
    ///The index of the current pair
    size_t current_pair_;

    ///An iterator over the basis functions in the shell pair
    class BasisFunctionPairItr{
    private:
        friend class ShellPairItr;
        ///Parent ShellPairItr instance
        const ShellPairItr& parent_;
        ///Current shell i
        const size_t i_;
        ///Current shell j
        const size_t j_;
        ///Current basis function pair
        array_t basis_pair_idx_;
        const BasisFunctionPairItr& next()
        {
            if(++basis_pair_idx_[1]<parent_.off_[j_+1])
                return *this;
            ++basis_pair_idx_[0];
            basis_pair_idx_[1]=parent_.off_[j_];
            return *this;
        }

        ///Constructs a new BasisFunctionPairItr, only callable by ShellPairItr
        BasisFunctionPairItr(const ShellPairItr& parent,bool end):
            parent_(parent),
            i_(parent_.shell_pairs_[parent_.current_pair_][0]),
            j_(parent_.shell_pairs_[parent_.current_pair_][1]),
            basis_pair_idx_(array_t({parent_.off_[i_+end],parent_.off_[j_]})){}
    public:
        void reset(){
            const auto& offs=parent_.off_;
            basis_pair_idx_[0]=offs[i_];
            basis_pair_idx_[1]=offs[j_];
        }

        ///Returns the size of the block
        size_t size()const{
            return (parent_.off_[i_+1]-parent_.off_[i_])*
                   (parent_.off_[j_+1]-parent_.off_[j_]);
        }
        ///Returns the current basis function pair (indices refer to full basis)
        const array_t& operator*()const{
            return basis_pair_idx_;
        }

        ///Increments the iterator
        const BasisFunctionPairItr& operator++(){return next();}

        //@{
        /** Comparison operators.  No check to ensure basis sets are the same. */
        bool operator==(const BasisFunctionPairItr& other)const{
            return basis_pair_idx_==other.basis_pair_idx_;
        }
        bool operator!=(const BasisFunctionPairItr& other)const
        {
            return basis_pair_idx_!=other.basis_pair_idx_;
        }
        bool operator<(const BasisFunctionPairItr& other)const
        {
            const bool use_2(basis_pair_idx_[0]==other.basis_pair_idx_[0]);
            return basis_pair_idx_[use_2]<other.basis_pair_idx_[use_2];
        }
        bool operator<=(const BasisFunctionPairItr& other)const
        {
            return *this==other || *this<other;
        }
        bool operator >(const BasisFunctionPairItr& other)const
        {
            return !(*this<=other);
        }
        bool operator>=(const BasisFunctionPairItr& other)const
        {
            return !(*this<other);
        }
        //@}
    };

    const ShellPairItr& next()
    {
        ++current_pair_;
        return *this;
   }

public:

    ///Makes an iterator over all the shell pairs in a particular basis set
    ShellPairItr(const pulsar::BasisSet& bs):
        off_([&](){
            std::vector<size_t> rv;
            size_t nfuncs=0;
            for(const auto& shell : bs)
                for(size_t i=0;i<shell.n_general_contractions();++i)
                {
                    rv.push_back(nfuncs);
                    nfuncs+=shell.general_n_functions(i);
                }
            rv.push_back(nfuncs);
            return rv;
        }()),
        nshells_(off_.size()-1),
        shell_pairs_([=](){
            //Number of shell pairs is N choose 2 with repeats
            std::vector<array_t> rv(nshells_*(nshells_+1)/2);
            for(size_t i=0,counter=0;i<nshells_;++i)
                for(size_t j=0;j<=i;++j,++counter)
                    rv[counter]=array_t({i,j});
            return rv;
        }()),
        current_pair_(0)
    {
    }
    ///Default destructor
    ~ShellPairItr()=default;
    ///No default construction 'cause of const members
    ShellPairItr()=delete;
    ///Default copy is fine
    ShellPairItr(const ShellPairItr&)=default;
    ///Can't copy 'cause of const members
    ShellPairItr& operator=(const ShellPairItr&)=delete;
    ///Can't move 'cause of const members
    ShellPairItr(ShellPairItr&&)=delete;

    ///Returns the number of shell pairs this iterator iterates over
    size_t size()const{return shell_pairs_.size();}

    ///Makes the iterator start over (very cheap, no reallocations or computation)
    void reset(){current_pair_=0;}

    ///True if we are done iterating
    bool done()const{return current_pair_==shell_pairs_.size();}

    ///True if this iterator is still good
    operator bool() const{return !done();}

    ///The shell pair being pointed to
    const array_t& operator*()const{return shell_pairs_[current_pair_];}

    ///Allows random access
    const array_t& operator[](size_t i)const{return shell_pairs_[i];}

    ///Type of an iterator over the basis functions in this shell pair
    using const_iterator=BasisFunctionPairItr;

    ///Returns an iterator to beginning of this shell pair
    BasisFunctionPairItr begin()const{
        return BasisFunctionPairItr(*this,false);
    }

    ///Returns an iterator just past the end of this shell pair
    BasisFunctionPairItr end()const{
        return BasisFunctionPairItr(*this,true);
    }

    ///Increments the iterator
    const ShellPairItr& operator++(){return next();}

    //@{
    /** Comparison operators.  All comparisons are lexographic and ignore check
     *  for same basis set
     */
    bool operator==(const ShellPairItr& other)const{
        return current_pair_==other.current_pair_;
    }
    bool operator<(const ShellPairItr& other)const
    {
        return current_pair_<other.current_pair_;
    }
    bool operator!=(const ShellPairItr& other)const
    {
        return !((*this)==other);
    }
    bool operator<=(const ShellPairItr& other)const
    {
        return *this==other || *this<other;
    }
    bool operator >(const ShellPairItr& other)const
    {
        return !(*this<=other);
    }
    bool operator>=(const ShellPairItr& other)const
    {
        return !(*this<other);
    }
    //@}
};

}//End namespace pulsar_scf
