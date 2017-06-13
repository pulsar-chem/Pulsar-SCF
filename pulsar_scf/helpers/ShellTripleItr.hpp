#pragma once
#include "pulsar_scf/helpers/ShellPairItr.hpp"

namespace pulsar_scf {

///This class currently assumes it is iterating over (P|mn)
///Usage and details of this class are nearly identical to ShellPairItr so
///consult its documentation
class ShellTripleItr{
    using array_t=std::array<size_t,3>;
    ///The starting index for the DF shell (last element is the number of basis functions)
    const std::vector<size_t> off_;
    ///The total number of DF shells
    const size_t nshells_;
    ///The index of the current DF shell
    size_t current_DF_;
    ///The current shell ket
    ShellPairItr ket_;

    ///An iterator over the basis functions
    class BasisFunctionTripleItr{
    private:
        friend class ShellTripleItr;
        ///Who made me
        const ShellTripleItr& parent_;
        ///Current shell
        const size_t P_;
        ///Current bf
        size_t bf_;
        ///Just past the last bf for P
        const size_t bf_end_;
        ///The current bfs in the shell pair
        ShellPairItr::const_iterator pair_;
        ///The end of the bfs in the current shell pair
        const ShellPairItr::const_iterator end_;
        ///Increments the iterator
        const BasisFunctionTripleItr& next()
        {
            ++pair_;
            if(pair_!=end_)return *this;
            pair_.reset();
            ++bf_;
            return *this;
        }

        ///Constructs a new BasisFunctionTripleItr, only callable by ShellTripleItr
        BasisFunctionTripleItr(const ShellTripleItr& parent,bool end):
            parent_(parent),
            P_(parent_.current_DF_),
            bf_(parent_.off_[P_+end]),
            bf_end_(parent_.off_[P_+1]),
            pair_(parent_.ket_.begin()),
            end_(parent_.ket_.end())
        {}
    public:
        void reset(){
            pair_.reset();
            bf_=parent_.off_[P_];
        }

        ///Returns the size of the block
        size_t size()const{
            return (parent_.off_[P_+1]-parent_.off_[P_])*
                    pair_.size();
        }
        ///Returns the current basis function pair (indices refer to full basis)
        const array_t operator*()const{
            return {bf_,(*pair_)[0],(*pair_)[1]};
        }

        ///Increments the iterator
        const BasisFunctionTripleItr& operator++(){return next();}

        //@{
        /** Comparison operators.  No check to ensure basis sets are the same. */
        bool operator==(const BasisFunctionTripleItr& other)const{
            return bf_==other.bf_ && pair_==other.pair_;
        }
        bool operator!=(const BasisFunctionTripleItr& other)const
        {
            return !((*this)==other);
        }
        bool operator<(const BasisFunctionTripleItr& other)const
        {
            return std::tie(bf_,pair_)<std::tie(other.bf_,other.pair_);
        }
        bool operator<=(const BasisFunctionTripleItr& other)const
        {
            return *this==other || *this<other;
        }
        bool operator >(const BasisFunctionTripleItr& other)const
        {
            return !(*this<=other);
        }
        bool operator>=(const BasisFunctionTripleItr& other)const
        {
            return !(*this<other);
        }
        //@}
    };

    const ShellTripleItr& next()
    {
        ++ket_;
        if(ket_)return *this;
        ket_.reset();
        ++current_DF_;
        return *this;
   }

public:

    ///Makes an iterator over all the shell pairs in a particular basis set
    ShellTripleItr(const pulsar::BasisSet& dfbs,const pulsar::BasisSet& bs):
        off_([&](){
            std::vector<size_t> rv;
            size_t nfuncs=0;
            for(const auto& shell : dfbs)
                for(size_t i=0;i<shell.n_general_contractions();++i)
                {
                    rv.push_back(nfuncs);
                    nfuncs+=shell.general_n_functions(i);
                }
            rv.push_back(nfuncs);
            return rv;
        }()),
        nshells_(off_.size()-1),
        current_DF_(0),
        ket_(bs)
    {
    }
    ///Default destructor
    ~ShellTripleItr()=default;
    ///No default construction 'cause of const members
    ShellTripleItr()=delete;
    ///Default copy is fine
    ShellTripleItr(const ShellTripleItr&)=default;
    ///Can't copy 'cause of const members
    ShellTripleItr& operator=(const ShellTripleItr&)=delete;
    ///Can't move 'cause of const members
    ShellTripleItr(ShellTripleItr&&)=delete;

    ///Returns the number of shell pairs this iterator iterates over
    size_t size()const{return nshells_*ket_.size();}

    ///Makes the iterator start over
    void reset()
    {
        current_DF_=0;
        ket_.reset();
    }

    ///True if we are done iterating
    bool done()const{return current_DF_==nshells_;}

    ///True if this iterator is still good
    operator bool() const{return !done();}

    ///The shell pair being pointed to
    const array_t operator*()const{return {current_DF_,(*ket_)[0],(*ket_)[1]};}

    ///Type of an iterator over the basis functions in this shell pair
    using const_iterator=BasisFunctionTripleItr;

    ///Returns an iterator to beginning of this shell pair
    BasisFunctionTripleItr begin()const{
        return BasisFunctionTripleItr(*this,false);
    }

    ///Returns an iterator just past the end of this shell pair
    BasisFunctionTripleItr end()const{
        return BasisFunctionTripleItr(*this,true);
    }

    ///Increments the iterator
    const ShellTripleItr& operator++(){return next();}

    //@{
    /** Comparison operators.  All comparisons are lexographic and ignore check
     *  for same basis set
     */
    bool operator==(const ShellTripleItr& other)const{
        return current_DF_==other.current_DF_ && ket_==other.ket_;
    }
    bool operator<(const ShellTripleItr& other)const
    {
        return std::tie(current_DF_,ket_)<std::tie(other.current_DF_,other.ket_);
    }
    bool operator!=(const ShellTripleItr& other)const
    {
        return !((*this)==other);
    }
    bool operator<=(const ShellTripleItr& other)const
    {
        return *this==other || *this<other;
    }
    bool operator >(const ShellTripleItr& other)const
    {
        return !(*this<=other);
    }
    bool operator>=(const ShellTripleItr& other)const
    {
        return !(*this<other);
    }
    //@}
};



}//End namespace pulsar_scf
