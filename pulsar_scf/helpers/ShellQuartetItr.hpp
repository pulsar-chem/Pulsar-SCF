#pragma once
#include "pulsar_scf/helpers/ShellPairItr.hpp"

namespace pulsar_scf{

/** \brief Iterates over all unique quartets of shells for a given basis set.
 *
 *   The usage for this class is similar to that of ShellPairItr so consult its
 *   documentation.
 */
class ShellQuartetItr{
private:
    using array_t=std::array<size_t,4>;
    using shell_quartet_t=std::pair<ShellPairItr,ShellPairItr>;
    ///Current shell quartet
    shell_quartet_t shell_quartet_;

    ///An iterator over the basis functions
    class BasisFunctionQuartetItr{
    private:
        friend class ShellQuartetItr;

        ///Parent ShellPairItr instance
        const ShellQuartetItr& parent_;

        using bf_quartet_t=std::pair<ShellPairItr::const_iterator,
                                     ShellPairItr::const_iterator>;
        ///The actual quartet
        bf_quartet_t bf_quartet_;

        const BasisFunctionQuartetItr& next()
        {
            const auto end=parent_.shell_quartet_.second.end();
            if(++bf_quartet_.second!=end)
                return *this;
            ++bf_quartet_.first;
            bf_quartet_.second.reset();
            return *this;
        }

        ///Only callable by ShellQuartetItr
        BasisFunctionQuartetItr(const ShellQuartetItr& parent,bool end):
            parent_(parent),
            bf_quartet_({!end?parent.shell_quartet_.first.begin():
                              parent.shell_quartet_.first.end()  ,
                         parent.shell_quartet_.second.begin()})
           {}
    public:

        ///Returns the current basis function pair (indices refer to full basis)
        array_t operator*()const{
            const auto& bra=*bf_quartet_.first;
            const auto& ket=*bf_quartet_.second;
            return {bra[0],bra[1],ket[0],ket[1]};
        }

        ///Increments the iterator
        const BasisFunctionQuartetItr& operator++(){return next();}

        bool operator==(const BasisFunctionQuartetItr& other)const{
            return bf_quartet_==other.bf_quartet_;
        }
        bool operator!=(const BasisFunctionQuartetItr& other)const
        {
            return bf_quartet_!=other.bf_quartet_;
        }
        bool operator<(const BasisFunctionQuartetItr& other)const
        {
            return bf_quartet_<other.bf_quartet_;
        }
        bool operator<=(const BasisFunctionQuartetItr& other)const
        {
            return *this==other || *this<other;
        }
        bool operator >(const BasisFunctionQuartetItr& other)const
        {
            return !(*this<=other);
        }
        bool operator>=(const BasisFunctionQuartetItr& other)const
        {
            return !(*this<other);
        }

    };

    ///Increments the iterator
    const ShellQuartetItr& next()
    {
        if(++shell_quartet_.second<=shell_quartet_.first)return *this;
        ++shell_quartet_.first;
        shell_quartet_.second.reset();
        return *this;
   }

public:

    ShellQuartetItr(const pulsar::BasisSet& bs):
        shell_quartet_(ShellPairItr(bs),ShellPairItr(bs))
    {
    }

    ~ShellQuartetItr()=default;
    ShellQuartetItr()=delete;
    ShellQuartetItr(const ShellQuartetItr&)=default;
    ShellQuartetItr& operator=(const ShellQuartetItr&)=default;
    ShellQuartetItr(ShellQuartetItr&&)=default;

    ///True if we are done iterating
    bool done()const{return shell_quartet_.first.done();}

    ///True if this iterator is still good
    operator bool() const{return !done();}

    ///The shell pair being pointed to
    array_t operator*()const{
        const auto& bra=*shell_quartet_.first;
        const auto& ket=*shell_quartet_.second;
        return {bra[0],bra[1],ket[0],ket[1]};
    }

    ///Returns an iterator to beginning of this shell pair
    BasisFunctionQuartetItr begin()const{
        return BasisFunctionQuartetItr(*this,false);
    }

    ///Returns an iterator just past the end of this shell pair
    BasisFunctionQuartetItr end()const{
        return BasisFunctionQuartetItr(*this,true);
    }

    ///Increments the iterator
    const ShellQuartetItr& operator++(){return next();}

    ///The degenercy of the current shell quartet
    double degeneracy()const{
        const auto& shell=*(*this);
        const double ij_deg(shell[0]==shell[1]?1.0:2.0);
        const double kl_deg(shell[2]==shell[3]?1.0:2.0);
        const double ij_kl_deg(shell[0]==shell[2] && shell[1]==shell[3]?1.0:2.0);
        return ij_deg*kl_deg*ij_kl_deg;
    }


    //@{
    /** Comparison operators.  All comparisons are lexographic and ignore check
     *  for same basis set
     */
    bool operator==(const ShellQuartetItr& other)const{
        return shell_quartet_==other.shell_quartet_;
    }
    bool operator!=(const ShellQuartetItr& other)const
    {
        return shell_quartet_!=other.shell_quartet_;
    }
    bool operator<(const ShellQuartetItr& other)const
    {
        return shell_quartet_<other.shell_quartet_;
    }
    bool operator<=(const ShellQuartetItr& other)const
    {
        return *this==other || *this<other;
    }
    bool operator >(const ShellQuartetItr& other)const
    {
        return !(*this<=other);
    }
    bool operator>=(const ShellQuartetItr& other)const
    {
        return !(*this<other);
    }
    //@}
};
}//End namespace Pulsar
