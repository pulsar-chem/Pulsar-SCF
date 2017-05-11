#pragma once
#include<pulsar/system/BasisSet.hpp>
#include<array>

namespace pulsar_scf {
class ShellPairItr;


/** \brief Class that iterates over the pairs of the shell functions
 *
 *  A BasisFunctionPairItr instance can only be made by a ShellPairItr instance.
 *  When made it will iterate over the the actual basis functions of the shell
 *  pair. E.G., assume our ShellPairItr is on shell pair 1 3, and assume shell
 *  1's basis functions start at index 3 and shell 3's start at index 8 then
 *  this class will return:
 *  ~~~
 *  3,8
 *  3,9
 *  ...
 *  3,m
 *  4,8
 *  4,9
 *  ...
 *  4,m
 *  5,8
 *  ...
 *  n,m
 *  ~~~
 *  where n is the length of shell 1 and m is the length of shell 3.
 *
 */
class BasisFunctionPairItr{
private:
     friend class ShellPairItr;
     using array_t=std::array<size_t,2>;
     array_t idx_;///<The current pair of shells
     const array_t nbf_;///<Last BF number of the shells plus 1
     const array_t& off_;///<BF number starting the shells

     ///Sets idx_ to the next pair of basis functions
     const BasisFunctionPairItr& next()
     {
         if(++idx_[1]==nbf_[1])
             idx_=array_t({++idx_[0],off_[1]});
         return *this;
     }

     ///Makes an iterator to the beginning of the block
     BasisFunctionPairItr(const array_t& off, const array_t& nbf):
        idx_({off[0],off[1]}),
        nbf_({nbf[0]+off[0],nbf[1]+off[1]}),
        off_(off)
     {}

     ///Makes an iterator just past the last shell pair
     BasisFunctionPairItr(const array_t& off, const array_t& nbf,bool):
         idx_({nbf[0]+off[0],off[1]}),
         nbf_({nbf[0]+off[0],nbf[1]+off[1]}),
         off_(off)
     {}

 public:
     ///Returns the current index
     const array_t& operator*()const{return idx_;}

     ///Default copy constructor, shallow copies off_
     BasisFunctionPairItr(const BasisFunctionPairItr&)=default;
     ///Default move constructor
     BasisFunctionPairItr(BasisFunctionPairItr&&)=default;
     ///Because of reference no default constructor
     BasisFunctionPairItr()=delete;
     ///Default destructor is fine
     ~BasisFunctionPairItr()=default;

     ///Returns true if this iterator equals other
     bool operator==(const BasisFunctionPairItr& other)const{
         return idx_==other.idx_ &&
                nbf_==other.nbf_ &&
                &off_==&other.off_;//Use fact that off is a reference
     }

     ///Negates operator==
     bool operator!=(const BasisFunctionPairItr& other)const{
         return !(*this==other);
     }

     ///Prefix increment only (no postfix increment owing to inefficiency)
     const BasisFunctionPairItr& operator++(){return next();}

 };

/** \brief Iterates over all unique pairs of shells
 *
 *   A ShellPairItr is only good so long as the BasisSet it is associated with
 *   is in scope.
 *
 */
class ShellPairItr{
protected:
    using array_t=std::array<size_t,2>;
    const pulsar::BasisSet bs_;
    array_t idx_;
    array_t nbf_;
    array_t off_;
    size_t max_;

    const ShellPairItr& next(){
        const bool is_good(idx_[1]<idx_[0]);
        if(is_good)++idx_[1];
        else idx_=array_t({++idx_[0],0});
        if(done())return *this;
        const size_t newj=idx_[1];
        nbf_[1]=bs_.shell(newj).n_functions();
        off_[1]=bs_.shell_start(newj);
        if(!is_good)
        {
           const size_t newi=idx_[0];
           nbf_[0]=bs_.shell(newi).n_functions();
           off_[0]=bs_.shell_start(newi);
        }
        return *this;
   }

public:
    ShellPairItr(const pulsar::BasisSet& bs):
        bs_(bs),
        idx_({0,0}),
        nbf_({bs.shell(0).n_functions(),bs.shell(0).n_functions()}),
        off_({0,0}),
        max_(bs.n_shell())
    {}

    ~ShellPairItr()=default;
    ShellPairItr()=delete;
    ShellPairItr(const ShellPairItr&)=default;
    ShellPairItr(ShellPairItr&&)=default;

    bool done()const{return idx_[0]==max_;}
    operator bool() const{return !done();}
    const array_t& operator*()const{return idx_;}


    ///Returns an iterator to beginning of this shell pair
    BasisFunctionPairItr begin()const{
        return !done()?BasisFunctionPairItr(off_,nbf_):end();
    }

    ///Returns an iterator just past the end of this shell pair
    BasisFunctionPairItr end()const{
        return BasisFunctionPairItr(off_,nbf_,true);
    }

    const ShellPairItr& operator++(){return next();}

};

}
