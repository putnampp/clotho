#ifndef IGENERATOR_HPP_
#define IGENERATOR_HPP_

template < class Element >
struct igenerator {
    typedef Element result_type;

    virtual const std::string name() const = 0;
    virtual result_type operator()() = 0;
};

template < class URNG, class Element >
struct iran_generator : virtual public igenerator< Element > {

    virtual URNG *  getRNG() = 0;
};

#endif  // IGENERATOR_HPP_
