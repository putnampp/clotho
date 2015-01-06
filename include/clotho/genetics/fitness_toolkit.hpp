#ifndef FITNESS_TOOLKIT_HPP_
#define FITNESS_TOOLKIT_HPP_

#include "clotho/genetics/ifitness_generator.hpp"

#include <unordered_map>
#include <memory>
#include <boost/property_tree/ptree.hpp>

#include <ostream>

extern const std::string FITNESS_BLOCK_K;

class fitness_toolkit {
public:

    typedef std::unordered_map< std::string, ifitness_generator * >  generator_map;
    typedef generator_map::iterator                                  generator_iterator;
    typedef generator_map::const_iterator                            generator_citerator;

    static fitness_toolkit * getInstance() {
        static fitness_toolkit inst;
        return &inst;
    }

    std::shared_ptr< ifitness_generator > get_tool( boost::property_tree::ptree & config );

    void tool_configurations( boost::property_tree::ptree & config );

    void register_tool( ifitness_generator * gen );

    friend std::ostream & operator<<( std::ostream & out, const fitness_toolkit & ftk );

    virtual ~fitness_toolkit();

protected:
    fitness_toolkit();

    generator_map   m_tools;
};

std::ostream & operator<<( std::ostream & out, const fitness_toolkit & ftk );

#endif  // FITNESS_TOOLKIT_HPP_
