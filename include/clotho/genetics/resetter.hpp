#ifndef RESETTER_HPP_
#define RESETTER_HPP_

template < class ObjectType >
struct resetter {
    static void reset( ObjectType & obj ) {}
};

template < class ObjectType >
struct resetter< std::shared_ptr< ObjectType > > {
    static void reset( std::shared_ptr< ObjectType > & ptr ) {
        ptr.reset();
    }
};

#endif  // RESETTER_HPP_
