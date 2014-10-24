#ifndef CLOTHO_CONFIG_H_
#define CLOTHO_CONFIG_H_

#ifdef CLOTHO_UNITTEST
#define CLOTHO_PROTECTED
#else
#define CLOTHO_PROTECTED protected:
#endif

#endif  // CLOTHO_CONFIG_H_
