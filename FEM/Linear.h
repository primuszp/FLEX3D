/**
 * @file Linear.h
 * Derived class from Analysis for linear elastic problems.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#ifndef Linear_h
#define Linear_h

#include "Analysis.h"

/* Derived class for solving linear elastic problems.
 */
class Linear : public Analysis
{
  public:
    /* See the documentation of base class Analysis.
     */
    Linear(Mesh & meshInfo); // ctor cannot be inherited, should explicitly call base class's ctor in derived class's ctor
    ~Linear();
    void solve();

};

#endif /* Linear_h */
