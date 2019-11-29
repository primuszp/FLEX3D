/**
 * @file ElementB20.h
 * Derived class from Element for the isoparametric B20 element.
 *
 * @author Haohang Huang
 * @date May 19, 2019
 */

#ifndef ElementB20_h
#define ElementB20_h

#include "Element.h"

/* Derived class for the isoparametric B20 element.
 */
class ElementB20 : public Element
{
    public:
        /* See the documentation of base class Element. */
        ElementB20();
        ElementB20(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material);
        ~ElementB20();

        Shape* shape() const;

    private:
        /** A static structure that manages all the static members used in this class */
        static staticMembers statics;

};

#endif /* ElementB20_h */
