//
//  Branch.hpp
//  SINGER
//
//  Created by Yun Deng on 3/31/22.
//

#ifndef Branch_hpp
#define Branch_hpp

#include <stdio.h>
#include "Node.hpp"

using Node_ptr = shared_ptr<Node>;

class Branch {
    
public:

    Node *lower_node = nullptr;
    Node *upper_node = nullptr;

    Branch();

    Branch(Node *l, Node *u);

    Branch(const Node_ptr &l, const Node_ptr &u);

    bool operator<(const Branch &other) const;
    
    bool operator==(const Branch &other) const;
    
    bool operator!=(const Branch &other) const;
};

#endif /* Branch_hpp */
