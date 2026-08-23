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
    
    /*
    Node *lower_node;
    Node *upper_node;
     */
    
    Node *lower_node = nullptr;
    Node *upper_node = nullptr;
    
    Branch();
    
    // Branch(Node *l, Node *u);
    
    Branch(Node *l, Node *u);

    Branch(const Node_ptr &l, const Node_ptr &u);

    Branch(Node *l, const Node_ptr &u);

    Branch(const Node_ptr &l, Node *u);
    
    double length();
    
    bool operator<(const Branch &other) const;
    
    bool operator==(const Branch &other) const;
    
    bool operator!=(const Branch &other) const;
};

struct branch_hash {
    std::size_t operator()(const Branch& b) const {
        std::hash<const Node *> hasher;
        std::size_t h1 = hasher(b.lower_node);
        std::size_t h2 = hasher(b.upper_node);
        return h1 ^ (h2 << 1);
    }
};

#endif /* Branch_hpp */
