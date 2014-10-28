// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology
//
// An implementation of bucked-based kd-tree branching with mean value
// The algorithm taken from:
// Jerome H. Friedman, Jon Bentley and Raphael Finkel
// "An Algorithm for Finding Best Matches in Logarithmic Expected Time"
// ACM Transactions on Mathematical Software, VOL.3, No.3, pp.209-226, 1977

#pragma once

#include "mesh.h"

namespace hbm {

/// KD-tree
/// coordinates are duplicated for acceleration.
class KdTree {
public:
    static const int DEFAULT_BUCKET_SIZE = 15 ;
    static const int DIMENSION  = 3;

    KdTree() = default;
    KdTree(const KdTree&) = delete;
    KdTree& operator=(const KdTree&) = delete;

    void Construct(std::vector<CCustomVertex>::const_iterator begin, std::vector<CCustomVertex>::const_iterator end, int bucket_size = DEFAULT_BUCKET_SIZE);

    /// update kd-tree
    void Reconstruct(int bucket_size = DEFAULT_BUCKET_SIZE);

    /// k-nearest-neighbor search
    /// @return true if found, otherwise false
    bool GetKNearest(const float *query, // query coordinates
                     float max_distance2, // upper bound of squared distance
                     int k,
                     size_t *indices, // index of the nearest neighbor
                     float *distance2s // squared distance to the nearest neighbor; distance2 < max_distance2 if found.
                    ) const;

private:
    // element in sorted array
    struct Data {
        float pos[DIMENSION];
        size_t id; // data ID
    };

    /// node
    struct Node {
        Node *lower = 0; //  lower coordinates; 0 if leaf
        union {
            // node
            struct {
                Node *upper; //  upper coordinates
                float median; // lower bound of the upper subtree
                int split; // dimension to split data
            };
            // leaf
            struct {
                Data *begin; // first data
                Data *end; // end of data
            };
        };
    };

    /// determines the dimension for data splitting
    int SelectDimension(const Data *begin, const Data *end) const;

    /// construct KD-tree
    void BuildTree(Data *begin, Data *end, int bucket_size);

    /// search for the data at the nearest position
    void SearchForKNearest(const Node *node,
                           const float *query, // query coordinates
                           int k,
                           size_t *indices, // iterator to the nearest neighbor
                           float *distance2s // squared distance to the nearest neighbor
                          ) const;

private:
    std::vector<Node> nodes_; // nodes of KD-tree
    std::vector<Data> data_; // input data
};

} // slib
