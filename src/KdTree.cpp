// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology
//
// An implementation of bucked-based kd-tree branching with mean value
// The algorithm taken from:
// Jerome H. Friedman, Jon Bentley and Raphael Finkel
// "An Algorithm for Finding Best Matches in Logarithmic Expected Time"
// ACM Transactions on Mathematical Software, VOL.3, No.3, pp.209-226, 1977

#include "precompile.h"

#include "KdTree.h"


/// compute the smallest power of two larger than or equal to an integer.
/// @see http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
size_t  NextPowerOfTwo(size_t n) {
    --n;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    if (sizeof(size_t) > 1) {
        n |= n >> 8;    // >= 16bit
    }
    if (sizeof(size_t) > 2) {
        n |= n >> 16;    // >= 32bit
    }
    if (sizeof(size_t) > 4) {
        n |= n >> 32;    // >= 64bit
    }
    return ++n;
}

/// the number of all nodes.
/// @see http://programmizm.sourceforge.net/blog/2011/a-practical-implementation-of-kd-trees
size_t NumNodes(size_t ndata, size_t bucket_size) {
    int nbuckets = ndata / bucket_size;
    if (ndata % bucket_size) {
        nbuckets++;
    }
    size_t m = NextPowerOfTwo(nbuckets);
    return std::min(2 * m - 1, 2 * ndata - (bucket_size - 1) * m - 1);
}

/// compute a square value
/// @return v*v
float sqr(float v) {
    return v * v;
}


float ssd(const float *query, const float *data) {
    return sqr(query[0] - data[0]) + sqr(query[1] - data[1]) + sqr(query[2] - data[2]);
}




void InsertResult(size_t vid,
                  float distance2,
                  int k,
                  size_t *indices, // iterator to the nearest neighbor
                  float *distance2s  // squared distance to the nearest neighbor
                 )  {
    for (int i = 0; i < k; i++) {
        if (distance2s[i] > distance2) {
            for (int j = k - 1; j > i; j--) {
                indices[j] = indices[j - 1];
                distance2s[j] = distance2s[j - 1];
            }
            indices[i] = vid;
            distance2s[i] = distance2;
            return;
        }
    }
}

namespace hbm {
void KdTree::Construct(std::vector<CCustomVertex>::const_iterator  begin, std::vector<CCustomVertex>::const_iterator end, int bucket_size /*= DEFAULT_BUCKET_SIZE*/) {
    size_t ndata = std::distance(begin, end);
    data_.resize(ndata);
    for (size_t i = 0; i < ndata; i++, ++begin) {
        //   CopyCoordinate < DIMENSION - 1 > (*begin, data_[i].pos);
        data_[i].pos[0] = begin->position[0];
        data_[i].pos[1] = begin->position[1];
        data_[i].pos[2] = begin->position[2];
        data_[i].id = i;
    }
    Reconstruct(bucket_size);
}

/// update kd-tree
void KdTree::Reconstruct(int bucket_size /*= DEFAULT_BUCKET_SIZE*/) {
    size_t ndata = data_.size();
    size_t nnodes = NumNodes(ndata, bucket_size);
    nodes_.clear();
    nodes_.reserve(nnodes);
    nodes_.push_back(Node()); // root node
    BuildTree(&data_.front(), &data_.back() + 1, bucket_size);
    assert(nnodes == nodes_.size());
}

/// k-nearest-neighbor search
/// @return true if found, otherwise false
bool KdTree::GetKNearest(const float *query, // query coordinates
                         float max_distance2, // upper bound of squared distance
                         int k,
                         size_t *indices, // index of the nearest neighbor
                         float *distance2s // squared distance to the nearest neighbor; distance2 < max_distance2 if found.
                        ) const {
    std::fill_n(distance2s, k, max_distance2);
    SearchForKNearest(&nodes_.front(), query, k, indices, distance2s);
    return distance2s[0] < max_distance2; // found at least one point?
}






/// determines the dimension for data splitting
int KdTree::SelectDimension(const Data *begin, const Data *end) const {
    float lower[DIMENSION], upper[DIMENSION];
    auto data = begin->pos;
    for (int dim = 0; dim < DIMENSION; ++dim) {
        lower[dim] = data[dim];
        upper[dim] = data[dim];
    }
    for (begin++; begin != end; begin++) {
        auto data = begin->pos;
        for (int dim = 0; dim < DIMENSION; dim++) {
            lower[dim] = std::min(lower[dim], data[dim]);
            upper[dim] = std::max(upper[dim], data[dim]);
        }
    }
    int selected = 0;
    float diff = upper[0] - lower[0];
    for (int dim = 1; dim < DIMENSION; ++dim) {
        float length = upper[dim] - lower[dim];
        if (diff < length) {
            diff = length;
            selected = dim;
        }
    }
    return selected;
}

/// construct KD-tree
void KdTree::BuildTree(Data *begin, Data *end, int bucket_size) {
    size_t ndata = std::distance(begin, end);
    if (ndata <= bucket_size) {
        nodes_.back().begin = begin;
        nodes_.back().end = end;
    } else {
        Node& current = nodes_.back();
        int split_dim = SelectDimension(begin, end);
        current.split = split_dim;
        Data *median = begin + ndata / 2;
        std::nth_element(
            begin,
            median,
            end,
        [&](const Data & it1, const Data & it2) {
            return it1.pos[split_dim] < it2.pos[split_dim];
        });
        current.median = median->pos[split_dim];

        // left
        nodes_.push_back(Node());
        current.lower = &nodes_.back();
        BuildTree(begin, median, bucket_size);

        // right
        nodes_.push_back(Node());
        current.upper = &nodes_.back();
        BuildTree(median, end, bucket_size);
    }
}

/// search for the data at the nearest position
void KdTree::SearchForKNearest(const KdTree::Node *node,
                               const float *query, // query coordinates
                               int k,
                               size_t *indices, // iterator to the nearest neighbor
                               float *distance2s // squared distance to the nearest neighbor
                              ) const {
    if (!node->lower) { // leaf
        // compare to all
        for (auto data = node->begin; data != node->end; ++data) {
            float dist2 = ssd(query, data->pos);
            if (dist2 < distance2s[k - 1]) {
                InsertResult(data->id, dist2, k, indices, distance2s);
            }
        }
    } else {
        float distance2split = sqr(node->median - query[node->split]);
        if (query[node->split] < node->median) {
            SearchForKNearest(node->lower, query, k, indices, distance2s);
            if (distance2split < distance2s[k - 1]) {
                SearchForKNearest(node->upper, query, k, indices, distance2s);
            }
        } else {
            SearchForKNearest(node->upper, query, k, indices, distance2s);
            if (distance2split < distance2s[k - 1]) {
                SearchForKNearest(node->lower, query, k, indices, distance2s);
            }
        }
    }
}


} // slib
