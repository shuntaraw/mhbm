// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology
//
// An implementation of bucked-based kd-tree branching with mean value
// The algorithm taken from:
// Jerome H. Friedman, Jon Bentley and Raphael Finkel
// "An Algorithm for Finding Best Matches in Logarithmic Expected Time"
// ACM Transactions on Mathematical Software, VOL.3, No.3, pp.209-226, 1977

#pragma once

#include <cassert>
#include <algorithm>
#include <functional>
#include <vector>

namespace slib {

// traits
namespace kdt {
template <typename Point>
struct dimension {
    // static const int value = ...;
};
template <typename Point>
struct coordinate_type {
    // typedef ... type;
};
template <typename Point, int D>
struct access {
    // static typename coordinate_type<Point>::type get(const Point& p) { return ...; }
};
}

/// KD-tree
/// coordinates are duplicated for acceleration.
template <typename Point>
class KdTree {
public:
    static const int DEFAULT_BUCKET_SIZE = 16;
    static const int DIMENSION = kdt::dimension<Point>::value;
    typedef typename kdt::coordinate_type<Point>::type T;

    KdTree() = default;
    KdTree(const KdTree&) = delete;
    KdTree& operator=(const KdTree&) = delete;

    template<typename FwdIt>
    void Construct(FwdIt begin, FwdIt end, int bucket_size = DEFAULT_BUCKET_SIZE) {
        size_t ndata = std::distance(begin, end);
        data_.resize(ndata);
        for (size_t i = 0; i < ndata; i++, ++begin) {
            CopyCoordinate < DIMENSION - 1 > (*begin, data_[i].pos);
            data_[i].id = i;
        }
        Reconstruct(bucket_size);
    }

    /// update kd-tree
    void Reconstruct(int bucket_size = DEFAULT_BUCKET_SIZE) {
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
    bool GetKNearest(const T *query, // query coordinates
                     T max_distance2, // upper bound of squared distance
                     int k,
                     size_t *indices, // index of the nearest neighbor
                     T *distance2s // squared distance to the nearest neighbor; distance2 < max_distance2 if found.
                    ) const {
        std::fill_n(distance2s, k, max_distance2);
        SearchForKNearest(&nodes_.front(), query, k, indices, distance2s);
        return distance2s[0] < max_distance2; // found at least one point?
    }

private:
    // element in sorted array
    struct Data {
        T pos[DIMENSION];
        size_t id; // data ID
    };

    /// node
    struct Node {
        Node *lower = 0; //  lower coordinates; 0 if leaf
        union {
            // node
            struct {
                Node *upper; //  upper coordinates
                T median; // lower bound of the upper subtree
                int split; // dimension to split data
            };
            // leaf
            struct {
                Data *begin; // first data
                Data *end; // end of data
            };
        };
    };

    template <int D>
    void CopyCoordinate(const Point& point, T *pos) {
        pos[D] = kdt::access<Point, D>::get(point);
        CopyCoordinate < D - 1 > (point, pos);
    }

    template <>
    void CopyCoordinate<0>(const Point& point, T *pos) {
        pos[0] = kdt::access<Point, 0>::get(point);
    }

    /// compute the smallest power of two larger than or equal to an integer.
    /// @see http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
    static size_t NextPowerOfTwo(size_t n) {
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
    static size_t NumNodes(size_t ndata, size_t bucket_size) {
        int nbuckets = ndata / bucket_size;
        if (ndata % bucket_size) {
            nbuckets++;
        }
        size_t m = NextPowerOfTwo(nbuckets);
        return std::min(2 * m - 1, 2 * ndata - (bucket_size - 1) * m - 1);
    }

    /// compute a square value
    /// @return v*v
    static T sqr(T v) {
        return v * v;
    }

    template <int D>
    static T ssd(const T *query, const T *data) {
        return sqr(query[D] - data[D]) + ssd < D - 1 > (query, data);
    }

    template <>
    static T ssd<0>(const T *query, const T *data) {
        return sqr(query[0] - data[0]);
    }

    static void InsertResult(size_t vid,
                             T distance2,
                             int k,
                             size_t *indices, // iterator to the nearest neighbor
                             T *distance2s  // squared distance to the nearest neighbor
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

    /// determines the dimension for data splitting
    int SelectDimension(const Data *begin, const Data *end) const {
        T lower[DIMENSION], upper[DIMENSION];
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
        T diff = upper[0] - lower[0];
        for (int dim = 1; dim < DIMENSION; ++dim) {
            T length = upper[dim] - lower[dim];
            if (diff < length) {
                diff = length;
                selected = dim;
            }
        }
        return selected;
    }

    /// construct KD-tree
    void BuildTree(Data *begin, Data *end, int bucket_size) {
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
    void SearchForKNearest(const Node *node,
                           const T *query, // query coordinates
                           int k,
                           size_t *indices, // iterator to the nearest neighbor
                           T *distance2s // squared distance to the nearest neighbor
                          ) const {
        if (!node->lower) { // leaf
            // compare to all
            for (auto data = node->begin; data != node->end; ++data) {
                T dist2 = ssd < DIMENSION - 1 > (query, data->pos);
                if (dist2 < distance2s[k - 1]) {
                    InsertResult(data->id, dist2, k, indices, distance2s);
                }
            }
        } else {
            T distance2split = sqr(node->median - query[node->split]);
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

private:
    std::vector<Node> nodes_; // nodes of KD-tree
    std::vector<Data> data_; // input data
};

} // slib
