#include "kernel.h"

namespace kernel {
Vector iteration(const Matrix& kernel, const Vector& start, double accuracy,
        facilities::OnOffStream status)
{
    Vector previous{start};
    Vector next{start + kernel*start};
    unsigned count{1};
    status<<count<<'\n';
    while (max_distance(previous,next)>accuracy) {
        previous = next;
        next = start + kernel*next;
        status<<++count<<'\n';
    }
    status<<"terminated\n";
    return next;
}

Vector inverse(const Matrix& kernel, const Vector& start)
{
    const auto n{kernel.rows()};
    const Matrix identity{Matrix::Identity(n,n)};
    return (identity-kernel).partialPivLu().solve(start);
}
} // kernel
