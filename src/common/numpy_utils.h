#ifndef NUMPY_UTILS_H
#define NUMPY_UTILS_H

#include <string>
#include <vector>
#include <npe.h>


bool assert_shape_and_dtype(const pybind11::array& arr, std::string name, pybind11::dtype dtype,
                            const std::vector<ssize_t>& shape) {
    if (!arr.dtype().is(dtype)) {
        throw pybind11::value_error("Invalid dtype for argument '" + name + "'. Expected '" +
                                    dtype.kind() + "' but got '" + arr.dtype().kind() + "'.");
    }
    if (shape.size() != arr.ndim()) {
        throw pybind11::value_error("Invalid number of dimensions for argument '" + name + "'. Expected " +
                                    std::to_string(shape.size()) + " but got " + std::to_string(arr.ndim()) + ".");
    }
    bool nonempty = true;
    for (int i = 0; i < shape.size(); i++) {
        if (arr.shape()[i] <= 0) {
            nonempty = false;
        }
        if (shape[i] < 0) {
            if (arr.shape()[i] == 0) {
                continue;
            } else if (arr.shape()[i] == -shape[i]) {
                continue;
            }
        } else if (shape[i] == arr.shape()[i]) {
            continue;
        }

        throw pybind11::value_error("Invalid  shape for argument '" + name + "' at dimension " +
                                    std::to_string(i) + ". Expected " + std::to_string(shape[i]) +
                                    " but got " + std::to_string(arr.shape()[i]) + ".");
    }

    return nonempty;
}

#endif // NUMPY_UTILS_H