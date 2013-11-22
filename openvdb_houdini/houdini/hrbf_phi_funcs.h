/*
 * This software is governed by the CeCILL-B license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-B
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info" or the LICENCE.txt file present in this project.
*/

#ifndef HRBF_PHI_FUNCS_HPP_
#define HRBF_PHI_FUNCS_HPP_

/** @brief Radial basis functions definitions (function phi)
  Here you can add more radial basis function definitions.
*/

/**
 * @class Rbf_pow3
 * Radial basis function phi(x) = x^3 first and second derivative
 *
 **/
template<typename Scalar>
struct Rbf_pow3
{
    // phi(x) = x^3
    static inline Scalar f  (const Scalar& x) { return x*x*x;             }
    // first derivative phi'(x) = 3x^2
    static inline Scalar df (const Scalar& x) { return Scalar(3) * x * x; }
    // second derivative phi''(x) = 6x
    static inline Scalar ddf(const Scalar& x) { return Scalar(6) * x;     }
};

// -----------------------------------------------------------------------------

template<typename Scalar>
struct Rbf_x_sqrt_x
{
    static inline Scalar f  (const Scalar& x2) {
        Scalar x = std::sqrt(x2); return x*x*x;
    }

    static inline Scalar df (const Scalar& x2) {
        return Scalar(3)/Scalar(2) * std::sqrt(x2);
    }

    static inline Scalar ddf(const Scalar& x2) {
        return Scalar(3)/Scalar(4) / std::sqrt(x2);
    }
};

// -----------------------------------------------------------------------------

#define ORDER_x_2 2

template<typename Scalar>
struct Rbf_thin_plate
{
    static inline Scalar f (const Scalar& x) {
        return pow( x, ORDER_x_2 ) * log(x);
    }

    static inline Scalar df (const Scalar& x) {
        return pow( x, ORDER_x_2 - Scalar(1) ) * ( ORDER_x_2 * log(x) + Scalar(1) );
    }

    static inline Scalar ddf(const Scalar& x) {
        return pow( x, ORDER_x_2 - Scalar(2) ) * ( ORDER_x_2 * (ORDER_x_2 - Scalar(1)) * log(x)
                                                               + ORDER_x_2 - Scalar(1) + ORDER_x_2 );
    }
};

#endif //HRBF_PHI_FUNCS_HPP_
