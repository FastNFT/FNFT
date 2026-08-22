/*
 * This file is part of FNFT.
 *
 * FNFT is free software; you can redistribute it and/or
 * modify it under the terms of the version 2 of the GNU General
 * Public License as published by the Free Software Foundation.
 *
 * FNFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 * Sander Wahls (TU Delft) 2017-2018.
 * Shrinivas Chimmalgi (TU Delft) 2017.
 * Peter J. Prins (TU Delft) 2018, 2020-2021.
 * Igor Chekhovskoy 2026.
 * Irina Vaseva 2026.
 */

/**
 * @brief Properties of the discretizations for the Korteweg-de Vries equation.
 *
 * @file fnft__akns_discretization.h
 * @ingroup akns
 */
#ifndef FNFT__AKNS_DISCRETIZATION_H
#define FNFT__AKNS_DISCRETIZATION_H

#include "fnft__akns_discretization_t.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"


/**
 * @brief This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 *
 * It defines the step size of the frequency grid
 * \f$z = \text{e}^{2j\xi\epsilon_t d}\f$ based on the discretization type.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft_kdv_discretization_t \endlink.
 * @returns polynomial degree, or 0 for discretizations not supported by
 * \link fnft__akns_fscatter \endlink.
 *
 * @ingroup akns
 */
FNFT_UINT fnft__akns_discretization_degree(fnft__akns_discretization_t
        discretization);

/**
 * @brief This routine returns the boundary coefficient based on the
 * discretization.
 *
 * The boundary coefficient is the fraction of the step size that a discretized
 * potential extends beyond the last sample. This routine returns this value
 * based on the discretization of type \link fnft__akns_discretization_t \endlink.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft__akns_discretization_t \endlink.
 * @returns the boundary coefficient, or NAN for discretizations not supported
 * by \link fnft__akns_fscatter \endlink.
 *
 * @ingroup akns
 */
FNFT_REAL fnft__akns_discretization_boundary_coeff(fnft__akns_discretization_t discretization);

/**
 * @brief This routine returns the scaling for effective number of samples based on the
 * discretization.
 *
 * Higher order methods use more than one sample per integration step. This routine returns
 * the value upsampling_factor based on the discretization of type \link fnft__akns_discretization_t \endlink.
 * D_effective = upsampling_factor * D.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft__akns_discretization_t \endlink.
 * @returns the upsampling_factor value, or 0 for unknown discretizations.
 *
 * @ingroup akns
 */
FNFT_UINT fnft__akns_discretization_upsampling_factor(fnft__akns_discretization_t discretization);

/**
 * @brief This routine returns the order of the method based on the
 * discretization.
 *
 * Different numerical methods have different orders of accuray. This routine returns
 * the order of the order based on the discretization of type \link fnft__akns_discretization_t \endlink.
 * When the step-size of the signal samples is reduced by a factor \f$s\f$, the error in the
 * computed values is expected to decrease by a factor \f$s^{order}\f$.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft__akns_discretization_t \endlink.
 * @returns the method_order value, or 0 for unknown discretization.
 *
 * @ingroup akns
 */
FNFT_UINT fnft__akns_discretization_method_order(fnft__akns_discretization_t discretization);

/**
 * @brief This routine maps \f$\lambda\f$ from continuous-time domain to
 * \f$z\f$ in the discrete-time domain based on the discretization. 
 * 
 * This routine maps continuous-time domain value \f$\lambda\f$ to discrete-time domain value
 * \f$z = e^{2j \lambda \epsilon_t degree1step}\f$, where degree1step is based on the discretization 
 * of type \link fnft__akns_discretization_t \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued continuous-time domain spectral parameter \f$\lambda\f$. The values are replaced with
 * discrete-time domain values \f$z\f$.
 * @param[in] discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_lambda_to_z(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft__akns_discretization_t discretization);

/**
 * @brief This routine maps \f$z\f$ from the discrete-time domain to
 * \f$\lambda\f$ in the continuous-time domain based on the discretization. 
 * 
 * This routine maps discrete-time domain value \f$z\f$ to continuous-time domain value
 * \f$\lambda = degree1step\log(z)/(2j\epsilon_t)\f$, where degree1step is based on the discretization 
 * of type \link fnft__akns_discretization_t \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued discrete-time domain spectral parameter \f$z\f$. The values are replaced with
 * continuous-time domain values \f$\lambda\f$.
 * @param[in] discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_z_to_lambda(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft__akns_discretization_t discretization);

/**
 * @brief This routine computes various weights required by some methods
 * based on the discretization. 
 * 
 * This routing computes the special weights required for the 
 * higher-order methods CF\f$^{[4]}_2\f$, CF\f$^{[4]}_3\f$, CF\f$^{[5]}_3\f$ 
 * and CF\f$^{[6]}_4\f$. The weights are used in \link fnft__nse_discretization_preprocess_signal \endlink,
 * \link fnft__akns_scatter_matrix \endlink and \link fnft__nse_scatter_bound_states \endlink.
 * The weights for CF\f$^{[4]}_3\f$ are taken from Alvermann and Fehske (<a href="https://doi.org/10.1016/j.jcp.2011.04.006">Journal of Computational Phys. 230, 2011</a>)
 * and the weights for the others are from Blanes, Casas and Thalhammer(<a href="https://doi.org/10.1016/j.cpc.2017.07.016">Computer Phys. Comm. 220, 2017</a>).
 * The weights are mentioned as matrices in the references. This routine returns 
 * them in row-major order.
 * @param[in,out] qr_weights_ptr Pointer to the starting location of potential weights.
 * @param[in,out] eps_t_weights_ptr Pointer to the starting location of step size weights.
 * @param[in] akns_discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_method_weights(FNFT_COMPLEX **qr_weights_ptr,
                                                  FNFT_COMPLEX **eps_t_weights_ptr,
                                                  fnft__akns_discretization_t const akns_discretization);

/**
 * Scaled five-point finite differences used by the sixth-order exponential
 * methods. The entries are eps_t times the potential and its dimensionless
 * first through fourth differences from Eqs. 13--18 of
 * https://doi.org/10.1016/j.jcp.2021.110764.
 */
typedef struct {
    FNFT_COMPLEX value;
    FNFT_COMPLEX first;
    FNFT_COMPLEX second;
    FNFT_COMPLEX first_low;
    FNFT_COMPLEX second_low;
    FNFT_COMPLEX third;
    FNFT_COMPLEX fourth;
} fnft__akns_es6_stencil_t;

static inline void fnft__akns_discretization_es6_stencil(
        FNFT_COMPLEX const samples[5], FNFT_REAL const eps_t,
        fnft__akns_es6_stencil_t * const stencil)
{
    const FNFT_COMPLEX qm2 = samples[0];
    const FNFT_COMPLEX qm1 = samples[1];
    const FNFT_COMPLEX q = samples[2];
    const FNFT_COMPLEX qp1 = samples[3];
    const FNFT_COMPLEX qp2 = samples[4];

    stencil->value = eps_t*q;
    stencil->first = eps_t*(-qp2+8.0*qp1-8.0*qm1+qm2)/12.0;
    stencil->second = eps_t*(-qp2+16.0*qp1-30.0*q
            +16.0*qm1-qm2)/12.0;
    stencil->first_low = eps_t*(qp1-qm1)/2.0;
    stencil->second_low = eps_t*(qp1-2.0*q+qm1);
    stencil->third = eps_t*(qp2-2.0*qp1+2.0*qm1-qm2)/2.0;
    stencil->fourth = eps_t*(qp2-4.0*qp1+6.0*q-4.0*qm1+qm2);
}

typedef struct {
    FNFT_COMPLEX value;
    FNFT_COMPLEX first;
    FNFT_COMPLEX second;
    FNFT_COMPLEX third;
    FNFT_COMPLEX fourth;
    FNFT_COMPLEX fifth;
    FNFT_COMPLEX sixth;
} fnft__akns_es8_stencil_t;

/* Seven-point central differences from Table 1 of arXiv:2608.11892v1.
 * Missing samples at the signal boundary are supplied as zero by the caller.
 * As in the paper, every result includes one factor eps_t and derivatives are
 * with respect to the dimensionless grid index. */
static inline void fnft__akns_es8_stencil(
        FNFT_COMPLEX const samples[7], FNFT_REAL const eps_t,
        fnft__akns_es8_stencil_t * const stencil)
{
    const FNFT_COMPLEX qm3 = samples[0], qm2 = samples[1];
    const FNFT_COMPLEX qm1 = samples[2], q = samples[3];
    const FNFT_COMPLEX qp1 = samples[4], qp2 = samples[5];
    const FNFT_COMPLEX qp3 = samples[6];

    stencil->value = eps_t*q;
    stencil->first = eps_t*(-qm3/60.0+3.0*qm2/20.0-3.0*qm1/4.0
            +3.0*qp1/4.0-3.0*qp2/20.0+qp3/60.0);
    stencil->second = eps_t*(qm3/90.0-3.0*qm2/20.0+3.0*qm1/2.0
            -49.0*q/18.0+3.0*qp1/2.0-3.0*qp2/20.0+qp3/90.0);
    stencil->third = eps_t*(qm3/8.0-qm2+13.0*qm1/8.0
            -13.0*qp1/8.0+qp2-qp3/8.0);
    stencil->fourth = eps_t*(-qm3/6.0+2.0*qm2-13.0*qm1/2.0
            +28.0*q/3.0-13.0*qp1/2.0+2.0*qp2-qp3/6.0);
    stencil->fifth = eps_t*(-qm3/2.0+2.0*qm2-5.0*qm1/2.0
            +5.0*qp1/2.0-2.0*qp2+qp3/2.0);
    stencil->sixth = eps_t*(qm3-6.0*qm2+15.0*qm1-20.0*q
            +15.0*qp1-6.0*qp2+qp3);
}

static inline void fnft__akns_es8_set_pauli_coefficient(
        FNFT_COMPLEX coeff[24], FNFT_UINT const k,
        FNFT_COMPLEX const a1, FNFT_COMPLEX const a2,
        FNFT_COMPLEX const a3)
{
    coeff[4*k] = a3;
    coeff[4*k+1] = a1-I*a2;
    coeff[4*k+2] = a1+I*a2;
    coeff[4*k+3] = -a3;
}

/* Coefficients of the degree-five matrix polynomial Z(eps_t*lambda),
 * Eqs. 51--60 of arXiv:2608.11892v1. q[0] and r[0] are eps_t times the
 * potentials; the remaining entries are the scaled differences above. */
static inline void fnft__akns_es8_z_coefficients(
        FNFT_COMPLEX const qv[7], FNFT_COMPLEX const rv[7],
        FNFT_COMPLEX coeff[24])
{
    const FNFT_COMPLEX q=qv[0], q1=qv[1], q2=qv[2], q3=qv[3];
    const FNFT_COMPLEX q4=qv[4], q5=qv[5], q6=qv[6];
    const FNFT_COMPLEX r=rv[0], r1=rv[1], r2=rv[2], r3=rv[3];
    const FNFT_COMPLEX r4=rv[4], r5=rv[5], r6=rv[6];
    const FNFT_COMPLEX qr=q*r, qp2=q*q, qp3=qp2*q, rp2=r*r;
    const FNFT_COMPLEX rp3=rp2*r, q1p2=q1*q1, q2p2=q2*q2;
    const FNFT_COMPLEX r1p2=r1*r1, r2p2=r2*r2;

    const FNFT_COMPLEX a1c0=(504.0*q4+3.0*q6+967680.0*r-8064.0*q1p2*r+80.0*q2p2*r-384.0*q1*q3*r-48.0*q4*rp2
        +128.0*q1p2*rp3+8064.0*q1*r*r1+144.0*q3*r*r1+40320.0*r2-144.0*q1p2*r2+144.0*q1*r1*r2
        +128.0*qp3*(r1p2+2.0*r*r2)-16.0*q2*(-2520.0+168.0*rp2-9.0*q1*r1+9.0*r1p2+5.0*r*r2)
        +240.0*q1*r*r3+504.0*r4-16.0*qp2*(16.0*q2*rp2+64.0*q1*r*r1-56.0*r*r1p2+168.0*r2+16.0*rp2*r2+3.0*r4)
        +16.0*q*(60480.0+3.0*q4*r+56.0*q1p2*rp2+504.0*q1*r1+15.0*q3*r1-64.0*q1*rp2*r1-504.0*r1p2
        +q2*(168.0*r+16.0*rp3-5.0*r2)+168.0*r*r2+5.0*r2p2+9.0*q1*r3-24.0*r1*r3+3.0*r*r4)+3.0*r6)/1935360.0;
    const FNFT_COMPLEX a1c1=I*(9.0*q5-48.0*q3*(-21.0+qr)-40320.0*r1+288.0*q1p2*r1+64.0*q*q2*r1
        +2688.0*qr*r1+272.0*q2*r*r1-256.0*qp2*rp2*r1+208.0*q*r1*r2
        +16.0*q1*(2520.0-168.0*qr-13.0*q2*r+16.0*qp2*rp2-18.0*r1p2-17.0*q*r2-4.0*r*r2)
        -1008.0*r3+48.0*qr*r3-9.0*r5)/483840.0;
    const FNFT_COMPLEX a1c2=(-8.0*qp2*r2+40.0*q*q1*r1+8.0*q2*(3.0*qr-rp2-21.0)+24.0*qr*r2
        -24.0*q*r1p2-24.0*q1p2*r+40.0*q1*r*r1-3.0*q4-168.0*r2-3.0*r4)/60480.0;
    const FNFT_COMPLEX a1c3=I*(3.0*q3+8.0*q1*(21.0-4.0*qr)+8.0*(-21.0+4.0*qr)*r1-3.0*r3)/30240.0;
    const FNFT_COMPLEX a1c4=-(q2+r2)/3780.0;
    const FNFT_COMPLEX a1c5=I*(q1-r1)/1890.0;

    const FNFT_COMPLEX a2c0=I*(128.0*qp3*(2.0*r*r2+r1p2)-16.0*qp2*(64.0*q1*r*r1+16.0*q2*rp2-16.0*rp2*r2+56.0*r*r1p2+168.0*r2+3.0*r4)
        +16.0*q*(56.0*q1p2*rp2+64.0*q1*rp2*r1+504.0*q1*r1+9.0*q1*r3+q2*(-16.0*rp3+168.0*r-5.0*r2)
        +15.0*q3*r1+3.0*q4*r-168.0*r*r2-3.0*r*r4+504.0*r1p2+24.0*r1*r3-5.0*r2p2+60480.0)
        -128.0*q1p2*rp3-8064.0*q1p2*r-144.0*q1p2*r2+16.0*q2*(9.0*(q1*r1+r1p2+280.0)+168.0*rp2+5.0*r*r2)
        -384.0*q1*q3*r-8064.0*q1*r*r1-240.0*q1*r*r3-144.0*q1*r1*r2+80.0*q2p2*r-144.0*q3*r*r1
        +48.0*q4*rp2+504.0*q4+3.0*q6-967680.0*r-40320.0*r2-504.0*r4-3.0*r6)/1935360.0;
    const FNFT_COMPLEX a2c1=(-9.0*q5+48.0*q3*(-21.0+qr)-40320.0*r1-288.0*q1p2*r1-64.0*q*q2*r1
        +2688.0*qr*r1+272.0*q2*r*r1-256.0*qp2*rp2*r1+208.0*q*r1*r2
        -16.0*q1*(2520.0-168.0*qr-13.0*q2*r+16.0*qp2*rp2+18.0*r1p2-17.0*q*r2+4.0*r*r2)
        -1008.0*r3+48.0*qr*r3-9.0*r5)/483840.0;
    const FNFT_COMPLEX a2c2=-I*(3.0*q4+24.0*q1p2*r-8.0*q2*(-21.0+3.0*qr+rp2)-40.0*q*q1*r1
        +40.0*q1*r*r1-24.0*q*r1p2-168.0*r2+8.0*qp2*r2+24.0*qr*r2-3.0*r4)/60480.0;
    const FNFT_COMPLEX a2c3=(-3.0*q3+8.0*q1*(-21.0+4.0*qr)-168.0*r1+32.0*qr*r1-3.0*r3)/30240.0;
    const FNFT_COMPLEX a2c4=I*(r2-q2)/3780.0;
    const FNFT_COMPLEX a2c5=-(q1+r1)/1890.0;

    const FNFT_COMPLEX a3c0=(-256.0*qp3*rp2*r1+q1*(256.0*qp2*rp3-16.0*rp2*(168.0*q+13.0*q2)
        -160.0*r*(q*r2-252.0)+9.0*(-32.0*q*r1p2+112.0*r2+r4))+2688.0*qp2*r*r1+48.0*qp2*r*r3
        +208.0*qp2*r1*r2+160.0*q*q2*r*r1-6.0*q3*(8.0*q*rp2-168.0*r-5.0*r2)-40320.0*q*r1
        -1008.0*q*r3-9.0*q*r5+288.0*q1p2*r*r1-1008.0*q2*r1-30.0*q2*r3-9.0*q4*r1+9.0*q5*r)/483840.0;
    const FNFT_COMPLEX a3c1=-I*(60480.0-3.0*q4*r+8.0*q1p2*rp2+1008.0*q1*r1+24.0*q3*r1
        -112.0*q*q1*r*r1+8.0*qp2*r1p2+2.0*q2*(-84.0*r+8.0*q*rp2-5.0*r2)-168.0*q*r2
        +16.0*qp2*r*r2+24.0*q1*r3-3.0*q*r4)/60480.0;
    const FNFT_COMPLEX a3c2=(3.0*q3*r-168.0*q*r1+11.0*q2*r1+32.0*qp2*r*r1
        +q1*(168.0*r-32.0*q*rp2-11.0*r2)-3.0*q*r3)/30240.0;
    const FNFT_COMPLEX a3c3=I*(q2*r-8.0*q1*r1+q*r2)/3780.0;
    const FNFT_COMPLEX a3c4=(q1*r-q*r1)/1890.0;

    fnft__akns_es8_set_pauli_coefficient(coeff,0,a1c0,a2c0,a3c0);
    fnft__akns_es8_set_pauli_coefficient(coeff,1,a1c1,a2c1,a3c1);
    fnft__akns_es8_set_pauli_coefficient(coeff,2,a1c2,a2c2,a3c2);
    fnft__akns_es8_set_pauli_coefficient(coeff,3,a1c3,a2c3,a3c3);
    fnft__akns_es8_set_pauli_coefficient(coeff,4,a1c4,a2c4,a3c4);
    fnft__akns_es8_set_pauli_coefficient(coeff,5,a1c5,a2c5,0.0);
}

/**
 * @brief  This routine preprocesses the signal by resampling and subsampling based on the discretization.
 * The preprocessing is necessary for higher-order methods.
 *
 * This routine preprocess q to generate q_preprocessed and r_preprocessed
 * based on the discretization. The preprocessing may involve resampling
 * and sub-sampling.
 * The routine is based on the following papers:
 *      - Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot; Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019.
 *      - Medvedev, Vaseva, Chekhovskoy and  Fedoruk, <a href="https://doi.org/10.1364/OE.377140">&quot; Exponential fourth order schemes for direct Zakharov-Shabat problem,&quot;</a> Optics Express, vol. 28, pp. 20--39, 2020.
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] r_from_q A function pointer array of length 3. `r_from_q[0](q)` has to return the value of an r-sample given the corresponding q-sample, before preprocessing. `r_from_q[1](q)` and `r_from_q[2](q)` have to return respectively the sample of the first/second derivative of r given the corresponding sample of the first/second derivative of q.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[out] q_preprocessed_ptr Pointer to the starting location of preprocessed signal q_preprocessed.
 * @param[out] r_preprocessed_ptr Pointer to the starting location of preprocessed signal r_preprocessed.
 * @param[in,out] Dsub_ptr Pointer to number of processed samples. Upon entry, *Dsub_ptr
 *             should contain a desired number of samples. Upon exit, *Dsub_ptr
 *             has been overwritten with the actual number of samples that the
 *             routine has chosen. It is usually close to the desired one.
 * @param[out] first_last_index Vector of length two. Upon exit, it contains
 *             the original index of the first and the last sample used to build
 *             q_preprocessed.
 * @param[in] discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
*/
FNFT_INT fnft__akns_discretization_preprocess_signal(FNFT_UINT const D,
                                                     FNFT_COMPLEX const * const q,
                                                     FNFT_COMPLEX (*r_from_q[3])(FNFT_COMPLEX const),
                                                     FNFT_REAL const eps_t,
                                                     FNFT_UINT * const Dsub_ptr,
                                                     FNFT_COMPLEX **q_preprocessed_ptr,
                                                     FNFT_COMPLEX **r_preprocessed_ptr,
                                                     FNFT_UINT * const first_last_index,
                                                     fnft__akns_discretization_t discretization);

/**
 * @brief This routine returns the change of basis matrix from the basis of the discretization to S.
 * @param[out] T 2x2 or 4x4 matrix. Left multiplication of a vector in the basis of the discretization by T changes it to the equivalent vector in S basis.
 * @param[in] xi spectral parameter \f$ \xi \f$.
 * @param[in] derivative_flag When 0, T is 2x2. When 1 T is 4x4, to include the derivatives.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in] akns_discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @param[in] vanilla_flag For calculations for the KdV equation, pass 1 for the original mapping to the AKNS framework with r=-1. Pass 0 for the alternative mapping with q=-1. Unused for NSE.
 * @param[in] PDE PDE of type \link fnft__akns_pde_t \endlink.
 *
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_change_of_basis_matrix_to_S(FNFT_COMPLEX * const T,
                                                               FNFT_COMPLEX const xi,
                                                               FNFT_UINT  const derivative_flag, // 0- > 2x2, 1->4x4
                                                               FNFT_REAL const eps_t,
                                                               fnft__akns_discretization_t const akns_discretization,
                                                               FNFT_UINT vanilla_flag,
                                                               fnft__akns_pde_t const PDE);

/**
 * @brief This routine returns the change of basis matrix from the S basis to the basis of the discretization.
 * @param[out] T 2x2 or 4x4 matrix. Left multiplication of a vector in the S basis of the discretization by T changes it to the equivalent vector in the basis of the discretization.
 * @param[in] xi spectral parameter \f$ \xi \f$.
 * @param[in] derivative_flag When 0, T is 2x2. When 1 T is 4x4, to include the derivatives.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in] akns_discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @param[in] vanilla_flag For calculations for the KdV equation, pass 1 for the original mapping to the AKNS framework with r=-1. Pass 0 for the alternative mapping with q=-1. Unused for NSE.
 * @param[in] PDE PDE of type \link fnft__akns_pde_t \endlink.
 * 
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_change_of_basis_matrix_from_S(FNFT_COMPLEX * const T,
                                                                FNFT_COMPLEX const xi,
                                                                FNFT_UINT  const derivative_flag, // 0- > 2x2, 1->4x4
                                                                FNFT_REAL const eps_t,
                                                                fnft__akns_discretization_t const akns_discretization,
                                                                FNFT_UINT vanilla_flag,
                                                                fnft__akns_pde_t const PDE);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define akns_discretization_degree(...) fnft__akns_discretization_degree(__VA_ARGS__)
#define akns_discretization_boundary_coeff(...) fnft__akns_discretization_boundary_coeff(__VA_ARGS__)
#define akns_discretization_upsampling_factor(...) fnft__akns_discretization_upsampling_factor(__VA_ARGS__)
#define akns_discretization_method_order(...) fnft__akns_discretization_method_order(__VA_ARGS__)
#define akns_discretization_lambda_to_z(...) fnft__akns_discretization_lambda_to_z(__VA_ARGS__)
#define akns_discretization_z_to_lambda(...) fnft__akns_discretization_z_to_lambda(__VA_ARGS__)
#define akns_discretization_method_weights(...) fnft__akns_discretization_method_weights(__VA_ARGS__)
#define akns_discretization_es6_stencil(...) fnft__akns_discretization_es6_stencil(__VA_ARGS__)
#define akns_discretization_preprocess_signal(...) fnft__akns_discretization_preprocess_signal(__VA_ARGS__)
#define akns_discretization_change_of_basis_matrix_to_S(...) fnft__akns_discretization_change_of_basis_matrix_to_S(__VA_ARGS__)
#define akns_discretization_change_of_basis_matrix_from_S(...) fnft__akns_discretization_change_of_basis_matrix_from_S(__VA_ARGS__)
#endif

#endif
