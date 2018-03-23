#!/usr/bin/python env
import os
import pandas as pd

desired_width = 320
pd.set_option('display.width', desired_width)


def qc_location(df, flag_value=128):
    df['VLOC'] = 1
    boolean = df['VFLG'] == flag_value

    df['VLOC'] = df['VLOC'].where(~boolean, other=4)
    return df


def qc_radial_count(df, min_radials=150, low_radials=300):
    num_radials = df.shape[0]
    if num_radials < min_radials:
        radial_count_flag = 4
    elif (num_radials >= min_radials) and (num_radials <= low_radials):
        radial_count_flag = 3
    elif num_radials > low_radials:
        radial_count_flag = 1

    return radial_count_flag


def qc_speed(df, threshold=250):
    df['MVEL'] = 1
    try:
        boolean = df['VELO'].abs() > threshold
    except TypeError:
        df['VELO'] = df['VELO'].astype(float)
        boolean = df['VELO'].abs() > threshold

    df['MVEL'] = df['MVEL'].where(~boolean, other=4)
    return df


def qc_syntax(header_dict, df):
    return header_dict, df


# def qc_average_radial_bearing(df):
#     return df


# def qc_spatial_median_filter():
    # """
    # """
    # return def


# def qc_temporal_gradient(df, deriv_limit):
#     """
#     Remove data points that have a backward derivative greater than the deriv_limit that is passed to the function
#     :param df:
#     :param deriv_limit:
#     :return:
#     """
#     # function[DATA, nout, ind] = backwardsTemporalDerivative(rad_vel, deriv_limit)
#     # deriv_limit is expressed in velocity change(cm / s) over the course of an hour
#
#     # Convert the derivative limit from cm / s * hour to cm / s * s
#     deriv_limit = deriv_limit / 3600  # units cm / s ^ 2
#
#     # use the diff function to calculate the derivative
#     deriv = diff(rad_vel(:, 2))./ (24 * 60 * 60 * diff(rad_vel(:, 1)))
#
#     # the backward derivative is the output of the diff calculation with the last value removed
#     bckwd_deriv = abs(deriv);
#     bckwd_deriv(length(bckwd_deriv)) = []
#
#     # Find the indices of the backward derivatives that are greater than the derivative limit
#     ind = find(bckwd_deriv > deriv_limit)
#
#     # find the number of outliers
#     nout = length(ind);
#
#     # add 1 to the indices so that they match the indices of the vectors that you are screening
#     ind = ind + 1;
#
#     # remove the outliers
#     rad_vel(ind,:)=[];
#
#     # assign the new matrix to the variable DATA
#     DATA = rad_vel;
#     return df

if __name__ == '__main__':
    file_path = 'data/radials/AMAG/RDLi_AMAG_2017_05_03_1800.ruv'

    headers = 'LOND LATD VELU VELV VFLG ESPC ETMP MAXV MINV ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC'

    headers = headers.split(' ')
    radial = os.path.join(os.path.dirname(__file__), '..', file_path)

    df = pd.read_csv(radial, comment='%', skipinitialspace=True, sep=' ', names=headers, na_values=['999.000'])
