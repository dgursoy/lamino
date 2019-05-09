#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for external library wrappers.
"""

import sys
import os.path
import ctypes
import numpy as np
import logging
import six

logger = logging.getLogger(__name__)


__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2015, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'
__all__ = ['project', 'art', 'mlem']


LIB_LAMINO = ctypes.cdll.LoadLibrary('liblamino.cpython-36m-darwin.so')


def project(obj, gridx, gridy, gridz, phi, theta, detgridx, detgridy):
    prj = np.zeros((phi.size, detgridx.size, detgridy.size), dtype='float32')
    LIB_LAMINO.project.restype = as_c_void_p()
    LIB_LAMINO.project(
        as_c_float_p(obj),
        as_c_int(obj.shape[0]),
        as_c_int(obj.shape[1]),
        as_c_int(obj.shape[2]),
        as_c_float_p(gridx),
        as_c_float_p(gridy),
        as_c_float_p(gridz),
        as_c_float_p(phi),
        as_c_float_p(theta),
        as_c_float_p(detgridx),
        as_c_float_p(detgridy),
        as_c_int(detgridx.size),
        as_c_int(detgridy.size),
        as_c_float_p(prj),
        as_c_int(prj.shape[0]),
        as_c_int(prj.shape[1]),
        as_c_int(prj.shape[2]))
    return prj


def art(prj, phi, theta, detgridx, detgridy, gridx, gridy, gridz, iters):
    rec = np.zeros((gridx.size-1, gridy.size-1, gridz.size-1), dtype='float32')
    LIB_LAMINO.art.restype = as_c_void_p()
    LIB_LAMINO.art(
        as_c_float_p(prj),
        as_c_int(prj.shape[0]),
        as_c_int(prj.shape[1]),
        as_c_int(prj.shape[2]),
        as_c_float_p(phi),
        as_c_float_p(theta),
        as_c_float_p(detgridx),
        as_c_float_p(detgridy),
        as_c_int(detgridx.size),
        as_c_int(detgridy.size),
        as_c_float_p(rec),
        as_c_int(rec.shape[0]),
        as_c_int(rec.shape[1]),
        as_c_int(rec.shape[2]),
        as_c_float_p(gridx),
        as_c_float_p(gridy),
        as_c_float_p(gridz),
        as_c_int(iters))
    return rec


def mlem(prj, phi, theta, detgridx, detgridy, gridx, gridy, gridz, iters):
    rec = np.zeros((gridx.size-1, gridy.size-1, gridz.size-1), dtype='float32')
    LIB_LAMINO.art.restype = as_c_void_p()
    LIB_LAMINO.art(
        as_c_float_p(prj),
        as_c_int(prj.shape[0]),
        as_c_int(prj.shape[1]),
        as_c_int(prj.shape[2]),
        as_c_float_p(phi),
        as_c_float_p(theta),
        as_c_float_p(detgridx),
        as_c_float_p(detgridy),
        as_c_int(detgridx.size),
        as_c_int(detgridy.size),
        as_c_float_p(rec),
        as_c_int(rec.shape[0]),
        as_c_int(rec.shape[1]),
        as_c_int(rec.shape[2]),
        as_c_float_p(gridx),
        as_c_float_p(gridy),
        as_c_float_p(gridz),
        as_c_int(iters))
    return rec


def as_ndarray(arr, dtype=None, copy=False):
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=dtype, copy=copy)
    return arr


def as_dtype(arr, dtype, copy=False):
    if not arr.dtype == dtype:
        arr = np.array(arr, dtype=dtype, copy=copy)
    return arr


def as_float32(arr):
    arr = as_ndarray(arr, np.float32)
    return as_dtype(arr, np.float32)


def as_int32(arr):
    arr = as_ndarray(arr, np.int32)
    return as_dtype(arr, np.int32)


def as_uint16(arr):
    arr = as_ndarray(arr, np.uint16)
    return as_dtype(arr, np.uint16)


def as_uint8(arr):
    arr = as_ndarray(arr, np.uint8)
    return as_dtype(arr, np.uint8)


def as_c_float_p(arr):
    c_float_p = ctypes.POINTER(ctypes.c_float)
    return arr.ctypes.data_as(c_float_p)


def as_c_int(arr):
    return ctypes.c_int(arr)


def as_c_int_p(arr):
    c_int_p = ctypes.POINTER(ctypes.c_int)
    return arr.ctypes.data_as(c_int_p)


def as_c_float(arr):
    return ctypes.c_float(arr)


def as_c_char_p(arr):
    return ctypes.c_char_p(six.b(arr))


def as_c_void_p():
    return ctypes.POINTER(ctypes.c_void_p)