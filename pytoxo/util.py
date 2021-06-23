# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python tool to calculate penetrance tables for 
# high-order epistasis models
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

"""PyToxo util module."""

import functools
from threading import Thread


def timeout(timeout):
    """Timeout wrapper to associate a given timeout with a function. If the
    time is exceeded, a `TimeoutError` is raised. This approach works both in
    Unix-like and Windows machines. These last are more complicated because
    do not support signals well, which are the normally workaround to
    achieve these stuff.

    Raises
    ------
    TimeoutError
        If the configured timeout is exceeded.
    """

    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = [TimeoutError()]

            def target_helper():
                """Overrides the `res` variable, if `func` works."""
                res[0] = func(*args, **kwargs)

            t = Thread(target=target_helper)
            t.daemon = True
            try:
                t.start()
                t.join(timeout)
            except:
                """Possible but very unusual error handling multiprocessing.
                Ignore them seems enough workaround"""
                pass

            """Check if `res` has been overridden by `func`. If not, it signs 
            that `func` has been interrupted, so raise the time exceeded 
            error"""
            if isinstance(res[0], BaseException):
                raise TimeoutError
            return res[0]

        return wrapper

    return deco
