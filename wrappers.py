import numpy as np
import os
from PIL import Image
import traceback
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy import stats
from collections import defaultdict
from csv import writer
from skimage.segmentation import random_walker
from skimage.morphology import closing, dilation
from scipy import ndimage as ndi
from scipy.signal import medfilt2d
from skimage.morphology import watershed
from skimage.feature import peak_local_max
from skimage.morphology import disk
from skimage.morphology import skeletonize, medial_axis
from skimage.filters import threshold_otsu, rank, median
from skimage.draw import line_aa
from scipy.stats import t
from scipy.stats import ttest_ind
from itertools import combinations
from functools import wraps
import types
import collections
import debug_renders as dbg


dtype2bits = {'uint8': 8,
              'uint16': 16,
              'uint32': 32}


def safe_dir_create(path):
    if not os.path.isdir(path):
        os.makedirs(path)


class PipeArgError(ValueError):
    pass


def pad_skipping_iterator(secondary_namespace):
    for key, value in secondary_namespace.iteritems():
        if key != '_pad':
            yield value

def doublewrap(f):
    """
    a decorator decorator, allowing the decorator to be used as:
    @decorator(with, arguments, and=kwargs)
    or
    @decorator

    credits: http://stackoverflow.com/questions/653368/how-to-create-a-python-decorator-that-can-be-used-either-with-or-without-paramet

    """
    @wraps(f)
    def new_dec(*args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
            # actual decorated function
            return f(args[0])
        else:
            # decorator arguments
            return lambda realf: f(realf, *args, **kwargs)

    return new_dec

def list_not_string(argument):
    """
    function that checks if a list is not a string

    credits: http://stackoverflow.com/questions/1055360/how-to-tell-a-variable-is-iterable-but-not-a-string

    :param argument:
    :return:
    """
    if isinstance(argument, collections.Iterable):
        if isinstance(argument, types.StringTypes):
            return False
        else:
            return True
    else:
        raise PipeArgError("Expected a name of channel or list of names. Found: '%s' " % argument)
#

@doublewrap
def generator_wrapper(f, in_dims=(3,), out_dims=None):

    if out_dims is None:
            out_dims = in_dims

    @wraps(f)
    def inner_wrapper(*args, **kwargs):
        """
        converts a function to accepting a generator of named dicts and adds channel selection logic
        """

        iterator = args[0]
        args = args[1:]

        if 'in_channel' in kwargs:
            #Start in/out channel logic

            in_chan = kwargs['in_channel']  # Multiple arguments
            del kwargs['in_channel']

            if not list_not_string(in_chan):  # convert string to list of strings
                in_chan = [in_chan]

            if 'out_channel' in kwargs:
                out_chan = kwargs['out_channel']  # output explicitely provided
                del kwargs['out_channel']
                if not list_not_string(out_chan):  # convert string to list of strings
                    out_chan = [out_chan]

            else:  # implicit output, bound to in_channel only if a single input is provided
                if len(in_chan) == 1:
                    print 'Input %s will be overwritten by function %s' % (in_chan[0], f.__name__)
                    out_chan = in_chan
                else:
                    print f.__name__
                    print in_chan, in_dims
                    raise PipeArgError('Please provide out_channel argument')

            if len(in_chan) != len(in_dims):
                print f.__name__
                print in_chan, in_dims
                print len(in_chan), len(in_dims)
                raise PipeArgError('%s inbound channels are piped, function allows %s' %
                                   (len(in_chan), len(in_dims)))

            if len(out_chan) != len(out_dims):
                print f.__name__
                print out_chan, out_dims
                print len(out_chan), len(out_dims)
                raise PipeArgError('%s outbound channels are piped, function allows %s' %
                                   (len(out_chan), len(out_dims)))
            # end in/out channel logic

            for name_space in iterator:
                # start args prepare
                args_puck = []

                for i, chan in enumerate(in_chan):
                    if in_dims[i] and len(name_space[chan].shape) != in_dims[i]:
                        print f.__name__
                        print chan, len(name_space[chan].shape), in_dims[i]
                        raise PipeArgError('Mismatched inbound channel dimension for channel. %s is of dim %s, expected %s'%
                                           (chan, len(name_space[chan].shape), in_dims[i]))
                    # print name_space.keys()
                    args_puck.append(name_space[chan])

                local_args = tuple(args_puck) + args
                # end args prepare
                return_puck = f(*local_args, **kwargs)

                if return_puck is None and out_chan[0] == '_':
                    yield name_space  # unlike return, yield is probably non-blocking....

                else:
                    # start output prepare
                    if not isinstance(return_puck, tuple):
                        return_puck = (return_puck, )

                    for i, chan in enumerate(out_chan):
                        if out_dims[i] and len(return_puck[i].shape) != out_dims[i]:
                            print f.__name__
                            print chan
                            raise PipeArgError('Mismatched outgoing channel dimension for channel. %s is of dim %s, expected %s' %
                                               (chan, len(return_puck[i].shape), out_dims[i]))
                        if chan != '_':
                            name_space[chan] = return_puck[i]
                    # end output prepare

                    yield name_space

        else:
            for name_space in iterator:

                local_args = (name_space,) + args
                name_space = f(*local_args, **kwargs)
                yield name_space

    return inner_wrapper