import time
import datetime


def loop_with_progress(iterator, n_iter, progress_step, printer_function, with_enumerate=False):
    """
    Generator which prints progress for loops

    :param iterator:
    :param n_iter:
    :param progress_step:
    :param printer_function:
    :return:
    """

    # Gather current time
    start_time = time.time()

    for i, obj in enumerate(iterator):

        if with_enumerate:

            yield i, obj

        else:

            yield obj

        if (i+1) % progress_step == 0:

            this_time = time.time()
            elapsed_time = (this_time - start_time)
            time_per_ts = elapsed_time / (i + 1)
            remaining_time = (n_iter - (i + 1)) * time_per_ts

            msg = "Processed %i out of %i\n" % (i + 1, n_iter)
            msg += "Elapsed time: %s, remaining time: %s" % (datetime.timedelta(seconds=elapsed_time),
                                                             datetime.timedelta(seconds=remaining_time))

            printer_function(msg)