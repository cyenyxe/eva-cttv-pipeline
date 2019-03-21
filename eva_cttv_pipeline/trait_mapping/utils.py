import logging
logger = logging.getLogger(__package__)


def request_retry_helper(function, retry_count: int, url: str):
    """
    Given a function make a number of attempts to call function for it to successfully return a
    non-None value, subsequently returning this value. Makes the number of tries specified in
    retry_count parameter.

    :param function: Function that could need multiple attempts to return a non-None value
    :param retry_count: Number of attempts to make
    :param url: String specifying the url to make a request.
    :return: Returned value of the function.
    """
    for retry_num in range(retry_count):
        return_value = function(url)
        if return_value is not None:
            return return_value
        logger.warning("attempt {}: failed running function {} with url {}".format(
            retry_num, function, url))
    logger.warning("error on last attempt, skipping")
    return None
