import signal
import functools

class TimeoutError(Exception):
    pass

def timeout(minutes):
    seconds = int(minutes * 60)
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            def timeout_handler(signum, frame):
                raise TimeoutError(f"{func.__name__} timed out after {minutes} minutes")
            
            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
            return result
        return wrapper
    return decorator