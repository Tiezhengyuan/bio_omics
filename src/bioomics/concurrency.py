import threading
import multiprocessing

# multiple threads
def multi_threads(data_iter, func):
    threads = []
    for _data in data_iter:
        thread = threading.Thread(target=func, args=(_data,))
        threads.append(thread)
        thread.start()
    
    for thread in threads:
        thread.join()


# multiple processess
def multi_processes(data_iter, func):
    processes = []
    for _data in data_iter:
        process = multiprocessing.Process(target=func, args=(_data,))
        process.start()
        processes.append(process)
    
    for p in processes:
        p.join()


def mp_pool(data_iter, func):
    pool = multiprocessing.Pool(processes=len(data_iter))

    results = pool.map(func, data_iter)
    return results