import os
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
    '''
    NB: the old version used processes=len(data_iter), which raises
    TypeError for generators and forks len(data_iter) processes
    '''
    data = list(data_iter)
    if not data:
        return []
    workers = min(len(data), os.cpu_count() or 1)
    with multiprocessing.Pool(processes=workers) as pool:
        results = pool.map(func, data)
    return results